#!/usr/bin/env python3
"""Extract the pair-distribution / structure-factor / diffraction block out of
turbogap.f90 into compute_exp_spectra in turbogap_exp.f90.

The block moves verbatim apart from one transformation: the four
`#ifdef _MPIF90` groups whose two branches differ only by the `this_` prefix on
an energies/forces/virial triple collapse to the plain name, because inside a
procedure the dummy has one name either way and the caller chooses once.

Every step asserts. Nothing is written unless all of them pass.
"""
import re
import sys

DRIVER = 'src/turbogap.f90'
EXPMOD = 'src/turbogap_exp.f90'
LO, HI = 1979, 2680          # inclusive, 1-based

ARGS = [
    # mirrors master's compute_exp_spectra order
    'params', 'n_sites', 'species', 'rjs', 'xyz', 'neighbors_list',
    'n_neigh', 'neighbor_species', 'indices', 'a_box', 'b_box', 'c_box',
    'i_beg', 'i_end', 'j_beg', 'j_end', 'rank', 'ntasks', 'ierr',
    'md_istep', 'mc_istep',
    'energies_sf', 'forces_sf', 'virial_sf',
    'energies_xrd', 'forces_xrd', 'virial_xrd',
    'energies_nd', 'forces_nd', 'virial_nd',
    'time_sf', 'time_xrd', 'time_nd', 'time_exp_batched',
    # GPU-only: the batched-device machinery
    'i_beg_list', 'i_end_list', 'j_beg_list', 'j_end_list',
    'n_omp', 'omp_task', 'this_i_beg', 'this_i_end', 'this_j_beg',
    'this_j_end', 'n_sites_temp', 'n_pairs_temp',
    'write_condition', 'overwrite_condition', 'temp_string',
    'species_types_actual', 'v_uc',
]

DUMMY_DECL = """
!   Read-only inputs.
    integer, intent(in) :: n_sites, i_beg, i_end, j_beg, j_end, rank, ntasks
    integer, intent(in) :: indices(1:3), md_istep, mc_istep
    real(dp), intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real(dp), intent(in), allocatable :: rjs(:), xyz(:,:)
    integer,  intent(in), allocatable :: neighbors_list(:), n_neigh(:), &
         neighbor_species(:), species(:)

!   params carries the do_* switches and the exp_data containers, which the
!   exp_interface routines update, so it cannot be intent(in).
    type(input_parameters), intent(inout) :: params
    integer, intent(inout) :: ierr

!   The three contribution families this routine produces. The caller passes
!   the this_-prefixed arrays under MPI and the plain ones otherwise; in here
!   there is one name either way, which is what removes four
!   preprocessor-interrupted argument lists from the driver.
    real(dp), allocatable, intent(inout) :: energies_sf(:),  forces_sf(:,:)
    real(dp), allocatable, intent(inout) :: energies_xrd(:), forces_xrd(:,:)
    real(dp), allocatable, intent(inout) :: energies_nd(:),  forces_nd(:,:)
    real(dp), intent(inout) :: virial_sf(1:3,1:3), virial_xrd(1:3,1:3), &
         virial_nd(1:3,1:3)

!   Timing buckets, accumulated across snapshots.
    real(dp), intent(inout) :: time_sf(1:3), time_xrd(1:3), time_nd(1:3), &
         time_exp_batched(1:3)

!   The batch decomposition. Built by the caller for the electrostatics and
!   SOAP paths and consumed here too, which is why it is passed rather than
!   held in gpu_context.
    integer, allocatable, intent(inout) :: i_beg_list(:), i_end_list(:), &
         j_beg_list(:), j_end_list(:)

!   Per-batch bounds and the OpenMP task index. omp_task is intended to be
!   thread-private, so it stays a variable the caller owns rather than moving
!   into gpu_context with the arrays it indexes.
    integer, intent(inout) :: n_omp, omp_task, this_i_beg, this_i_end, &
         this_j_beg, this_j_end, n_sites_temp, n_pairs_temp

!   Output bookkeeping shared with the XPS path.
    logical, intent(inout) :: write_condition, overwrite_condition
    character*1024, intent(inout) :: temp_string
    character*8, allocatable, intent(inout) :: species_types_actual(:)
    real(dp), intent(inout) :: v_uc
"""

LOCAL_DECL = """
!   Block-private. All of these were driver declarations that nothing outside
!   this block referenced.
    real(dp), parameter :: pi = acos(-1.0)

!   i, j, k and f were driver scratch. Each is written before it is read here,
!   and the driver rewrites i and j (do-loops) before next reading them; k and f
!   it never reads again. f is only ever used inside this block.
    integer :: i, j, k
    real(dp) :: f
    integer :: n_dim_idx, n_dim_partial, n_species_actual, q_beg, q_end

    real(dp), allocatable, target :: dV(:), n_atoms_of_species(:), prefactor(:)
    real(dp), allocatable, target :: sinc_factor_matrix(:,:)
    real(dp), allocatable, target :: pair_distribution_der(:,:), &
         pair_distribution_partial(:,:), pair_distribution_partial_temp(:,:), &
         pair_distribution_partial_der(:,:,:), &
         pair_distribution_partial_temp_der(:,:,:)
    real(dp), allocatable, target :: structure_factor_partial(:,:), &
         structure_factor_partial_temp(:,:)
    real(dp), allocatable, target :: x_pair_distribution(:), &
         y_pair_distribution(:), y_pair_distribution_temp(:), &
         x_structure_factor(:), x_structure_factor_temp(:), &
         y_structure_factor(:), y_structure_factor_temp(:), &
         x_xrd(:), x_xrd_temp(:), y_xrd(:), y_xrd_temp(:), &
         x_nd(:), x_nd_temp(:), y_nd(:), y_nd_temp(:)

    integer, allocatable :: nk(:)

!   Device-side scratch, private to this block.
    type(c_ptr) :: xpdf_d, dV_d, prefactor_d, sinc_factor_matrix_d
    type(c_ptr), allocatable :: nk_d(:), k_index_d(:), j2_index_d(:), &
         xyz_k_d(:), pair_distribution_partial_d(:), &
         pair_distribution_partial_der_d(:), all_scattering_factors_d(:)
    integer(c_size_t) :: st_x_d, st_sinc_factor_matrix_d, st_prefactor_d
    integer(c_size_t), allocatable :: st_nk_d(:), st_k_index_d(:), &
         st_j2_index_d(:), st_pair_distribution_partial_d(:), &
         st_pair_distribution_partial_der_d(:)

    type( gpu_host_storage_type ), allocatable :: gpu_host_exp_storage(:)
"""


def collapse_ifdefs(body):
    """Replace each #ifdef/#else/#endif group with its #else branch.

    Asserts that every group has an #else and that stripping `this_` from the
    #ifdef branch reproduces the #else branch, so a group that guards anything
    other than the this_/plain choice cannot be silently swallowed.
    """
    out, i, collapsed, kept = [], 0, 0, 0
    while i < len(body):
        s = body[i].strip()
        if s.startswith('#ifdef') or s.startswith('#if '):
            depth, j, else_at = 1, i + 1, None
            while j < len(body) and depth:
                t = body[j].strip()
                if t.startswith('#ifdef') or t.startswith('#if '):
                    depth += 1
                elif t.startswith('#else') and depth == 1:
                    else_at = j
                elif t.startswith('#endif'):
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            assert depth == 0, f'unterminated #ifdef at body line {i+1}'
            if else_at is None:
                # No #else: an ordinary conditional over whole statements.
                # It must survive the move untouched.
                out.extend(body[i:j + 1])
                kept += 1
            else:
                a = body[i + 1:else_at]
                b = body[else_at + 1:j]

                def norm(ls):
                    t = ' '.join(x.split('!', 1)[0].strip() for x in ls)
                    return re.sub(r'\s+', '', t.replace('&', ' ')).lower()

                stripped = re.sub(r'this_(energies|forces|virial)_', r'\1_', norm(a))
                assert stripped == norm(b), (
                    f'body line {i+1}: #ifdef branch is not the #else branch '
                    f'modulo the this_ prefix -- refusing to collapse')
                out.extend(b)
                collapsed += 1
            i = j + 1
            continue
        out.append(body[i])
        i += 1
    return out, collapsed, kept


def main():
    lines = open(DRIVER, encoding='utf-8').read().splitlines()
    body = lines[LO - 1:HI]
    assert 'Pair distribution functions and XRD' in body[1], \
        f'block does not start at the expected banner: {body[1]!r}'
    assert 'deallocate(species_types_actual)' in body[-1], \
        f'block does not end at the expected statement: {body[-1]!r}'

    new_body, collapsed, kept = collapse_ifdefs(body)
    print(f'#ifdef groups: collapsed {collapsed}, kept {kept}')
    assert collapsed == 4 and kept == 0, 'unexpected directive tally'
    assert not any(l.strip().startswith('#') for l in new_body), \
        'preprocessor directives remain in the moved body'

    # Three references to the this_ arrays sit OUTSIDE any #ifdef: the block
    # allocates and zeroes this_forces_xrd/this_virial_xrd unconditionally, then
    # hands either those or the plain pair to collect_batched_forces. Since the
    # dummy has one name, these must follow it. Both trees build with
    # `-D _MPIF90` unconditionally (see makefiles/Makefile.*), so the #ifdef
    # branch is the one that compiles and this reproduces it exactly.
    renamed = 0
    for n, line in enumerate(new_body):
        if re.search(r'\bthis_(energies|forces|virial)_(pdf|sf|xrd|nd)\b', line):
            new_body[n] = re.sub(r'\bthis_(energies|forces|virial)_(pdf|sf|xrd|nd)\b',
                                 r'\1_\2', line)
            renamed += 1
    print(f'bare this_ references renamed to the dummy name: {renamed}')
    assert renamed == 3, f'expected 3 bare this_ lines, found {renamed}'
    assert not re.search(r'\bthis_(energies|forces|virial)_', '\n'.join(new_body)), \
        'this_ contribution names remain in the moved body'

    # --- build the subroutine ---
    sig = 'compute_exp_spectra( ' + ', '.join(ARGS) + ' )'
    wrapped, cur = [], '  subroutine ' + sig
    while len(cur) > 96:
        cut = cur.rfind(', ', 0, 92) + 1
        wrapped.append(cur[:cut] + ' &')
        cur = '       ' + cur[cut:].lstrip()
    wrapped.append(cur)

    sub = []
    sub.append('!**************************************************************************')
    sub.append('! The (partial) pair distribution functions, structure factors and')
    sub.append('! diffraction patterns, and the energies and forces that fitting them')
    sub.append('! against experiment produces.')
    sub.append('!')
    sub.append(f'! Lifted from turbogap.f90 (lines {LO}-{HI} at the time of the move). The')
    sub.append('! only change to the body is that four argument lists split across')
    sub.append('! #ifdef _MPIF90 collapsed: the two branches differed only by the this_')
    sub.append('! prefix on an energies/forces/virial triple, and in here the dummy has')
    sub.append('! one name either way, so the caller chooses once at the single call.')
    sub.append('!')
    sub.append('! The device context this block shares with the electrostatics and SOAP')
    sub.append('! paths -- streams, cuBLAS handles, gpu_exp, gpu_neigh, gpu_batch_storage')
    sub.append('! -- comes from gpu_context by USE, which is what lets the body move')
    sub.append('! unchanged instead of growing eight more arguments.')
    sub.append('!**************************************************************************')
    sub.extend(wrapped)
    sub.append('')
    sub.append('    implicit none')
    sub.append('')
    sub.append('!   ---- dummy arguments ----')
    sub.extend(DUMMY_DECL.strip('\n').splitlines())
    sub.append('')
    sub.append('!   ---- local variables ----')
    sub.extend(LOCAL_DECL.strip('\n').splitlines())
    sub.append('')
    sub.extend(new_body)
    sub.append('')
    sub.append('  end subroutine compute_exp_spectra')
    sub.append('!**************************************************************************')

    exp = open(EXPMOD, encoding='utf-8').read().splitlines()
    end_idx = [i for i, l in enumerate(exp) if l.strip().lower().startswith('end module')]
    assert len(end_idx) == 1, 'could not find a unique `end module`'
    exp = exp[:end_idx[0]] + [''] + sub + [''] + exp[end_idx[0]:]

    pub = [i for i, l in enumerate(exp) if l.strip() == 'public :: compute_exp_xps']
    assert len(pub) == 1, 'could not find the public statement'
    exp[pub[0]] = '  public :: compute_exp_xps, compute_exp_spectra'

    for mod, after in (('use gpu_context', '  use exp_interface'),
                       ('use F_B_C', '  use exp_interface'),
                       ('use iso_c_binding', '  use exp_interface')):
        if not any(l.strip().lower() == mod.lower() for l in exp):
            at = [i for i, l in enumerate(exp) if l.rstrip() == after]
            assert len(at) == 1, f'anchor {after!r} not unique'
            exp.insert(at[0] + 1, '  ' + mod)

    open(EXPMOD, 'w', encoding='utf-8').write('\n'.join(exp) + '\n')

    # --- patch the driver ---
    def call_lines(prefix):
        actual = [f'{prefix}{a}' if a in THIS_ABLE else a for a in ARGS]
        text, res, cur = 'call compute_exp_spectra( ' + ', '.join(actual) + ' )', [], ''
        cur = '         ' + text
        while len(cur) > 96:
            cut = cur.rfind(', ', 0, 92) + 1
            res.append(cur[:cut] + ' &')
            cur = '              ' + cur[cut:].lstrip()
        res.append(cur)
        return res

    THIS_ABLE = {'energies_sf', 'forces_sf', 'virial_sf',
                 'energies_xrd', 'forces_xrd', 'virial_xrd',
                 'energies_nd', 'forces_nd', 'virial_nd'}

    repl = []
    repl.append('         !##############################################################!')
    repl.append('         !###---   (Partial) Pair distribution functions and XRD   ---###!')
    repl.append('         !##############################################################!')
    repl.append('         !')
    repl.append('         ! Moved to src/turbogap_exp.f90. The #ifdef is here, at the one')
    repl.append('         ! call, rather than inside four continued argument lists where')
    repl.append('         ! nothing Fortran-aware could parse it.')
    repl.append('#ifdef _MPIF90')
    repl.extend(call_lines('this_'))
    repl.append('#else')
    repl.extend(call_lines(''))
    repl.append('#endif')
    repl.append('')

    lines = lines[:LO - 1] + repl + lines[HI:]
    open(DRIVER, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

    print(f'driver: replaced {HI-LO+1} lines with {len(repl)}')
    print(f'turbogap_exp.f90: +{len(sub)} lines, {len(ARGS)} arguments')


if __name__ == '__main__':
    sys.exit(main())
