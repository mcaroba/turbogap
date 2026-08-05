#!/usr/bin/env python3
"""Split read_input_file's keyword chain into eleven themed subroutines.

Groundwork this relies on, all measured rather than assumed:

  * The chain's own `else if`/`else`/`end if` sit at one indent. Handlers contain
    nested `else` and `end do` at deeper indents, so every boundary test here is
    indent-anchored -- a naive scan finds 59 handlers instead of 229.
  * Of read_input_file's locals, `bsf j nw iostatus2 long_line long_line_items`
    are used only inside the chain, and `i i2 k keyword_notrim valid_choice`
    cross but carry nothing: `i` is rewritten by `do i = 1, n_species` after the
    chain, three are never read again, and `valid_choice` is written by
    get_atomic_mass (intent(out)) before its next read. All become locals.
  * `are_vdw_refs_read` and `masses_in_input_file` genuinely carry a value out of
    the chain, so they are passed to the two subroutines that set them.
  * implemented_thermostats/barostats/mc_types are constant lists, so they move
    to module parameters instead of being rebuilt per call.
"""
import re
import sys
from collections import OrderedDict

PATH = 'src/read_files.f90'

EXACT = {
    'atoms_file': 'general', 'species': 'general', 'masses': 'general',
    'radii': 'general', 'e0': 'general', 'which_atom': 'general',
    'max_gbytes_per_process': 'general', 'neighbors_buffer': 'general',
    'random_seed': 'general', 'timing': 'general', 'verbosity': 'general',
    'do_md': 'control', 'do_mc': 'control', 'do_prediction': 'control',
    'do_forces': 'control', 'do_derivatives': 'control',
    'do_derivatives_fd': 'control', 'write_soap': 'control',
    'write_derivatives': 'control', 'optimize': 'control', 'gamma0': 'control',
    'max_opt_step': 'control', 'max_opt_step_eps': 'control',
    'e_tol': 'control', 'f_tol': 'control', 'p_tol': 'control',
    'target_pos_step': 'control', 'tau_dt': 'control',
    't_beg': 'md', 't_end': 'md', 'p_beg': 'md', 'p_end': 'md',
    'tau_t': 'md', 'tau_p': 'md', 'md_step': 'md', 'md_nsteps': 'md',
    'thermostat': 'md', 'barostat': 'md', 'barostat_sym': 'md',
    'randomize_velocities': 'md', 'n_t_hold': 'md', 't_hold': 'md',
    'gamma_p': 'md', 'scale_box': 'md', 'box_scaling_factor': 'md',
    'n_nested': 'nested', 't_extra': 'nested', 'p_nested': 'nested',
    'nested_max_strain': 'nested', 'nested_max_volume_change': 'nested',
    'scale_box_nested': 'nested',
    'accessible_volume': 'mc',
    'poly_cut_xmin': 'vdw', 'poly_cut_xmax': 'vdw',
    'mbd_correction_freq': 'vdw', 'do_nnls': 'vdw', 'print_vdw_forces': 'vdw',
    'n_local_properties': 'local_properties',
    'compute_local_properties': 'local_properties',
    'core_pot_cutoff': 'local_properties', 'core_pot_buffer': 'local_properties',
    'print_lp_forces': 'local_properties',
    'adaptive_time': 'stopping', 'electronic_stopping': 'stopping',
    'nonadiabatic_processes': 'stopping', 'estop_filename': 'stopping',
    'model_eph': 'stopping',
    'write_xyz': 'output', 'write_thermo': 'output', 'write_lv': 'output',
    'print_progress': 'output',
    'print_estat_forces': 'estat',
    # GPU branch only
    'n_batches': 'gpu',
}
PREFIX = [('mc_', 'mc'), ('n_mc_', 'mc'), ('vdw_', 'vdw'), ('estat_', 'estat'),
          ('gpu_', 'gpu'),
          ('eph_', 'stopping'), ('eel_', 'stopping'), ('adapt_', 'stopping'),
          ('xps_', 'exp'), ('xrd_', 'exp'), ('nd_', 'exp'), ('sf_', 'exp'),
          ('pair_distribution_', 'exp'), ('structure_factor_', 'exp'),
          ('exp_', 'exp'), ('q_', 'exp'), ('r_range_', 'exp'), ('write_', 'output')]
EXP_EXTRA = {'do_pair_distribution', 'do_structure_factor', 'do_xrd', 'do_nd',
             'do_exp', 'n_exp', 'do_xps', 'structure_factor_window',
             'structure_factor_matrix', 'structure_factor_from_pdf'}
ORDER = ['general', 'control', 'md', 'nested', 'mc', 'vdw', 'estat', 'exp',
         'stopping', 'local_properties', 'output', 'gpu']
TITLE = {
    'general': 'the system: files, species, masses, and global limits',
    'control': 'what the run does: prediction, forces, optimisation',
    'md': 'molecular dynamics: thermostat, barostat, step and duration',
    'nested': 'nested sampling',
    'mc': 'Monte Carlo: move types, limits and acceptance',
    'vdw': 'van der Waals: TS, MBD, and the Hirshfeld model',
    'estat': 'electrostatics',
    'exp': 'experimental observables: pdf, structure factor, XRD, ND, XPS',
    'stopping': 'electronic stopping and the electron-phonon model',
    'local_properties': 'local properties and the core potential',
    'output': 'what gets written, and how often',
    'gpu': 'GPU execution: batching and device memory (GPU branch only)',
}
# locals a group may need, declared only when its handlers actually mention them
OPTIONAL_LOCALS = OrderedDict([
    ('i', '      integer :: i\n'),
    ('j', '      integer :: j\n'),
    ('k', '      integer :: k\n'),
    ('i2', '      integer :: i2\n'),
    ('nw', '      integer :: nw\n'),
    ('iostatus2', '      integer :: iostatus2\n'),
    ('bsf', '      real(dp) :: bsf\n'),
    ('long_line', '      character*1024 :: long_line\n'),
    ('long_line_items', '      character*128, allocatable :: long_line_items(:)\n'),
    ('keyword_notrim', '      character*64 :: keyword_notrim\n'),
    ('valid_choice', '      logical :: valid_choice\n'),
])
EXTRA_ARG = {'general': ('masses_in_input_file', '      logical, intent(inout) :: masses_in_input_file\n'),
             'vdw': ('are_vdw_refs_read', '      logical, intent(inout) :: are_vdw_refs_read(1:3)\n')}


def group_of(kws):
    for k in kws:
        if k in EXACT:
            return EXACT[k]
    for k in kws:
        if k in EXP_EXTRA:
            return 'exp'
    for k in kws:
        for pre, g in PREFIX:
            if k.startswith(pre):
                return g
    return None


src = open(PATH, encoding='utf-8').read()
if 'read_options_general' in src:
    print('ABORT: already split', file=sys.stderr)
    sys.exit(1)

start = src.index('subroutine read_input_file(')
mstop = re.search(r'^ *end subroutine\b', src[start + 1:], re.M)
sub_end = start + 1 + mstop.start()
head, body, tail = src[:start], src[start:sub_end], src[sub_end:]
lines = body.splitlines(keepends=True)

C = re.compile(r"^( *)if \(keyword_first == '#'")
ci = next(i for i, l in enumerate(lines) if C.match(l))
IND = C.match(lines[ci]).group(1)
KW = re.compile(r"^( *)(else if|if) \(keyword == ['\"]([a-z0-9_]+)['\"](.*)\) then\s*$")
ALT = re.compile(r"""keyword == ['"]([a-z0-9_]+)['"]""")

handlers, cur, chain_end = [], None, None
for i in range(ci + 1, len(lines)):
    l = lines[i]
    m = KW.match(l)
    if m and m.group(1) == IND:
        if cur:
            handlers.append((cur[0], i, cur[1]))
        cur = (i, ALT.findall(l))
    elif l.rstrip('\n') == IND + 'else' and cur:
        handlers.append((cur[0], i, cur[1]))
        chain_end = i
        break
if chain_end is None:
    print('ABORT: chain end not found', file=sys.stderr)
    sys.exit(1)
after_chain = next(i for i in range(chain_end, len(lines))
                   if lines[i].rstrip('\n') == IND + 'end if') + 1
print(f'parsed {len(handlers)} handlers')

groups = OrderedDict((g, []) for g in ORDER)
for s, e, kws in handlers:
    g = group_of(kws)
    if g is None:
        print(f'ABORT: unassigned {kws}', file=sys.stderr)
        sys.exit(1)
    groups[g].append((s, e, kws))


def rewrite(chunk):
    out = []
    for l in chunk:
        l = l.replace('backspace (10)', 'backspace (unit)')
        l = l.replace('read (10, *,', 'read (unit, *,')
        l = l.replace('read_parameters(10,', 'read_parameters(unit,')
        out.append(l)
    return out


ACTIVE = [g for g in ORDER if groups[g]]

subs = []
for g in ACTIVE:
    hs = groups[g]
    text = ''.join(''.join(rewrite(lines[s:e])) for s, e, _ in hs)
    extra = EXTRA_ARG.get(g)
    args = 'unit, iostatus, rank, keyword, mode, params, n_species, keyword_found'
    if extra:
        args += ', ' + extra[0]
    b = ['!**************************************************************************\n',
         f'!  Input keywords for {TITLE[g]}.\n', '!\n',
         '!  Sets keyword_found when it recognises the keyword and leaves it alone\n',
         '!  otherwise, so read_input_file can offer the line to the next family.\n',
         f'   subroutine read_options_{g}({args})\n', '\n',
         '      implicit none\n', '\n',
         '      integer, intent(in) :: unit, rank, n_species\n',
         '      integer, intent(inout) :: iostatus\n',
         '      character*64, intent(in) :: keyword\n',
         '      character(len=*), intent(in) :: mode\n',
         '      type(input_parameters), intent(inout) :: params\n',
         '      logical, intent(inout) :: keyword_found\n']
    if extra:
        b.append(extra[1])
    b.append('\n      character*64 :: cjunk\n')
    for name, decl in OPTIONAL_LOCALS.items():
        if re.search(r'\b' + name + r'\b', text):
            b.append(decl)
    b.append('\n')
    first = True
    for s, e, _ in hs:
        chunk = rewrite(lines[s:e])
        h = chunk[0]
        h = re.sub(r'^( *)else if ', r'\1if ', h) if first else \
            (h if h.lstrip().startswith('else if') else re.sub(r'^( *)if ', r'\1else if ', h, count=1))
        first = False
        b.append(h)
        b.extend(chunk[1:])
    b += ['      else\n', '         return\n', '      end if\n', '\n',
          '      keyword_found = .true.\n', '\n',
          f'   end subroutine read_options_{g}\n', '\n']
    subs.append(''.join(b))

d = [f"{IND}if (keyword_first == '#' .or. keyword_first == '!' .or. keyword == 'pot_file' &\n",
     f"{IND}    .or. keyword == 'n_species') cycle\n", '\n',
     f'{IND}! Offer the line to each family in turn. The first that recognises the\n',
     f'{IND}! keyword sets keyword_found; falling off the end means nothing claimed it.\n',
     f'{IND}keyword_found = .false.\n']
for g in ACTIVE:
    extra = EXTRA_ARG.get(g)
    a = '10, iostatus, rank, keyword, mode, params, n_species, keyword_found'
    if extra:
        a += ', ' + extra[0]
    d.append(f'{IND}call read_options_{g}({a})\n')
    d.append(f'{IND}if (keyword_found) cycle\n')
d += ['\n',
      f'{IND}write (*, *) "ERROR: I do not recognize the input file keyword ", trim(keyword)\n',
      f'{IND}call turbogap_abort()\n']

new_body = ''.join(lines[:ci]) + ''.join(d) + ''.join(lines[after_chain:])

# keyword_found declaration
m = re.search(r'^ *logical :: masses_in_input_file[^\n]*\n', new_body, re.M)
new_body = new_body[:m.end()] + '      logical :: keyword_found\n' + new_body[m.end():]

# the constant lists become module parameters
for name, decl in [('implemented_thermostats', r'^ *character\*32 :: implemented_thermostats[^\n]*\n'),
                   ('implemented_barostats', r'^ *character\*32 :: implemented_barostats[^\n]*\n'),
                   ('implemented_mc_types', r'^ *character\*32 :: implemented_mc_types[^\n]*\n')]:
    new_body = re.sub(decl, '', new_body, count=1, flags=re.M)
    new_body = re.sub(r'^ *' + name + r'\(\d+\) = [^\n]*\n', '', new_body, flags=re.M)

MODCONST = '''
!  Constant lists of what the code implements. Module parameters rather than
!  locals rebuilt on every call, so the keyword-family subroutines can validate
!  against them without being handed the lists.
   character*32, parameter :: implemented_thermostats(1:3) = &
      [character*32 :: "none", "berendsen", "bussi"]
   character*32, parameter :: implemented_barostats(1:2) = &
      [character*32 :: "none", "berendsen"]
   character*32, parameter :: implemented_mc_types(1:8) = &
      [character*32 :: "none", "move", "insertion", "removal", "relax", "md", &
       "swap", "volume"]

'''
# Anchor on `contains`, not on `implicit none`: this module has no module-level
# implicit none, so the first one belongs to a subroutine and the parameters
# would land inside its declaration section.
# The dispatch calls turbogap_abort on an unrecognised keyword.
if 'use error, only: turbogap_abort' not in head:
    mu = re.search(r'^ *use read_utils[^\n]*\n', head, re.M)
    if not mu:
        print('ABORT: no read_utils import to anchor on', file=sys.stderr)
        sys.exit(1)
    head = head[:mu.end()] + '   use error, only: turbogap_abort\n' + head[mu.end():]

mmod = re.search(r'^contains\s*$', head, re.M)
if not mmod:
    print('ABORT: no module-level contains', file=sys.stderr)
    sys.exit(1)
head = head[:mmod.start()] + MODCONST.lstrip('\n') + '\n' + head[mmod.start():]

nl = tail.index('\n') + 1
open(PATH, 'w', encoding='utf-8').write(head + new_body + tail[:nl] + '\n' + ''.join(subs) + tail[nl:])
print(f'wrote {len(subs)} subroutines; chain -> {len(ACTIVE)}-call dispatch')
for g in ACTIVE:
    print(f'   read_options_{g:18s} {len(groups[g]):3d} handlers')
