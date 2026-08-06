#!/usr/bin/env python3
"""Move the MD driver block out of turbogap.f90 into turbogap_md.f90.

The block runs from the "Do MD stuff here" banner to the end of the MPI
position broadcast, and includes its own rank guard, so the whole thing moves
and the driver is left with one call.

i, i2, j, j2 and k2 become locals: each is written before it is read in the
moved code, the driver rewrites i and j in do-loops afterwards, and i2, j2 and
k2 are never read again -- the leaked-value check from refactor_handoff.md 8.

cum_eel, gd_box_do_pos, gd_istep, restart_box_optim, target_temp and
time_step_prev are arguments even though nothing outside the block touches
them. The block runs inside the MD loop, so a driver local that is private to
it is still *per-run* state that persists between steps; a subroutine local
would reset on every call. Each of these is read before it is written within a
single step, which is what identifies them.
"""
import re
import sys

DRIVER = 'src/turbogap.f90'
MODULE = 'src/turbogap_md.f90'
# Boundaries discovered rather than hardcoded: the block sits at different line
# numbers on this branch.
_L = open(DRIVER, encoding='utf-8').read().splitlines()
_s = next(i for i, l in enumerate(_L) if 'Do MD stuff here' in l) - 1
_e = next(i for i in range(_s, len(_L))
          if _L[i].rstrip() == '#endif' and 'Nested sampling' in '\n'.join(_L[i:i + 8]))
LO, HI = _s + 1, _e + 1

ARGS = """params rank ierr n_sites md_istep md_time time_step
positions positions_prev positions_diff velocities forces forces_prev masses
xyz xyz_species a_box b_box c_box indices v_uc virial
energy energy_prev energies energies_soap energies_2b energies_3b energies_core_pot
energies_estat energies_vdw energies_lp energies_exp energies_pdf energies_sf energies_xrd
energies_nd local_properties local_property_labels
instant_temp instant_pressure instant_pressure_prev e_kin e_kinetic kb evpera3tobar
fix_atom exit_loop rebuild_neighbors_list i_image i_nested n_pos string
time_md time_mpi_positions
gd_box_do_pos gd_istep restart_box_optim target_temp time_step_prev""".split()

# filename and string are used as internal files (write (filename, ...)), which
# the standard forbids for intent(in), so they stay inout.
IN = set("""e_kin energies energy energy_prev evpera3tobar fix_atom i_image i_nested
instant_pressure_prev kb local_properties local_property_labels masses rank v_uc virial
xyz xyz_species""".split())

DECL = {
    'a_box': 'real(dp)', 'b_box': 'real(dp)', 'c_box': 'real(dp)',
    'allelstopdata': 'real(dp), allocatable', 'e_kin': 'real(dp)',
    'e_kinetic': 'real(dp)', 'energies': 'real(dp), allocatable',
    'energy': 'real(dp)', 'energy_prev': 'real(dp)',
    'ephbeta': 'type(EPH_Beta_class)', 'ephfdm': 'type(EPH_FDM_class)',
    'ephlsc': 'type(EPH_LangevinSpatialCorrelation_class)',
    'evpera3tobar': 'real(dp)', 'exit_loop': 'logical',
    'filename': 'character*1024', 'fix_atom': 'logical, allocatable',
    'forces': 'real(dp), allocatable', 'forces_prev': 'real(dp), allocatable',
    'i_image': 'integer', 'i_nested': 'integer', 'ierr': 'integer',
    'indices': 'integer', 'instant_pressure': 'real(dp)',
    'instant_pressure_prev': 'real(dp)', 'instant_temp': 'real(dp)',
    'kb': 'real(dp)', 'local_properties': 'real(dp), allocatable',
    'local_property_labels': 'character*1024, allocatable',
    'masses': 'real(dp), allocatable', 'masses_types': 'real(dp), allocatable',
    'md_istep': 'integer', 'md_time': 'real(dp)', 'n_pos': 'integer',
    'n_sites': 'integer', 'n_species': 'integer', 'nrows': 'integer',
    'params': 'type(input_parameters)', 'positions': 'real(dp), allocatable',
    'positions_diff': 'real(dp), allocatable',
    'positions_prev': 'real(dp), allocatable', 'rank': 'integer',
    'rebuild_neighbors_list': 'logical', 'string': 'character*1024',
    'time_md': 'real(dp)', 'time_mpi_positions': 'real(dp)',
    'time_step': 'real(dp)', 'v_uc': 'real(dp)',
    'velocities': 'real(dp), allocatable', 'virial': 'real(dp)',
    'xyz': 'real(dp), allocatable', 'xyz_species': 'character*8, allocatable',
    'cum_eel': 'real(dp)', 'gd_box_do_pos': 'logical', 'gd_istep': 'integer',
    'restart_box_optim': 'logical', 'target_temp': 'real(dp)',
    'time_step_prev': 'real(dp)',
}
DIMS = {'a_box': '(1:3)', 'b_box': '(1:3)', 'c_box': '(1:3)', 'indices': '(1:3)',
        'virial': '(1:3, 1:3)', 'time_md': '(1:3)', 'time_mpi_positions': '(1:3)',
        'allelstopdata': '(:)', 'energies': '(:)', 'fix_atom': '(:, :)',
        'forces': '(:, :)', 'forces_prev': '(:, :)', 'local_properties': '(:, :)',
        'local_property_labels': '(:)', 'masses': '(:)', 'masses_types': '(:)',
        'positions': '(:, :)', 'positions_diff': '(:, :)', 'positions_prev': '(:, :)',
        'velocities': '(:, :)', 'xyz': '(:, :)', 'xyz_species': '(:)'}
for e in ['soap', '2b', '3b', 'core_pot', 'vdw', 'lp', 'exp', 'pdf', 'sf', 'xrd', 'nd',
          'estat']:
    DECL[f'energies_{e}'] = 'real(dp), allocatable'
    DIMS[f'energies_{e}'] = '(:)'

LOCALS = """      character*64 :: cjunk
      real(dp) :: instant_pressure_tensor(1:3, 1:3)
      real(dp) :: lv(1:3, 1:3)
      character*1024 :: filename
!     Loop scratch. Written before read here; the driver rewrites i and j in
!     do-loops after the call and never reads i2, j2 or k2 again.
      integer :: i, i2, j, j2, k2
"""

src = open(DRIVER, encoding='utf-8').read()
if 'compute_md' in src:
    print('ABORT: already extracted', file=sys.stderr)
    sys.exit(1)
lines = src.splitlines(keepends=True)
body = lines[LO - 1:HI]
assert 'Do MD stuff here' in body[1], body[1]
assert 'endif' in body[-1].lower() or '#endif' in body[-1], repr(body[-1])

def wrap(prefix, items, width=92):
    out, cur = [], prefix
    for n, it in enumerate(items):
        piece = it + (', ' if n < len(items) - 1 else '')
        if len(cur) + len(piece) > width:
            out.append(cur.rstrip() + ' &')
            cur = ' ' * 8 + piece
        else:
            cur += piece
    out.append(cur)
    return out

sig = wrap('   subroutine compute_md(', ARGS)
sig[-1] += ')'

decls = []
for a in ARGS:
    t = DECL[a]
    intent = 'intent(in)' if a in IN else 'intent(inout)'
    if a == 'params':
        intent = 'intent(inout)'
    decls.append(f'      {t}, {intent} :: {a}{DIMS.get(a, "")}\n')

mod = []
mod.append('''! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_md.f90, is copyright (c) 2019-2026, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!  One MD step: the thermostat and barostat, the velocity-Verlet update, the
!  box scaling, the electronic-stopping and electron-phonon coupling, the
!  trajectory and thermo output, and the broadcast of the new positions.
!
!  Lifted verbatim from turbogap.f90. The block carried its own
!  `#ifdef _MPIF90 / IF (rank == 0)` guard and the position broadcast that
!  follows it, so both moved too and the driver is left with a single call.
module turbogap_md

   use kinds
   use types
   use md
   use bussi
   use xyz_module
   use timing
   use mpi_helper
!  No adaptive_time / electronic_stopping / eph_* here: this branch compiles
!  SRC_STOP but never links OBJ_STOP.
#ifdef _MPIF90
   use mpi
#endif

   implicit none

   private
   public :: compute_md

contains

!**************************************************************************
''')
mod.extend(l + '\n' for l in sig)
mod.append('\n      implicit none\n\n')
mod.extend(decls)
mod.append('\n')
mod.append(LOCALS)
mod.append('\n')
mod.extend(body)
mod.append('\n   end subroutine compute_md\n')
mod.append('!**************************************************************************\n')
mod.append('\nend module turbogap_md\n')
open(MODULE, 'w', encoding='utf-8').write(''.join(mod))

call = wrap('      call compute_md(', ARGS)
call[-1] += ')'
repl = ['      !**************************************************************************\n',
        '      !   Do MD stuff here. Moved to src/turbogap_md.f90; the rank guard and the\n',
        '      !   position broadcast moved with it.\n']
repl += [l + '\n' for l in call]
repl.append('\n')

lines = lines[:LO - 1] + repl + lines[HI:]

# The driver must USE the new module, or the call resolves as an external
# procedure and the link fails on compute_md_.
drv = ''.join(lines)
if 'use turbogap_md' not in drv:
    mu = re.search(r'^ *use turbogap_exp\s*$', drv, re.M)
    if not mu:
        print('ABORT: no `use turbogap_exp` to anchor on', file=sys.stderr)
        sys.exit(1)
    drv = drv[:mu.end()] + '\n   use turbogap_md' + drv[mu.end():]
lines = [drv]
open(DRIVER, 'w', encoding='utf-8').write(''.join(lines))
print(f'turbogap_md.f90 written; driver block of {HI-LO+1} lines -> {len(repl)}-line call')
print(f'compute_md takes {len(ARGS)} arguments')
