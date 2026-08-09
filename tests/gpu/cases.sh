#!/usr/bin/env bash
#
# The GPU test cases, as data.
#
# This file is sourced, never executed. It exists because two things need the
# same case definitions and they must not drift apart:
#
#   * tests/gpu/run_regression.sh  -- runs every case and checks correctness
#   * tools/profile_gpu.sh         -- runs ONE case under nsys/ncu
#
# A profile is only worth reading if it profiles the thing the suite tests. When
# the input body lived inside run_regression.sh, the only way to profile a case
# was to paste its input somewhere else, and the pasted copy started ageing the
# moment it was made.
#
# Interface
# ---------
#   tg_case_list                 prints every case name, one per line
#   tg_case_def <name>           sets TG_ATOMS TG_MODE TG_DATA TG_BODY TG_EXTRAS
#   tg_case_xfail <name>         prints the xfail reason, empty if the case
#                                is expected to pass
#
# tg_case_def reads two variables the caller must set first:
#   TG_TESTS_DIR   the turbogap_tests clone
#   TG_DATA        the default (CO) data directory inside it
#
# TG_EXTRAS is a bash ARRAY, so a case that needs no extra files leaves it
# empty rather than contributing an empty word to a command line.

tg_case_list() {
  printf '%s\n' CO_predict CO_md vdw_ts XRD_mad estat_gsf
}

# The large systems, listed separately and NOT returned by tg_case_list.
#
# They are not tests: there is no reference to compare a 500,000-atom result
# against, and one of them takes longer than the whole regression suite. The
# suite must never pick them up by iterating the list, which is the entire
# reason for the split -- profile_gpu.sh names them explicitly.
#
# Their data lives outside both repos (like cpu_vs_gpu_tests/input/) because a
# 101 MB xyz does not belong in a git history.
TG_LARGE_DIR=${TURBOGAP_LARGE_DIR:-${TG_TESTS_DIR%/*}/large_systems/diamond_1M}

# Each rung three times: a single point (input.gap), ten steps of MD
# (input.md), and the same ten steps with the experimental-observable machinery
# on (input.mad).
#
# The set is the point, and each pair isolates one cost. A single point pays the
# neighbour build, the descriptor upload and the device allocation once; MD pays
# them once per step, so a profile of one says little about the other. MD and
# MAD then differ only by the pair distribution, structure factor, XRD and
# experimental forces -- and by the cutoff those carry, which is what makes the
# neighbour build and the host arrays grow. Subtracting one from the other is
# the only way to say what the MAD phases actually cost.
tg_case_list_large() {
  local n
  for n in 14k 33k 125k 262k 512k 1M; do
    printf 'diamond_%s\ndiamond_%s_md\ndiamond_%s_mad\n' "$n" "$n" "$n"
  done
}

# The keywords every CO case shares. Kept as one string so a change to the
# potential or the seed cannot reach one case and miss another.
tg_case_common='atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345'

tg_case_def() {
  local name=$1
  TG_ATOMS=; TG_MODE=; TG_BODY=; TG_DATA_DIR=$TG_DATA; TG_EXTRAS=()

  case $name in

  # Single point: exercises the SOAP, 2b and 3b energy, force and virial paths.
  CO_predict)
    TG_ATOMS=atoms_7176.xyz
    TG_MODE=predict
    TG_BODY="$tg_case_common
write_xyz = 1"
    ;;

  # Short MD from explicit velocities, so no randomization is involved and the
  # trajectory is reproducible between the CPU and GPU builds.
  CO_md)
    TG_ATOMS=atoms_7176_vel.xyz
    TG_MODE=md
    TG_BODY="$tg_case_common
md_nsteps = 5
thermostat = berendsen
write_xyz = 1"
    ;;

  # Tkatchenko-Scheffler pairwise vdW on the P4 dimer. Mirrors the CPU branch's
  # vdw_ts case and shares its data, so the two suites cover the same code.
  #
  # Two things make this case worth more than its size. It is the only one here
  # that reaches the vdW path -- which could not run at all on this branch until
  # the deprecated has_vdw GAP form was migrated to a local property (6j) -- and
  # it is the only one that uses a *small* cell. The get_gap_soap overrun in 6i
  # corrupted the heap silently on the 7176-atom CO system and aborted instantly
  # on this one; a suite of a single system size cannot see that class of bug.
  vdw_ts)
    TG_ATOMS=p4_dimer.xyz
    TG_MODE=predict
    TG_DATA_DIR=$TG_TESTS_DIR/vdw_P
    TG_BODY='atoms_file = "atoms.xyz"
pot_file = "gap_files/phosphorus.gap"
n_species = 1
species = P
e0 = -0.52375977
masses = 30.97
random_seed = 12345

vdw_rcut = 25.
vdw_r0_ref = 2.12
vdw_alpha0_ref = 3.7046
vdw_c6_ref = 110.54
vdw_buffer = 0.5

vdw_type = ts
write_xyz = 1'
    ;;

  # Molecular Augmented Dynamics against an experimental XRD pattern. This is
  # the only case that reaches the pair-distribution, structure-factor and XRD
  # paths at all -- the CO cases exercise none of them. It mirrors the CPU
  # branch's xrd_mad case so the two suites cover the same code, and it shares
  # that case's data directory.
  #
  # It is also the only case that turns on gpu_batched for the experimental
  # observables, so it is the one place the per-batch stream machinery in
  # turbogap_exp.f90 runs at all -- which makes it the case to profile when
  # asking whether those streams overlap.
  XRD_mad)
    TG_ATOMS=atoms.xyz
    TG_MODE=md
    TG_DATA_DIR=$(dirname "$TG_DATA")/xrd_mad
    TG_EXTRAS=(xrd_glassy_carbon_zeng_2017.fq)
    TG_BODY='atoms_file = "atoms.xyz"
pot_file = "gap_files/CO.gap"
n_species = 2
species = C O
masses = 12.01 15.99
e0 = -.16138053 0.
random_seed = 12345

do_pair_distribution        = .true.
pair_distribution_kde_sigma =   0.1
pair_distribution_n_samples =  201
pair_distribution_partial   = .true.
pair_distribution_rcut      =   6.0
r_range_min                 =   0.1
r_range_max                 =   6.0

do_structure_factor         = .true.
structure_factor_from_pdf   = .true.
structure_factor_window     = .true.

do_xrd                      = .true.
q_range_min                 =   1.0
q_range_max                 =  10.0
xrd_output                  = "q*F(q)"

n_exp = 1
exp_labels = "xrd"
exp_data_files = "xrd_glassy_carbon_zeng_2017.fq"
exp_n_samples = 201
exp_energy_scales = 100.0
exp_energies = .true.
exp_forces   = .true.

md_nsteps = 5
md_step = 1.0
thermostat = berendsen
t_beg = 1000
t_end = 1000
write_xyz = 1'
    ;;

  # Damped-shifted-force electrostatics with the BATCHED device path, against
  # the CPU build's ordinary one. This is the only case that reaches
  # electrostatics at all -- the CCLi potential is the only one in the test data
  # declaring an atomic_charge local property, and until it was added the path
  # had never run on either branch (docs/refactor_strategy.md section 5).
  #
  # gpu_batched = .true. is the point: calculate_batched_electrostatics is
  # GPU-only, so what this compares is the batched device implementation against
  # the CPU's compute_coulomb_lamichhane. Nothing else in either suite does that.
  #
  # 897-atom carbon cell against the C+CLi potential. Pure carbon is deliberate
  # -- a real equilibrated structure already in the test data, whose C
  # descriptor carries atomic_charge, giving an estat term around -3.2 eV:
  # small, and far enough from zero that a regression in it shows.
  estat_gsf)
    TG_ATOMS=atoms_897.xyz
    TG_MODE=predict
    TG_DATA_DIR=$(dirname "$TG_DATA")/CLi
    TG_BODY='atoms_file = "atoms.xyz"
pot_file = "gap_files/CCLi_estat_ljrep.gap"
n_species = 2
species = C Li
masses = 12.01 6.94
e0 = 0. 0.
random_seed = 12345
estat_method = "gsf"
estat_rcut = 10.0
estat_dsf_alpha = 0.12
estat_damped = .true.
estat_tsf = .true.
estat_sp = .true.
estat_gsf = .true.
estat_self_energy_correction = .false.
gpu_batched = .true.
gpu_n_batches = 4
write_xyz = 1'
    ;;

  # ---------------------------------------------------- the large systems
  #
  # One rung of the diamond ladder each. The atom counts are the ACTUAL counts
  # in the files, not the ideal n^3*8 -- these are cut from a 300 K snapshot, so
  # a few atoms land either side of the cut plane and the file name has to say
  # what is in it rather than what was asked for.
  #
  # The input body is read from input.gap rather than written here. It is the
  # file someone editing this ladder by hand will change, and a copy in this
  # table would be a second thing to keep in step for no benefit.
  diamond_*)
    # diamond_<size>       -> single point, input.gap
    # diamond_<size>_md    -> ten steps of MD, input.md
    # diamond_<size>_mad   -> the same ten steps with pdf / S(q) / XRD /
    #                         experimental forces, input.mad
    #
    # _mad is matched before _md. The two globs happen not to overlap, so the
    # order changes nothing today; it is there to say which one is meant to win
    # if a later name makes them overlap.
    local size=${name#diamond_} file=input.gap
    TG_MODE=predict
    case $size in
    *_mad) size=${size%_mad}; file=input.mad; TG_MODE=md ;;
    *_md) size=${size%_md}; file=input.md; TG_MODE=md ;;
    esac
    case $size in
    14k)  TG_ATOMS=atoms_13824.xyz ;;
    33k)  TG_ATOMS=atoms_32785.xyz ;;
    125k) TG_ATOMS=atoms_124959.xyz ;;
    262k) TG_ATOMS=atoms_262202.xyz ;;
    512k) TG_ATOMS=atoms_511984.xyz ;;
    1M)   TG_ATOMS=atoms_1000000.xyz ;;
    *) printf 'tg_case_def: unknown diamond size %s\n' "$size" >&2; return 1 ;;
    esac
    TG_DATA_DIR=$TG_LARGE_DIR
    if [ -f "$TG_LARGE_DIR/$file" ]; then
      TG_BODY=$(cat "$TG_LARGE_DIR/$file")
    else
      TG_BODY=''
    fi
    ;;

  *)
    printf 'tg_case_def: unknown case %s\n' "$name" >&2
    return 1
    ;;
  esac
}

# tg_stage_case <name> <parent-dir>
#
# Build a runnable directory for one case under <parent-dir> and print its path.
# The data lives in the turbogap_tests clone and is linked, never copied: the CO
# system alone is 7176 atoms and the gap_files directory is shared by three
# cases.
tg_stage_case() {
  local name=$1 parent=$2 dir=$2/$1 extra
  tg_case_def "$name" || return 1
  mkdir -p "$dir"
  ln -sf "$TG_DATA_DIR/$TG_ATOMS" "$dir/atoms.xyz"
  ln -sf "$TG_DATA_DIR/gap_files" "$dir/gap_files"
  for extra in ${TG_EXTRAS+"${TG_EXTRAS[@]}"}; do
    ln -sf "$TG_DATA_DIR/$extra" "$dir/$extra"
  done
  printf '%s\n' "$TG_BODY" >"$dir/input"
  printf '%s' "$dir"
}

# Cases listed here are known to fail for a reason recorded in
# docs/gpu_fixes_handoff.md. They still run and still report, but they do not
# turn the suite red; if one starts passing, the suite says XPASS so the marker
# gets removed rather than quietly masking a regression later.
tg_case_xfail() {
  case $1 in
  XRD_mad) printf '%s' "local_energy drift ~1e-5 by frame 5; everything else agrees; see docs/gpu_fixes_handoff.md 6g" ;;
  estat_gsf) printf '%s' "batched device electrostatics disagrees with the CPU implementation: forces to 1.1 eV/A on |F|max 20.5, virial 0.7%, local_energy 0.16 eV -- every other energy component agrees exactly. Found the first time the path was ever exercised; see the commit that added this case" ;;
  *) printf '' ;;
  esac
}
