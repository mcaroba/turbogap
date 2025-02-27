mkdir weak_scaling
mkdir strong_scaling

for i in 1 2 4 8 16 32 64; do

echo $i

if [ $i -le 4 ]; then
natoms=$(($i*96))
string=weak_scaling/n_1_g_${i}_at_${natoms}
mkdir $string
cd $string
cp -r ../../gap_files .
cp ../../BP${natoms}.xyz .
cat <<EOF > input
atoms_file = 'BP${natoms}.xyz'
pot_file = 'gap_files/phosphorus.gap'

! Species info
n_species = 1 !2
species = P !C !H
!e0 = -.16138053
e0 = -0.52375977
masses = 30.97 !12.01 !1.00784

! van der Waals info
!vdw_type = ts
vdw_type = mbd
vdw_sr = 0.83 !0.97
vdw_d = 6. !20.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 17.
vdw_r0_ref = 2.12 !1.900 !1.64
vdw_alpha0_ref = 3.7046 !1.778 !0.667
vdw_c6_ref = 110.54 !27.8 !3.88
! The order of the cutoffs is vdw_scs_rcut < vdw_mbd_rcut < vdw_2b_rcut
! Also: vdw_mbd_rcut => vdw_mbd_rcut2, vdw_2b_rcut => vdw_2b_rcut2
! If two cutoffs are equal, there is no buffer region between them!
! If you are doing 2b with SCS polarizabilities (basically TS-SCS) and need only energies,
! you can set vdw_mbd_rcut = vdw_mbd_rcut2 = vdw_loc_rcut = vdw_2b_rcut2 = 0.
! Then you get away with neighbor lists of same size as normal TS.
! For forces you need vdw_2b_rcut2 > 0 and vdw_loc_rcut > 0 which require larger
! neighbor lists: vdw_loc_rcut = vdw_scs_rcut and vdw_2b_rcut2 = vdw_2b_rcut
! to get results that are equivalent to normal TS.
vdw_buffer = 0.5          ! Buffer for transitions between cut-off regions. Type: REAL
vdw_scs_rcut = 5.         ! SCS cut-off for polarizabilities. Type: REAL
vdw_loc_rcut = 5.         ! Cut-off for including polarizability gradients in force calculation (vdw_loc_rcut <= vdw_scs_rcut).
vdw_mbd_rcut = 10.0         ! Cut-off for atoms to include for local MBD energy (cut-off for central atom). Type: REAL
vdw_mbd_rcut2 = 6.0        ! Secondary MBD cut-off (cut-off sphere for neighbors)
vdw_2b_rcut = 10.0          ! Primary 2b cut-off (cut-off for central atom). Type: REAL
vdw_2b_rcut2 = 6.0         ! Secondary 2b cut-off (cut-off sphere for neighbors).
vdw_mbd_nfreq = 13        ! Number of frequency values for MBD integration. Type: INT
vdw_mbd_norder = 6        ! Contributions up to n-body interactions (i.e. cut-off degree for Taylor expansion of ln(I-AT)). Type: INT
vdw_mbd_grad = .true.     ! Calculate MBD forces. Type: LOGICAL
vdw_hirsh_grad = .true.   ! Include Hirshfeld gradients in the forces. Type: LOGICAL
vdw_polynomial = .false.  ! Use polynomial approximation for inverse matrices. Type: LOGICAL. OBSOLETE!!!
vdw_omega_ref = 1.3d0     ! Reference frequency value used for SCS characteristic frequency calculation. Might depend on the material (for now).
do_nnls = .false.         ! NNLS for frequency fitting to avoid calculating many frequency points. Use at least vdw_mbd_nfreq = 2*vdw_mbd_norder+1 or it might be unstable (overfitting)
vdw_mbd_cent_appr = .true.
vdw_mbd_gpu = .true.
poly_cut_xmin = 3.d0
poly_cut_xmax = 10.d0
EOF
cat <<EOF > qscript.sh
#!/bin/bash -l
#SBATCH --job-name=n_1_g_${i}_at_${natoms}
#SBATCH --account=project_2002883
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:${i}
#SBATCH --ntasks-per-node=${i}
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -t 0-2:00:00
#SBATCH -o out_files/out_n_1_g_${i}_at_${natoms}
#SBATCH --array=1


module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop


srun /projappl/project_2002883/turbogap_mbd/bin/turbogap predict
EOF
sbatch qscript.sh
cd ../..
string=strong_scaling/n_1_g_${i}
mkdir $string
cd $string
cp -r ../../gap_files .
cp ../../BP96.xyz .
cat <<EOF > input
atoms_file = 'BP96.xyz'
pot_file = 'gap_files/phosphorus.gap'

! Species info
n_species = 1 !2
species = P !C !H
!e0 = -.16138053
e0 = -0.52375977
masses = 30.97 !12.01 !1.00784

! van der Waals info
!vdw_type = ts
vdw_type = mbd
vdw_sr = 0.83 !0.97
vdw_d = 6. !20.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 17.
vdw_r0_ref = 2.12 !1.900 !1.64
vdw_alpha0_ref = 3.7046 !1.778 !0.667
vdw_c6_ref = 110.54 !27.8 !3.88
! The order of the cutoffs is vdw_scs_rcut < vdw_mbd_rcut < vdw_2b_rcut
! Also: vdw_mbd_rcut => vdw_mbd_rcut2, vdw_2b_rcut => vdw_2b_rcut2
! If two cutoffs are equal, there is no buffer region between them!
! If you are doing 2b with SCS polarizabilities (basically TS-SCS) and need only energies,
! you can set vdw_mbd_rcut = vdw_mbd_rcut2 = vdw_loc_rcut = vdw_2b_rcut2 = 0.
! Then you get away with neighbor lists of same size as normal TS.
! For forces you need vdw_2b_rcut2 > 0 and vdw_loc_rcut > 0 which require larger
! neighbor lists: vdw_loc_rcut = vdw_scs_rcut and vdw_2b_rcut2 = vdw_2b_rcut
! to get results that are equivalent to normal TS.
vdw_buffer = 0.5          ! Buffer for transitions between cut-off regions. Type: REAL
vdw_scs_rcut = 5.         ! SCS cut-off for polarizabilities. Type: REAL
vdw_loc_rcut = 5.         ! Cut-off for including polarizability gradients in force calculation (vdw_loc_rcut <= vdw_scs_rcut).
vdw_mbd_rcut = 10.0         ! Cut-off for atoms to include for local MBD energy (cut-off for central atom). Type: REAL
vdw_mbd_rcut2 = 6.0        ! Secondary MBD cut-off (cut-off sphere for neighbors)
vdw_2b_rcut = 10.0          ! Primary 2b cut-off (cut-off for central atom). Type: REAL
vdw_2b_rcut2 = 6.0         ! Secondary 2b cut-off (cut-off sphere for neighbors).
vdw_mbd_nfreq = 13        ! Number of frequency values for MBD integration. Type: INT
vdw_mbd_norder = 6        ! Contributions up to n-body interactions (i.e. cut-off degree for Taylor expansion of ln(I-AT)). Type: INT
vdw_mbd_grad = .true.     ! Calculate MBD forces. Type: LOGICAL
vdw_hirsh_grad = .true.   ! Include Hirshfeld gradients in the forces. Type: LOGICAL
vdw_polynomial = .false.  ! Use polynomial approximation for inverse matrices. Type: LOGICAL. OBSOLETE!!!
vdw_omega_ref = 1.3d0     ! Reference frequency value used for SCS characteristic frequency calculation. Might depend on the material (for now).
do_nnls = .false.         ! NNLS for frequency fitting to avoid calculating many frequency points. Use at least vdw_mbd_nfreq = 2*vdw_mbd_norder+1 or it might be unstable (overfitting)
vdw_mbd_cent_appr = .true.
vdw_mbd_gpu = .true.
poly_cut_xmin = 3.d0
poly_cut_xmax = 10.d0
EOF
cat <<EOF > qscript.sh
#!/bin/bash -l
#SBATCH --job-name=n_1_g_${i}
#SBATCH --account=project_2002883
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:${i}
#SBATCH --ntasks-per-node=${i}
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -t 0-2:00:00
#SBATCH -o out_files/out_n_1_g_${i}
#SBATCH --array=1


module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop


srun /projappl/project_2002883/turbogap_mbd/bin/turbogap predict
EOF
sbatch qscript.sh
cd ../..
fi

if [ $i -ge 8 ]; then
natoms=$(($i*96))
n=$(($i/4))
string=weak_scaling/n_${n}_g_4_at_${natoms}
mkdir $string
cd $string
cp -r ../../gap_files .
cp ../../BP${natoms}.xyz .
cat <<EOF > input
atoms_file = 'BP${natoms}.xyz'
pot_file = 'gap_files/phosphorus.gap'

! Species info
n_species = 1 !2
species = P !C !H
!e0 = -.16138053
e0 = -0.52375977
masses = 30.97 !12.01 !1.00784

! van der Waals info
!vdw_type = ts
vdw_type = mbd
vdw_sr = 0.83 !0.97
vdw_d = 6. !20.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 17.
vdw_r0_ref = 2.12 !1.900 !1.64
vdw_alpha0_ref = 3.7046 !1.778 !0.667
vdw_c6_ref = 110.54 !27.8 !3.88
! The order of the cutoffs is vdw_scs_rcut < vdw_mbd_rcut < vdw_2b_rcut
! Also: vdw_mbd_rcut => vdw_mbd_rcut2, vdw_2b_rcut => vdw_2b_rcut2
! If two cutoffs are equal, there is no buffer region between them!
! If you are doing 2b with SCS polarizabilities (basically TS-SCS) and need only energies,
! you can set vdw_mbd_rcut = vdw_mbd_rcut2 = vdw_loc_rcut = vdw_2b_rcut2 = 0.
! Then you get away with neighbor lists of same size as normal TS.
! For forces you need vdw_2b_rcut2 > 0 and vdw_loc_rcut > 0 which require larger
! neighbor lists: vdw_loc_rcut = vdw_scs_rcut and vdw_2b_rcut2 = vdw_2b_rcut
! to get results that are equivalent to normal TS.
vdw_buffer = 0.5          ! Buffer for transitions between cut-off regions. Type: REAL
vdw_scs_rcut = 5.         ! SCS cut-off for polarizabilities. Type: REAL
vdw_loc_rcut = 5.         ! Cut-off for including polarizability gradients in force calculation (vdw_loc_rcut <= vdw_scs_rcut).
vdw_mbd_rcut = 10.0         ! Cut-off for atoms to include for local MBD energy (cut-off for central atom). Type: REAL
vdw_mbd_rcut2 = 6.0        ! Secondary MBD cut-off (cut-off sphere for neighbors)
vdw_2b_rcut = 10.0          ! Primary 2b cut-off (cut-off for central atom). Type: REAL
vdw_2b_rcut2 = 6.0         ! Secondary 2b cut-off (cut-off sphere for neighbors).
vdw_mbd_nfreq = 13        ! Number of frequency values for MBD integration. Type: INT
vdw_mbd_norder = 6        ! Contributions up to n-body interactions (i.e. cut-off degree for Taylor expansion of ln(I-AT)). Type: INT
vdw_mbd_grad = .true.     ! Calculate MBD forces. Type: LOGICAL
vdw_hirsh_grad = .true.   ! Include Hirshfeld gradients in the forces. Type: LOGICAL
vdw_polynomial = .false.  ! Use polynomial approximation for inverse matrices. Type: LOGICAL. OBSOLETE!!!
vdw_omega_ref = 1.3d0     ! Reference frequency value used for SCS characteristic frequency calculation. Might depend on the material (for now).
do_nnls = .false.         ! NNLS for frequency fitting to avoid calculating many frequency points. Use at least vdw_mbd_nfreq = 2*vdw_mbd_norder+1 or it might be unstable (overfitting)
vdw_mbd_cent_appr = .true.
vdw_mbd_gpu = .true.
poly_cut_xmin = 3.d0
poly_cut_xmax = 10.d0
EOF
cat <<EOF > qscript.sh
#!/bin/bash -l
#SBATCH --job-name=n_${n}_g_4_at_${natoms}
#SBATCH --account=project_2002883
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:4
#SBATCH --ntasks-per-node=4
#SBATCH --exclusive
#SBATCH --nodes=${n}
#SBATCH -t 0-2:00:00
#SBATCH -o out_files/out_n_${n}_g_4_at_${natoms}
#SBATCH --array=1


module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop


srun /projappl/project_2002883/turbogap_mbd/bin/turbogap predict
EOF
sbatch qscript.sh
cd ../..
string=strong_scaling/n_${n}_g_4
mkdir $string
cd $string
cp -r ../../gap_files .
cp ../../BP96.xyz .
cat <<EOF > input
atoms_file = 'BP96.xyz'
pot_file = 'gap_files/phosphorus.gap'

! Species info
n_species = 1 !2
species = P !C !H
!e0 = -.16138053
e0 = -0.52375977
masses = 30.97 !12.01 !1.00784

! van der Waals info
!vdw_type = ts
vdw_type = mbd
vdw_sr = 0.83 !0.97
vdw_d = 6. !20.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 17.
vdw_r0_ref = 2.12 !1.900 !1.64
vdw_alpha0_ref = 3.7046 !1.778 !0.667
vdw_c6_ref = 110.54 !27.8 !3.88
! The order of the cutoffs is vdw_scs_rcut < vdw_mbd_rcut < vdw_2b_rcut
! Also: vdw_mbd_rcut => vdw_mbd_rcut2, vdw_2b_rcut => vdw_2b_rcut2
! If two cutoffs are equal, there is no buffer region between them!
! If you are doing 2b with SCS polarizabilities (basically TS-SCS) and need only energies,
! you can set vdw_mbd_rcut = vdw_mbd_rcut2 = vdw_loc_rcut = vdw_2b_rcut2 = 0.
! Then you get away with neighbor lists of same size as normal TS.
! For forces you need vdw_2b_rcut2 > 0 and vdw_loc_rcut > 0 which require larger
! neighbor lists: vdw_loc_rcut = vdw_scs_rcut and vdw_2b_rcut2 = vdw_2b_rcut
! to get results that are equivalent to normal TS.
vdw_buffer = 0.5          ! Buffer for transitions between cut-off regions. Type: REAL
vdw_scs_rcut = 5.         ! SCS cut-off for polarizabilities. Type: REAL
vdw_loc_rcut = 5.         ! Cut-off for including polarizability gradients in force calculation (vdw_loc_rcut <= vdw_scs_rcut).
vdw_mbd_rcut = 10.0         ! Cut-off for atoms to include for local MBD energy (cut-off for central atom). Type: REAL
vdw_mbd_rcut2 = 6.0        ! Secondary MBD cut-off (cut-off sphere for neighbors)
vdw_2b_rcut = 10.0          ! Primary 2b cut-off (cut-off for central atom). Type: REAL
vdw_2b_rcut2 = 6.0         ! Secondary 2b cut-off (cut-off sphere for neighbors).
vdw_mbd_nfreq = 13        ! Number of frequency values for MBD integration. Type: INT
vdw_mbd_norder = 6        ! Contributions up to n-body interactions (i.e. cut-off degree for Taylor expansion of ln(I-AT)). Type: INT
vdw_mbd_grad = .true.     ! Calculate MBD forces. Type: LOGICAL
vdw_hirsh_grad = .true.   ! Include Hirshfeld gradients in the forces. Type: LOGICAL
vdw_polynomial = .false.  ! Use polynomial approximation for inverse matrices. Type: LOGICAL. OBSOLETE!!!
vdw_omega_ref = 1.3d0     ! Reference frequency value used for SCS characteristic frequency calculation. Might depend on the material (for now).
do_nnls = .false.         ! NNLS for frequency fitting to avoid calculating many frequency points. Use at least vdw_mbd_nfreq = 2*vdw_mbd_norder+1 or it might be unstable (overfitting)
vdw_mbd_cent_appr = .true.
vdw_mbd_gpu = .true.
poly_cut_xmin = 3.d0
poly_cut_xmax = 10.d0
EOF
cat <<EOF > qscript.sh
#!/bin/bash -l
#SBATCH --job-name=n_${n}_g_4
#SBATCH --account=project_2002883
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:4
#SBATCH --ntasks-per-node=4
#SBATCH --exclusive
#SBATCH --nodes=${n}
#SBATCH -t 0-2:00:00
#SBATCH -o out_files/out_n_${n}_g_4
#SBATCH --array=1


module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop


srun /projappl/project_2002883/turbogap_mbd/bin/turbogap predict
EOF
sbatch qscript.sh
cd ../..
fi

done

