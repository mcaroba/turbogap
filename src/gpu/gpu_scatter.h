// Deterministic scatter of per-pair contributions onto atoms.
//
// Both force kernels used to close by adding each pair's 3-vector straight into
// forces[j2] with a double atomicAdd, and the virial into nine addresses the
// same way. Floating-point addition is not associative and the order blocks
// reach an address is not fixed, so the same binary on the same input produced
// different last bits from run to run -- invisible at the writer's F16.8, and
// about 1e-12 on the velocities after one MD step, which is where it starts to
// show. See docs/gpu_fixes_handoff.md, "6c revisited".
//
// This replaces the scatter with a sum whose order is a property of the
// neighbour list rather than of the scheduler. Implemented in gpu_scatter.cu.
#ifndef TURBOGAP_GPU_SCATTER_H
#define TURBOGAP_GPU_SCATTER_H

#include "gpu_common.h"

// n_pairs        contributions, one per pair
// n_sites        length of forces_d in atoms; j2_index_d indexes into it
// j2_index_d     destination atom of each pair, 1-based as the Fortran has it
// pair_force_d   3 * n_pairs, the per-pair force, xyz contiguous within a pair
// pair_xyz_d     3 * n_pairs, the pair separation, same layout. Only read when
//                virial_d is given; pass anything when it is not
// forces_d       3 * n_sites, accumulated into
// virial_d       9, accumulated into, or nullptr to skip the virial
// virial_weight  0.5 where the caller only symmetrises, 0.25 where it also
//                halves for the second visit to each unordered pair
void gpu_pair_scatter_reduce(int n_pairs, int n_sites, const int* j2_index_d, const double* pair_force_d,
                             const double* pair_xyz_d, double* forces_d, double* virial_d, double virial_weight,
                             hipStream_t* stream);

#endif // TURBOGAP_GPU_SCATTER_H
