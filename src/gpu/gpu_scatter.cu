// Deterministic scatter of per-pair contributions onto atoms. See
// gpu_scatter.h for why, and docs/gpu_fixes_handoff.md "6c revisited" for how
// the two call sites were found.
//
// The shape of it. A pair's contribution is computed by whichever block owns
// the pair and stored at its own index -- a plain store, so no ordering
// question arises there. Summing those onto atoms then needs, for each atom,
// the list of pairs that name it, in an order that does not depend on the
// scheduler. That list is a counting sort of j2_index:
//
//   count   how many pairs name each atom. Integer atomics, so the counts are
//           exact whatever order they arrive in
//   scan    exclusive prefix sum of the counts, giving each atom's slice
//   place   drop each pair into its atom's slice. The slot comes from an
//           integer atomic, so the slice holds the right pairs but in an
//           arbitrary order
//   order   rank each slice by pair index, which is what removes the
//           arbitrariness
//   gather  sum each slice in that order
//
// The virial is a global sum rather than a per-atom one, so it goes through a
// two-level ordered reduction instead: fixed contiguous chunks of pairs summed
// in index order, then the chunk totals summed in chunk order.
//
// Cost. Five short passes over the pairs, against a force kernel that already
// does an n_soap-long contraction per pair, plus 2 * n_pairs ints and
// 3 * n_pairs doubles of scratch.
#include "gpu_common.h"
#include "gpu_scatter.h"

#define TPB_SCAT 128

// Pairs per block in the scan and in the virial reduction. Both walk their
// chunk sequentially in one thread, so this trades that serial length against
// the number of chunks the second level has to walk.
#define SCAT_CHUNK 256
#define VIRIAL_CHUNK 256

__global__ void kernel_count_incoming(int n_pairs, const int* j2_index_d, int* counts) {
  int l = blockIdx.x * blockDim.x + threadIdx.x;
  if (l < n_pairs)
    atomicAdd(&counts[j2_index_d[l] - 1], 1);
}

// Level one of the scan: each block's total, summed in one thread so the
// arithmetic is plainly what it looks like. Integers, so order is immaterial
// here -- this is written for clarity rather than to protect anything.
__global__ void kernel_chunk_sums(int n_sites, const int* counts, int* chunk_sum) {
  if (threadIdx.x != 0)
    return;
  int beg = blockIdx.x * SCAT_CHUNK;
  int end = min(beg + SCAT_CHUNK, n_sites);
  int s = 0;
  for (int j = beg; j < end; j++)
    s += counts[j];
  chunk_sum[blockIdx.x] = s;
}

// Level two: an exclusive scan of the chunk totals, one thread. n_sites /
// SCAT_CHUNK iterations -- a few thousand for a million atoms.
__global__ void kernel_scan_chunk_sums(int n_chunks, int* chunk_sum) {
  if (threadIdx.x != 0)
    return;
  int running = 0;
  for (int b = 0; b < n_chunks; b++) {
    int c = chunk_sum[b];
    chunk_sum[b] = running;
    running += c;
  }
}

// Level three: the exclusive scan within each chunk, offset by its base.
__global__ void kernel_scan_within(int n_sites, const int* counts, const int* chunk_base, int* begin) {
  if (threadIdx.x != 0)
    return;
  int beg = blockIdx.x * SCAT_CHUNK;
  int end = min(beg + SCAT_CHUNK, n_sites);
  int running = chunk_base[blockIdx.x];
  for (int j = beg; j < end; j++) {
    begin[j] = running;
    running += counts[j];
  }
}

__global__ void kernel_bucket_place(int n_pairs, const int* j2_index_d, const int* begin, int* cursor, int* bucket) {
  int l = blockIdx.x * blockDim.x + threadIdx.x;
  if (l >= n_pairs)
    return;
  int j = j2_index_d[l] - 1;
  bucket[begin[j] + atomicAdd(&cursor[j], 1)] = l;
}

// Rank each atom's slice by pair index. The slice holds the right pairs after
// the placement above but in whatever order the atomics handed out slots;
// ordering it is the whole point of the exercise.
//
// Quadratic in the slice length, which is the atom's neighbour count -- of
// order fifty, and the work is spread over the block's threads.
__global__ void kernel_bucket_order(const int* begin, const int* counts, const int* bucket, int* ordered) {
  int j = blockIdx.x;
  int k = counts[j];
  if (k == 0)
    return;
  int b0 = begin[j];
  for (int t = threadIdx.x; t < k; t += blockDim.x) {
    int v = bucket[b0 + t];
    int rank = 0;
    for (int u = 0; u < k; u++)
      rank += (bucket[b0 + u] < v);
    ordered[b0 + rank] = v;
  }
}

__global__ void kernel_gather_forces(int n_sites, const int* begin, const int* counts, const int* ordered,
                                     const double* pair_force, double* forces) {
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= n_sites)
    return;
  int k = counts[j];
  if (k == 0)
    return;
  int b0 = begin[j];
  double fx = 0.0;
  double fy = 0.0;
  double fz = 0.0;
  for (int r = 0; r < k; r++) {
    int l = ordered[b0 + r];
    fx += pair_force[3 * l];
    fy += pair_force[3 * l + 1];
    fz += pair_force[3 * l + 2];
  }
  forces[3 * j] += fx;
  forces[3 * j + 1] += fy;
  forces[3 * j + 2] += fz;
}

// One thread per virial component, each walking a fixed chunk of pairs in
// index order.
__global__ void kernel_virial_partials(int n_pairs, const double* pair_force, const double* pair_xyz, double w,
                                       double* partials) {
  int c = threadIdx.x;
  if (c >= 9)
    return;
  int k1 = c / 3;
  int k2 = c - 3 * k1;
  int beg = blockIdx.x * VIRIAL_CHUNK;
  int end = min(beg + VIRIAL_CHUNK, n_pairs);
  double s = 0.0;
  for (int l = beg; l < end; l++)
    s += w * (pair_force[3 * l + k1] * pair_xyz[3 * l + k2] + pair_force[3 * l + k2] * pair_xyz[3 * l + k1]);
  partials[9 * blockIdx.x + c] = s;
}

__global__ void kernel_virial_finish(int n_chunks, const double* partials, double* virial) {
  int c = threadIdx.x;
  if (c >= 9)
    return;
  double s = 0.0;
  for (int b = 0; b < n_chunks; b++)
    s += partials[9 * b + c];
  virial[c] += s;
}

void gpu_pair_scatter_reduce(int n_pairs, int n_sites, const int* j2_index_d, const double* pair_force_d,
                             const double* pair_xyz_d, double* forces_d, double* virial_d, double virial_weight,
                             hipStream_t* stream) {
  if (n_pairs <= 0 || n_sites <= 0)
    return;

  int n_chunks = (n_sites + SCAT_CHUNK - 1) / SCAT_CHUNK;

  int* counts;
  int* chunk_sum;
  int* begin;
  int* cursor;
  int* bucket;
  int* ordered;
  gpuErrchk(hipMallocAsync(&counts, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&chunk_sum, (size_t) n_chunks * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&begin, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cursor, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&bucket, (size_t) n_pairs * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&ordered, (size_t) n_pairs * sizeof(int), stream[0]));

  gpuErrchk(hipMemsetAsync(counts, 0, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMemsetAsync(cursor, 0, (size_t) n_sites * sizeof(int), stream[0]));

  dim3 tb(TPB_SCAT);
  dim3 g_pairs((n_pairs + TPB_SCAT - 1) / TPB_SCAT);
  dim3 g_sites((n_sites + TPB_SCAT - 1) / TPB_SCAT);
  dim3 g_chunks(n_chunks);

  kernel_count_incoming<<<g_pairs, tb, 0, stream[0]>>>(n_pairs, j2_index_d, counts);
  kernel_chunk_sums<<<g_chunks, dim3(1), 0, stream[0]>>>(n_sites, counts, chunk_sum);
  kernel_scan_chunk_sums<<<dim3(1), dim3(1), 0, stream[0]>>>(n_chunks, chunk_sum);
  kernel_scan_within<<<g_chunks, dim3(1), 0, stream[0]>>>(n_sites, counts, chunk_sum, begin);

  kernel_bucket_place<<<g_pairs, tb, 0, stream[0]>>>(n_pairs, j2_index_d, begin, cursor, bucket);
  kernel_bucket_order<<<dim3(n_sites), tb, 0, stream[0]>>>(begin, counts, bucket, ordered);
  kernel_gather_forces<<<g_sites, tb, 0, stream[0]>>>(n_sites, begin, counts, ordered, pair_force_d, forces_d);

  if (virial_d != nullptr) {
    int n_vchunks = (n_pairs + VIRIAL_CHUNK - 1) / VIRIAL_CHUNK;
    double* partials;
    gpuErrchk(hipMallocAsync(&partials, (size_t) n_vchunks * 9 * sizeof(double), stream[0]));
    kernel_virial_partials<<<dim3(n_vchunks), dim3(WARP_SIZE), 0, stream[0]>>>(n_pairs, pair_force_d, pair_xyz_d,
                                                                              virial_weight, partials);
    kernel_virial_finish<<<dim3(1), dim3(WARP_SIZE), 0, stream[0]>>>(n_vchunks, partials, virial_d);
    gpuErrchk(hipFreeAsync(partials, stream[0]));
  }

  gpuErrchk(hipFreeAsync(ordered, stream[0]));
  gpuErrchk(hipFreeAsync(bucket, stream[0]));
  gpuErrchk(hipFreeAsync(cursor, stream[0]));
  gpuErrchk(hipFreeAsync(begin, stream[0]));
  gpuErrchk(hipFreeAsync(chunk_sum, stream[0]));
  gpuErrchk(hipFreeAsync(counts, stream[0]));
}
