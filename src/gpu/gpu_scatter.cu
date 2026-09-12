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
// does an n_soap-long contraction per pair, plus 3 * n_pairs ints and
// 3 * n_pairs doubles of scratch.
//
// TWO BACKENDS. Every pass here is a flat map -- one index, one output slot --
// written through tg_parallel_for, so it compiles to a <<<>>> launch or to
// Kokkos::parallel_for from one body. See src/gpu/gpu_backend.h.
//
// That conversion is safe for the determinism argument above precisely because
// none of the ordered sums are parallel reductions: each runs SEQUENTIALLY
// INSIDE ONE THREAD, over a range fixed by the neighbour list. Changing how a
// thread is handed its index cannot reorder them. Do not be tempted to replace
// the gather or the virial with a parallel_reduce -- that would reassociate,
// and hand back exactly the run-to-run wobble this file exists to remove.
//
// The one pass that changed shape is the ordering. It used to be one block per
// atom with the block's threads striding over that atom's slice, which is a
// team pattern; it is now a flat map over slots, which needs to know which
// atom owns each slot. That is the `owner` array, filled for free during the
// placement pass, and it is why the int scratch went from two to three per
// pair. Same comparisons, same writes, same result.
#include "gpu_backend.h"
#include "gpu_common.h"
#include "gpu_scatter.h"

#define TPB_SCAT 128

// Pairs per block in the scan and in the virial reduction. Both walk their
// chunk sequentially in one thread, so this trades that serial length against
// the number of chunks the second level has to walk.
#define SCAT_CHUNK 256
#define VIRIAL_CHUNK 256

void gpu_pair_scatter_reduce(int n_pairs, int n_sites, const int* j2_index_d, const double* pair_force_d, const double* pair_xyz_d,
                             double* forces_d, double* virial_d, double virial_weight, hipStream_t* stream) {
  if (n_pairs <= 0 || n_sites <= 0)
    return;

  const int n_chunks = (n_sites + SCAT_CHUNK - 1) / SCAT_CHUNK;

  int* counts;
  int* chunk_sum;
  int* begin;
  int* cursor;
  int* bucket;
  int* ordered;
  int* owner;
  gpuErrchk(hipMallocAsync(&counts, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&chunk_sum, (size_t) n_chunks * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&begin, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&cursor, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&bucket, (size_t) n_pairs * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&ordered, (size_t) n_pairs * sizeof(int), stream[0]));
  gpuErrchk(hipMallocAsync(&owner, (size_t) n_pairs * sizeof(int), stream[0]));

  gpuErrchk(hipMemsetAsync(counts, 0, (size_t) n_sites * sizeof(int), stream[0]));
  gpuErrchk(hipMemsetAsync(cursor, 0, (size_t) n_sites * sizeof(int), stream[0]));

  // count: how many pairs name each atom. Integer atomics, so exact whatever
  // order they arrive in.
  tg_parallel_for(
      "turbogap_scatter_count", n_pairs, stream, TG_LAMBDA(const int l) { TG_ATOMIC_ADD(&counts[j2_index_d[l] - 1], 1); },
      TPB_SCAT);

  // scan, level one: each chunk's total, summed in one thread so the
  // arithmetic is plainly what it looks like.
  tg_parallel_for(
      "turbogap_scatter_chunk_sums", n_chunks, stream,
      TG_LAMBDA(const int b) {
        const int beg = b * SCAT_CHUNK;
        const int end = (beg + SCAT_CHUNK < n_sites) ? beg + SCAT_CHUNK : n_sites;
        int s = 0;
        for (int j = beg; j < end; j++)
          s += counts[j];
        chunk_sum[b] = s;
      },
      TPB_SCAT);

  // level two: an exclusive scan of the chunk totals, one thread. n_sites /
  // SCAT_CHUNK iterations -- a few thousand for a million atoms.
  tg_parallel_for(
      "turbogap_scatter_scan_chunks", 1, stream,
      TG_LAMBDA(const int) {
        int running = 0;
        for (int b = 0; b < n_chunks; b++) {
          const int c = chunk_sum[b];
          chunk_sum[b] = running;
          running += c;
        }
      },
      1);

  // level three: the exclusive scan within each chunk, offset by its base.
  tg_parallel_for(
      "turbogap_scatter_scan_within", n_chunks, stream,
      TG_LAMBDA(const int b) {
        const int beg = b * SCAT_CHUNK;
        const int end = (beg + SCAT_CHUNK < n_sites) ? beg + SCAT_CHUNK : n_sites;
        int running = chunk_sum[b];
        for (int j = beg; j < end; j++) {
          begin[j] = running;
          running += counts[j];
        }
      },
      TPB_SCAT);

  // place: drop each pair into its atom's slice. The slot comes from an
  // integer atomic, so the slice holds the right pairs in an arbitrary order.
  // owner records which atom each slot belongs to, which is what lets the
  // ordering pass below be a flat map rather than one block per atom.
  tg_parallel_for(
      "turbogap_scatter_place", n_pairs, stream,
      TG_LAMBDA(const int l) {
        const int j = j2_index_d[l] - 1;
        const int slot = begin[j] + TG_ATOMIC_FETCH_ADD(&cursor[j], 1);
        bucket[slot] = l;
        owner[slot] = j;
      },
      TPB_SCAT);

  // order: rank each atom's slice by pair index, which is what removes the
  // arbitrariness the atomics introduced. Quadratic in the slice length --
  // the atom's neighbour count, of order fifty.
  tg_parallel_for(
      "turbogap_scatter_order", n_pairs, stream,
      TG_LAMBDA(const int s) {
        const int j = owner[s];
        const int b0 = begin[j];
        const int k = counts[j];
        const int v = bucket[s];
        int rank = 0;
        for (int u = 0; u < k; u++)
          rank += (bucket[b0 + u] < v);
        ordered[b0 + rank] = v;
      },
      TPB_SCAT);

  // gather: sum each slice in that order.
  tg_parallel_for(
      "turbogap_scatter_gather", n_sites, stream,
      TG_LAMBDA(const int j) {
        const int k = counts[j];
        if (k == 0)
          return;
        const int b0 = begin[j];
        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;
        for (int r = 0; r < k; r++) {
          const int l = ordered[b0 + r];
          fx += pair_force_d[3 * l];
          fy += pair_force_d[3 * l + 1];
          fz += pair_force_d[3 * l + 2];
        }
        forces_d[3 * j] += fx;
        forces_d[3 * j + 1] += fy;
        forces_d[3 * j + 2] += fz;
      },
      TPB_SCAT);

  if (virial_d != nullptr) {
    const int n_vchunks = (n_pairs + VIRIAL_CHUNK - 1) / VIRIAL_CHUNK;
    double* partials;
    gpuErrchk(hipMallocAsync(&partials, (size_t) n_vchunks * 9 * sizeof(double), stream[0]));

    // One thread per (chunk, component), each walking its chunk in index
    // order. Flattened from the old grid.x = chunk, thread.x = component.
    tg_parallel_for(
        "turbogap_scatter_virial_partials", n_vchunks * 9, stream,
        TG_LAMBDA(const int idx) {
          const int b = idx / 9;
          const int c = idx - 9 * b;
          const int k1 = c / 3;
          const int k2 = c - 3 * k1;
          const int beg = b * VIRIAL_CHUNK;
          const int end = (beg + VIRIAL_CHUNK < n_pairs) ? beg + VIRIAL_CHUNK : n_pairs;
          double s = 0.0;
          for (int l = beg; l < end; l++)
            s += virial_weight *
                 (pair_force_d[3 * l + k1] * pair_xyz_d[3 * l + k2] + pair_force_d[3 * l + k2] * pair_xyz_d[3 * l + k1]);
          partials[9 * b + c] = s;
        },
        TPB_SCAT);

    tg_parallel_for(
        "turbogap_scatter_virial_finish", 9, stream,
        TG_LAMBDA(const int c) {
          double s = 0.0;
          for (int b = 0; b < n_vchunks; b++)
            s += partials[9 * b + c];
          virial_d[c] += s;
        },
        9);

    gpuErrchk(hipFreeAsync(partials, stream[0]));
  }

  gpuErrchk(hipFreeAsync(owner, stream[0]));
  gpuErrchk(hipFreeAsync(ordered, stream[0]));
  gpuErrchk(hipFreeAsync(bucket, stream[0]));
  gpuErrchk(hipFreeAsync(cursor, stream[0]));
  gpuErrchk(hipFreeAsync(begin, stream[0]));
  gpuErrchk(hipFreeAsync(chunk_sum, stream[0]));
  gpuErrchk(hipFreeAsync(counts, stream[0]));
}
