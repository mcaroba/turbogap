// Parallel primitives over per-pair flag arrays, shared by the pdf and the
// electrostatics neighbour counting. Implemented in gpu_scan.cu.
#ifndef TURBOGAP_GPU_SCAN_H
#define TURBOGAP_GPU_SCAN_H

#include "gpu_common.h"

// Sum n ints into d_out, reducing recursively until one block remains.
void recursiveReduce(int* d_in, int* d_out, int n, hipStream_t* stream);

// In-place inclusive prefix sum over n ints.
void inclusiveScan(int* d_data_out, int n, hipStream_t* stream);

// Zero every prefix-sum entry whose flag is zero, so the result is a compaction
// index for the flagged pairs and nothing else.
//
// This is a host launcher rather than a kernel the callers launch themselves.
// The kernel lives in gpu_scan.cu, and a <<<>>> launch cannot cross a
// translation unit without -rdc=true; both call sites computed the same
// geometry inline, so it moved here with them.
void gpu_multiply_flags(int n_pairs, int* nk_flags_d, int* nk_sum_flags_d, hipStream_t* stream);

// Print the pending error on a stream without clearing the sticky one.
void gpu_peek_stream_error(hipStream_t* stream);

extern "C" {
void gpu_inclusive_scan_int(int size, int* n_neigh_index_d, hipStream_t* stream);
}

#endif // TURBOGAP_GPU_SCAN_H
