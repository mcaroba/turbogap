// Device and pinned-host memory, and device selection.
//
// Every one of these is called from Fortran through bind(C) in
// src/fortran_cuda_interfaces.f90; the declarations exist so the C++ side
// cannot drift from the definitions in gpu_memory.cu without a compile error.

#ifndef TURBOGAP_GPU_MEMORY_H
#define TURBOGAP_GPU_MEMORY_H

#include "gpu_common.h"

extern "C" {

// ---- gpu_memory.cu
void gpu_mem_report(const char* label);
void gpu_mem_stats(size_t* dev_free, size_t* dev_total, size_t* live, size_t* peak);
void gpu_mem_set_hint(double max_gb, int n_batches);
void gpu_malloc_async(void** a_d, size_t Np, hipStream_t* stream);
void cuda_memset_async(void* a_d, int value, size_t Np, hipStream_t* stream);
void cuda_malloc_all_blocking(void** a_d, size_t Np);
void cuda_device_reset();
void cuda_free(void** a_d);
void cuda_free_async(void** a_d, hipStream_t* stream);
void cuda_cpy_htod(void* a, void* a_d, size_t N, hipStream_t* stream);
void cuda_cpy_htod_blocking(void* a, void* a_d, size_t N);
void cuda_cpy_dtod(void* b_d, void* a_d, size_t N, hipStream_t* stream);
void cuda_cpy_dtoh(void* a_d, void* a, size_t N, hipStream_t* stream);
void cuda_cpy_dtoh_event(void* a_d, void* a, size_t N, hipStream_t* stream);
void cuda_cpy_dtoh_blocking(void* a_d, void* a, size_t N);
int gpu_device_count();
void cuda_set_device(int my_rank);
void gpu_device_sync();
void gpu_stream_sync(hipStream_t* stream);
int* host_alloc(size_t size);
void host_free(void* ptr);
void cpy_htoh_pinned(void* src, void* dest, size_t size);
void gpu_check_error();
void gpu_meminfo();
void gpu_print_pointer_int(int* p);
void gpu_print_pointer_double(double* p);

} // extern "C"

#endif // TURBOGAP_GPU_MEMORY_H
