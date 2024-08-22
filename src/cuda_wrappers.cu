#include "hip/hip_runtime.h"
// wrappers file
// compile with:
// rm cuda_wrappers.o; nvcc -lcublas -lcurand -arch=sm_80 src/cuda_wrappers.cu -c;
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <cmath>

#include <hip/hip_runtime.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
//#include <hipsolver.h>
#include <hiprand/hiprand.h>
#include <assert.h>
#include <hip/hip_complex.h>

#define tpb 64 // optimize for best performance & check the effect on each kernel which used tpb for the shared memory
#define tpb_get_soap_der_one 128
#define tpbcnk 64 // this is because k_max is 45???
#define static_alphamax 8

#define mode_polynomial 1

int counter=0;
/*__device__ double atomicDoubleAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}*/


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(hipError_t code, const char *file, int line, bool abort=true)
{
   if (code != hipSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ double N_a (double rcut, int a) {
   const int b = 2*a + 5;
   return sqrt( rcut / static_cast<double>(b) );
}


__device__ double gpu_spline(int nx, double r, double *x, double *y, double *y2, double rcut, 
		           double yp1, double ypn){
   int j;
   double s = 0.0;
   if( r < rcut ){
     if( r < x[0] ){
       s = y[0] + (r - x[0]) * yp1;
     } else if( r > x[nx-1] ){
       s = y[nx-1] + (r - x[nx-1]) * ypn;
     } else {
       for(j=0;j<nx-1; j++)
         if( r < x[j+1] ) break;
       double h = x[j+1] - x[j];
       double h26 = pow(h,2) / 6.0;
       double A = (x[j+1] - r) / h;
       double B = 1.0 - A;
       double C = (pow(A,3) - A) * h26;
       double D = (pow(B,3) - B) * h26;
       s = A*y[j] + B*y[j+1] + C*y2[j] + D*y2[j+1];
     }
   }
   return s;
}

__device__ double gpu_spline_der(int nx, double r, double *x, double *y, double *y2, double rcut,
                           double yp1, double ypn){
   int j;
   double ds = 0.0;
   if( r < rcut ){
     if( r < x[0] ){
       ds = yp1;
     } else if( r > x[nx-1] ){
       ds = ypn;
     } else {
       for(j=0;j<nx-1; j++)
         if( r < x[j+1] ) break;
       double h = x[j+1] - x[j];
       double h6 = h / 6.0;
       double A = (x[j+1] - r) / h;
       double B = 1.0 - A;
       double dAdx = -1.0 / h;
       double dBdx = -dAdx;
       double dCdx = (1.0 - 3.0*pow(A,2)) * h6;
       double dDdx = (3.0*pow(B,2) - 1.0) * h6;
       ds = dAdx*y[j] + dBdx*y[j+1] + dCdx*y2[j] + dDdx*y2[j+1];
     }
   }
   return ds;
}
   
__global__ void vect_dble(double *a, int N)
{
   int idx = threadIdx.x+blockIdx.x*gridDim.x;
   if (idx<N)printf(" %lf \n", a[idx]);
}


extern "C" void cuda_malloc_all(void **a_d, size_t Np, hipStream_t *stream )
{
  

  gpuErrchk(hipMallocAsync((void **) a_d,  Np ,stream[0]));
  //gpuErrchk(hipMalloc((void **) a_d,  Np ));
  hipError_t err;
  hipDeviceSynchronize();
  err = hipGetLastError();
//  if (err != hipSuccess) {  
//} 
   return;
}

//  extern "C" void cuda_malloc_all(void **a_d, size_t Np, hipStream_t *stream )
// {
  

//   //gpuErrchk(hipMallocAsync((void **) a_d,  Np ,stream[0]));
//   gpuErrchk(hipMalloc((void **) a_d,  Np ));
//   hipError_t err;
//   hipDeviceSynchronize();
//   err = hipGetLastError();
// //  if (err != hipSuccess) {  
// //} 
//    return;
// }

extern "C" void cuda_memset_async(void *a_d, int value,  size_t Np, hipStream_t *stream )
{
  hipMemsetAsync( a_d, value , Np ,stream[0]);
}
extern "C" void cuda_malloc_all_blocking(void **a_d, size_t Np)
{
  gpuErrchk(hipMalloc( (void **) a_d,  Np ));
   return;
}
extern "C" void cuda_device_reset(){
  hipDeviceReset();
}
extern "C" void cuda_free(void **a_d)
{
  gpuErrchk(hipFree(*a_d));
   //printf("GPU memory freed \n");
   //hipDeviceReset();
   return;
}


// extern "C" void cuda_free_async(void **a_d, hipStream_t *stream )
// {
//   //gpuErrchk(hipFreeAsync(*a_d, stream[0]));
//   gpuErrchk(hipFree(*a_d));
//    return;
// }

extern "C" void cuda_free_async(void **a_d, hipStream_t *stream )
{
  gpuErrchk(hipFreeAsync(*a_d, stream[0]));
   return;
}

/*extern "C" void GPU_fill_rand(double *A, int N, int ccc) {
	// Create a pseudo-random number generator
	hiprandGenerator_t prng;
	hiprandCreateGenerator(&prng, HIPRAND_RNG_PSEUDO_DEFAULT);

	// Set the seed for the random number generator using the system clock
	hiprandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock() + (unsigned long long)  ccc * N);

	// Fill the array with random numbers on the device
	hiprandGenerateUniformDouble(prng, A,N);
  //vect_dble<<<(N+128-1)/128,128>>>(A,N);
  //hipDeviceSynchronize();
  printf("\n Filled \n");
}
*/

extern "C" void cuda_cpy_htod(void *a, void *a_d, size_t N, hipStream_t *stream )
{
  gpuErrchk(hipMemcpyAsync(a_d, a, N, hipMemcpyHostToDevice,stream[0] ));
  //gpuErrchk(hipMemcpy(a_d, a, N, hipMemcpyHostToDevice));
   return;
}
extern "C" void cuda_cpy_htod_blocking(void *a, void *a_d, size_t N)
{;
  gpuErrchk(hipMemcpy(a_d, a, N, hipMemcpyHostToDevice));
   return;
}

extern "C" void cuda_cpy_dtod(void *b_d, void *a_d,size_t N, hipStream_t* stream )
{
  gpuErrchk(hipMemcpyAsync( a_d, b_d, N, hipMemcpyDeviceToDevice,stream[0]));
  //gpuErrchk(hipMemcpy( a_d, b_d, N, hipMemcpyDeviceToDevice));
   return;
}

extern "C" void cuda_cpy_dtoh(void *a_d, void *a, size_t N, hipStream_t *stream )
{
  gpuErrchk(hipMemcpyAsync(a, a_d,  N, hipMemcpyDeviceToHost,stream[0]));
  //gpuErrchk(hipMemcpy(a, a_d,  N, hipMemcpyDeviceToHost));
   return;
}

extern "C" void cuda_cpy_dtoh_blocking(void *a_d, void *a, size_t N)
{
  gpuErrchk(hipMemcpy(a, a_d,  N, hipMemcpyDeviceToHost));
   return;
}


/* extern "C" void cuda_cpy_double_htod(double *a, double *a_d, int N)
{
  gpuErrchk(hipMemcpy(a_d, a, sizeof(double) * N, hipMemcpyHostToDevice));
   return;
} */
/* 
extern "C" void cuda_cpy_bool_htod(bool *a, double *a_d, int N)
{
  gpuErrchk(hipMemcpy(a_d, a, sizeof(bool) * N, hipMemcpyHostToDevice));
   //gpuErrchk(hipMemcpyAsync(a_d, a, sizeof(bool) * N, hipMemcpyHostToDevice));
   return;
} */
/* 
extern "C" void cuda_cpy_bool_dtoh(bool *a_d, bool *a, int N)
{
  gpuErrchk(hipMemcpy(a, a_d, sizeof(bool) * N, hipMemcpyDeviceToHost));
   return;
} */

/* extern "C" void cuda_cpy_double_complex_htod(hipDoubleComplex *a, hipDoubleComplex *a_d, int N)
{
   
   gpuErrchk(hipMemcpy(a_d, a, sizeof(hipDoubleComplex) * N, hipMemcpyHostToDevice));
   //gpuErrchk(hipMemcpyAsync(a_d, a, sizeof(hipDoubleComplex) * N, hipMemcpyHostToDevice));
   return;
} */


/* 
extern "C" void cuda_cpy_double_complex_dtoh(hipDoubleComplex *a_d, hipDoubleComplex *a ,int N)
{
  //hipMemcpyAsync( a, a_d, sizeof(double) * N, hipMemcpyDeviceToHost );
  gpuErrchk(hipMemcpy( a, a_d, sizeof(hipDoubleComplex) * N, hipMemcpyDeviceToHost ));
   //printf("\nTest cpy D to H \n");
   
   return;
} */

/* extern "C" void cuda_cpy_int_htod(int *a, int *a_d, int N)
{

   gpuErrchk(hipMemcpy(a_d, a, sizeof(int) * N, hipMemcpyHostToDevice ));
   //gpuErrchk(hipMemcpyAsync(a_d, a, sizeof(int) * N, hipMemcpyHostToDevice ));
   return;
}
 */
/* 
extern "C" void cuda_cpy_double_dtoh(double *a_d, double *a ,int N)
{
  gpuErrchk(hipMemcpy( a, a_d, sizeof(double) * N, hipMemcpyDeviceToHost ));
   
   return;
} */
/* 
extern "C" void cuda_cpy_double_dtod(double *b_d, double *a_d,int N)
{
  gpuErrchk(hipMemcpy( a_d, b_d, sizeof(double) * N, hipMemcpyDeviceToDevice ));
   return;
}
 */

extern "C" void create_cublas_handle(hipblasHandle_t *handle,hipStream_t *stream )
{
 	  hipblasCreate(handle);
    hipStreamCreate(stream);
    hipblasSetStream(*handle, *stream);
//    hipsolverCreate(handle);
//    hipsolverSetStream(*handle,*stream);
    /*printf("\n cublas handle created \n");
    exit(0);*/

   return;
}

extern "C" void destroy_cublas_handle(hipblasHandle_t *handle,hipStream_t *stream )
{
 	 // Destroy the handle
   hipblasDestroy(*handle);
   hipStreamDestroy(*stream);
   //printf("\n cublas handle destroyed. \n The End? \n");
   return;
}

__global__ void gpu_pow(double *a,double *b, double zeta, int N)
{
   int idx = threadIdx.x+blockIdx.x*blockDim.x;
   if (idx<N){
   double loca=a[idx];
   b[idx]=pow(loca,zeta);

 }
}

extern "C" void gpu_kernels_pow(double *a,double *b, double zeta, int size, hipStream_t *stream )
{
  int ntpb=256;
  int nblocks=(size+ntpb-1)/ntpb;
  gpu_pow<<<nblocks, ntpb,0, stream[0]>>>(a,b,zeta, size);
  // gpuErrchk( hipPeekAtLastError() );
  // gpuErrchk( hipDeviceSynchronize() );
  return;
}

extern "C" void gpu_blas_mmul_t_n(hipblasHandle_t handle, const double *Qs_d, const double *soap_d, double *kernels_d, const int n_sparse, const int n_soap, const int n_sites)
//                                                           const double *A,     const double *B,         double *C,       const int nAx,
// const int nAy,      const int nBy,double *b, double zeta, int N)
{
// (hipblasHandle_t handle, const double *Qs_d, const double *soap_d, double *kernels_d, const int n_sparse, const int n_soap, const int n_sites,double *b, double zeta, int N)
	const double alf = 1;
	const double bet = 0;

// soap(n_soap,n_sites)
// Qs(1:n_soap, 1:n_sparse)
// kernels(1:n_sites, 1:n_sparse)
// call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, kernels, n_sites)

	// Do the actual multiplication
  hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n_sites, n_sparse, n_soap, &alf, soap_d, n_soap, Qs_d, n_soap, &bet, kernels_d, n_sites);

  return;
}

extern "C" void gpu_blas_mvmul_n(hipblasHandle_t handle,  double *kernels_copy_d, const double *alphas_d, double *energies_d, const int n_sites, const int n_sparse)
{

	const double alf = 1;
	const double bet = 0;
	const double *alpha = &alf;
	const double *beta = &bet;

	// Do the actual multiplication
  hipblasDgemv(handle, HIPBLAS_OP_N, n_sites,n_sparse, alpha, kernels_copy_d, n_sites, alphas_d, 1, beta, energies_d, 1);
 return;
}



__global__ void gpu_simpleaxpc(double *a, double dccc, double e0, int N)
{
   int idx = threadIdx.x+blockIdx.x*blockDim.x;
   if (idx<N){
   double loca=a[idx];
   a[idx]=dccc*loca+e0;
 }
}

extern "C" void gpu_axpc(double *a, double dccc, double e0, int size, hipStream_t *stream )
{
  int ntpb=256;
  int nblocks=(size+ntpb-1)/ntpb;
  gpu_simpleaxpc<<<nblocks, ntpb,0, stream[0]>>>(a,dccc,e0, size);
  /*gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() );*/
  return;

}

// extern "C" void wrappers_all(double *soap, double *kernels, double *kernels_copy, double *Qs, double *energies, double delta, double zeta, double e0, int n_sites, int n_soap, int n_sparse, int size_kernels, int size_soap, int size_Qs, int size_alphas, int  size_energies)
// {
//   int ntpb=256;
//   int nblocks=(size_kernels+ntpb-1)/ntpb;
//   // Create a handle for CUBLAS
// 	hipblasHandle_t handle;
// 	hipblasCreate(&handle);
//   double *kernels_d, *kernels_copy_d, *soap_d, *Qs_d, *energies_d;
//   hipMalloc( &kernels_d, sizeof(double) * size_kernels );
//   hipMalloc( &kernels_copy_d, sizeof(double) * size_kernels );
//   hipMalloc( &soap_d, sizeof(double) * size_soap );
//   hipMalloc( &Qs_d, sizeof(double) * size_Qs );
//   hipMalloc( &energies_d, sizeof(double)*size_energies);


//   const double alf = 1;
//   const double bet = 0;

//   hipMemcpy(kernels_d, kernels, sizeof(double) * size_kernels, hipMemcpyHostToDevice );
//   hipMemcpy(soap_d, soap, sizeof(double) * size_soap, hipMemcpyHostToDevice );
//   hipMemcpy(Qs_d, Qs, sizeof(double) * size_Qs, hipMemcpyHostToDevice );
//   // Do the actual multiplication

//   hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n_sites, n_sparse, n_soap, &alf, soap_d, n_soap, Qs_d, n_soap, &bet, kernels_d, n_sites);
// //hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N,  nBy, nAx, nAy, alpha, B, nAy, A, nAy, beta, C, nBy);
//     //printf("\n hipblasDgemm \n");
//   // gpu_blas_mmul_t_n(cubhandle,     A,     B,      C,         nAx,      nAy,       nBy,             bb, zeta, N)
//   // gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites, kernels_copy_d, zeta, size_kernels)

//   hipMemcpy( kernels, kernels_d, sizeof(double) * size_kernels, hipMemcpyDeviceToHost );
//   hipLaunchKernelGGL(gpu_pow, dim3(nblocks,1,1), dim3(ntpb,1,1), 0, 0, kernels_d,kernels_copy_d, zeta, size_kernels);
//   hipMemcpy( kernels_copy, kernels_copy_d, sizeof(double) * size_kernels, hipMemcpyDeviceToHost );
// 	// Destroy the handle
// 	hipblasDestroy(handle);
//   hipFree(kernels_d);
//   hipFree(kernels_copy_d);
//   hipFree(soap_d);
//   hipFree(Qs_d);
//   hipFree(energies_d);
//   //printf("\n %d %d %d %d %d %d %d %d  \n", n_sites, n_soap, n_sparse, size_kernels,  size_soap,  size_Qs,  size_alphas,  size_energies);
//   //printf("\n %d %d %d\n", nblocks,ntpb, size_kernels);
//   //exit(0);
//  return;
// }

extern "C" void cuda_set_device( int my_rank)
{

  int  num_gpus=0;
  int mygpuid;
  gpuErrchk(hipGetDeviceCount(&num_gpus));
  gpuErrchk(hipSetDevice(my_rank%num_gpus));
  gpuErrchk(hipGetDevice(&mygpuid));
  /*gpuErrchk(hipSetDevice(0));*/
  printf("\n Seta Aset at %d %d %d %d\n", num_gpus, my_rank%num_gpus,my_rank, mygpuid);
  //exit(0);
  return;
}


__global__ void matvect_kernels(double *kernels_d, double *alphas_d,int  n_sites, int n_sparse)
{
   int idx = threadIdx.x+blockIdx.x*blockDim.x;
   int ispa=idx/n_sites;
   int isite=idx%n_sites;
   if (ispa<n_sparse && isite<n_sites){
     double lock=kernels_d[idx]*alphas_d[ispa];
     kernels_d[idx]=lock;
 }
}

extern "C" void cuda_matvect_kernels(double *kernels_d, double *alphas_d,int  n_sites, int n_sparse, hipStream_t *stream )
{
  int  ntpb=256;
  int nblocks=(n_sites*n_sparse+ntpb-1)/ntpb;
  matvect_kernels<<<nblocks, ntpb,0, stream[0]>>>(kernels_d,alphas_d,n_sites,n_sparse);
  /*gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() );*/
  return;
}



__global__ void matvect_qs(double *qs_d,double *qs_copy_d, double *alphas_d,int  n_soap, int n_sparse)
{
   int idx = threadIdx.x+blockIdx.x*blockDim.x;
   int ispa=idx/n_soap;
   int isoap=idx%n_soap;
   if (ispa<n_sparse && isoap<n_soap){
     double lock=qs_d[idx]*alphas_d[ispa];
     qs_copy_d[idx]=lock;
 }
}

extern "C" void cuda_matvect_qs(double *qs_d,double *qs_copy_d, double *alphas_d,int  n_soap, int n_sparse, hipStream_t *stream )
{
  /*
  alphas(n_sparse)
  allocate( Qs_copy(1:n_soap, 1:n_sparse) )
  do i = 1, n_soap
    Qs_copy(i,:) = Qs(i,:)*alphas(:)
  end do
  */
  int  ntpb=256;
  int nblocks=(n_soap*n_sparse+ntpb-1)/ntpb;
  matvect_qs<<<nblocks, ntpb,0 , stream[0]>>>(qs_d,qs_copy_d,alphas_d,n_soap,n_sparse);
  /*gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() );*/
  return;
}


// gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, n_soap, n_sites, cdelta)
extern "C" void gpu_blas_mmul_n_t(hipblasHandle_t handle, const double *kernels_der_d, const double *Qs_copy_d, 
                       double *Qss_d, const int n_sparse, const int n_soap, const int n_sites, double cdelta)
{

	const double alf = cdelta;
	const double bet = 0;
	const double *alpha = &alf;
	const double *beta = &bet;
// soap(n_soap,n_sites)
// Qs(1:n_soap, 1:n_sparse)
// kernels(1:n_sites, 1:n_sparse)
// call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, kernels, n_sites)
// hipblasDgemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n_sites, n_sparse, n_soap, alpha, soap_d, n_soap, Qs_d, n_soap, beta, kernels_d, n_sites);

// allocate( kernels_der(1:n_sites, 1:n_sparse)
// allocate( Qs_copy(1:n_soap, 1:n_sparse) ))
// allocate( Qss(1:n_sites, 1:n_soap) )
// call dgemm("n", "t", n_sites, n_soap, n_sparse, cdelta, kernels_der, n_sites, Qs_copy, n_soap, 0.d0, Qss, n_sites)
  hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_T, n_sites,  n_soap, n_sparse, alpha,  kernels_der_d, n_sites,
                         Qs_copy_d, n_soap, beta, Qss_d, n_sites);
}


extern "C" void gpu_dgemm_n_n(int m, int n, int k, double alpha,
			      double* A, int lda, double* B, int ldb,
			      double beta, double* C, int ldc, hipblasHandle_t handle){

	const double *alpha_a = &alpha;
	const double *beta_a = &beta;

	hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, m, n, k, alpha_a, A, lda, B, ldb, beta_a, C, ldc);
}




__global__ void  cuda_soap_forces_virial_two(int n_sites,
                                        double *Qss_d,int n_soap,
                                        int *l_index_d, int *j2_index_d,
                                        double3 *soap_der_d,
                                        double3 *xyz_d, double *virial_d,
                                        int n_sites0, double *forces_d)
{
  int l_nn=blockIdx.x;
  int tid=threadIdx.x;
  int i_site=l_index_d[l_nn]-1;
  
  __shared__ double shxthis_block_force[tpb];
  __shared__ double shythis_block_force[tpb];
  __shared__ double shzthis_block_force[tpb];
  
  shxthis_block_force[tid]=0;
  shythis_block_force[tid]=0;
  shzthis_block_force[tid]=0;
  
  double locx_this_force=0;
  double locy_this_force=0;
  double locz_this_force=0;

  for(int ii=tid; ii < n_soap; ii=ii+tpb)
  {
    int i_Qss=i_site+ii*n_sites; // --> (i, 1:n_soap) 
    double loc_this_Qss=Qss_d[i_Qss];// this read  seems OK
    int in_soap_der=(l_nn*n_soap+ii); // (k,:,l) l is pair index, soap_der(3,n_soap,n_pairs)
    double3 loc_soap_der=soap_der_d[in_soap_der];
/*     if(isnan( loc_soap_der.x)|| isnan( loc_soap_der.y)||isnan( loc_soap_der.z)){
      printf("\n loc_soap_der is nan\n");
    } */
    // if(isnan( loc_this_Qss)){
    //   printf("\n loc_this_Qss is nan  %lf %lf %lf %lf\n", loc_this_Qss, loc_soap_der.x, loc_soap_der.y, loc_soap_der.z);
    // }
    locx_this_force+=loc_this_Qss*loc_soap_der.x;
    locy_this_force+=loc_this_Qss*loc_soap_der.y;
    locz_this_force+=loc_this_Qss*loc_soap_der.z;
  }
  
  shxthis_block_force[tid]=locx_this_force;
  shythis_block_force[tid]=locy_this_force;
  shzthis_block_force[tid]=locz_this_force;
/*   if(isnan(locx_this_force)||isnan(locy_this_force)||isnan(locz_this_force)){
    printf("\n loc_this_force is nan\n");
  }   */
  __syncthreads();
  //reduction
  for (int s=tpb/2; s>0; s>>=1) // s=s/2'
  {
    if (tid < s)
    {
      shxthis_block_force[tid] +=shxthis_block_force[tid + s];
      shythis_block_force[tid] +=shythis_block_force[tid + s];
      shzthis_block_force[tid] +=shzthis_block_force[tid + s];
    }
    __syncthreads();
  }
  
  //  at this point this_force is computed
  if(tid==0)
  {
    int j2=j2_index_d[l_nn]-1;
    atomicAdd(&forces_d[j2*3]  , shxthis_block_force[0]);
    atomicAdd(&forces_d[j2*3+1], shythis_block_force[0]);
    atomicAdd(&forces_d[j2*3+2], shzthis_block_force[0]);
    
    // now the virial
    double this_force[3];
    this_force[0]=shxthis_block_force[0];
    this_force[1]=shythis_block_force[0];
    this_force[2]=shzthis_block_force[0];
/*     if(isnan(shxthis_block_force[0])||isnan(shythis_block_force[0])||isnan(shzthis_block_force[0])){
      printf("\n this_force is nan\n");
    } */
    /* if(isnan(this_force[0])||isnan(this_force[1])||isnan(this_force[2])){
      printf("\n this_force is nan\n");
    } */
    
    double3 tmp_xyz;
    tmp_xyz=xyz_d[l_nn];
    double this_xyz[3];
    this_xyz[0]=tmp_xyz.x;
    this_xyz[1]=tmp_xyz.y;
    this_xyz[2]=tmp_xyz.z;
    
    for(int k1=0;k1<3;k1++){
      for(int k2=0;k2<3;k2++){
        double loc_viri=0.5*(this_force[k1]*this_xyz[k2]+this_force[k2]*this_xyz[k1]); 
/*         if(isnan(loc_viri)){
          printf("\n locviri is nan\n");
        } */
        atomicAdd(&virial_d[k2+3*k1], loc_viri);
      }
    }
/*     if(isnan(tmp_xyz.x)||isnan(tmp_xyz.y)||isnan(tmp_xyz.z)){
      printf("\n tmp is nan\n");
    } */
  }
}
 
extern "C" void gpu_final_soap_forces_virial(int n_sites,
                                             double *Qss_d,int n_soap, int *l_index_d, int *j2_index_d,
                                             double3 *soap_der_d,
                                             double3 *xyz_d, double *virial_d,
                                             int n_sites0, 
                                             double *forces_d, int n_pairs, hipStream_t *stream )
{
  dim3 nblocks(n_pairs,1);

  /*double *this_force_d; 
  hipMalloc((void**)&this_force_d,sizeof(double)*n_pairs*3);*/
  hipMemsetAsync(forces_d,0, 3*n_sites0*sizeof(double), stream[0]);
  hipMemsetAsync(virial_d,0, 9*sizeof(double), stream[0]);
     
  cuda_soap_forces_virial_two<<< nblocks, tpb,0, stream[0]>>>(n_sites,
                                              Qss_d,n_soap, l_index_d, j2_index_d,
                                              soap_der_d, xyz_d, virial_d,
                                              n_sites0, forces_d);

  /*gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() );*/

  return;
}



__global__ void  cuda_local_property_derivatives(int n_sites,
                                        double *Qss_d, int n_soap,
                                        int *l_index_d,
                                        double3 *soap_der_d,
                                        double *local_property_cart_der_d)
{
  int l_nn=blockIdx.x;
  int tid=threadIdx.x;
  int i_site=l_index_d[l_nn]-1;
  
  __shared__ double shxthis_block_local_property_cart_der[tpb];
  __shared__ double shythis_block_local_property_cart_der[tpb];
  __shared__ double shzthis_block_local_property_cart_der[tpb];
  
  shxthis_block_local_property_cart_der[tid]=0;
  shythis_block_local_property_cart_der[tid]=0;
  shzthis_block_local_property_cart_der[tid]=0;
  
  double locx_this_local_property_cart_der=0;
  double locy_this_local_property_cart_der=0;
  double locz_this_local_property_cart_der=0;

  // Here, every thread acts on a different pair
  // > We sum over the thread id as long as its less than the n_soap
  // > Then we increment by the number of threads in the block
  // > So, i_site is the pair index which comes from blockidx, this is over n_pairs
  // > Each thread goes over n_soap here, such that we can compute
  //   the dot_product( this_Qss(i,1:n_soap), soap_der( cart, 1:n_soap, n_pairs ))  
  for(int ii=tid; ii < n_soap; ii=ii+tpb)
  {
    int i_Qss=i_site+ii*n_sites; // --> (i, 1:n_soap) 
    double loc_this_Qss=Qss_d[i_Qss];
    int in_soap_der=(l_nn*n_soap+ii); // (k,:,l) l is pair index, soap_der(3,n_soap,n_pairs)
    double3 loc_soap_der=soap_der_d[in_soap_der];

    locx_this_local_property_cart_der+= loc_this_Qss*loc_soap_der.x;
    locy_this_local_property_cart_der+= loc_this_Qss*loc_soap_der.y;
    locz_this_local_property_cart_der+= loc_this_Qss*loc_soap_der.z;
  }
  
  shxthis_block_local_property_cart_der[tid]=locx_this_local_property_cart_der;
  shythis_block_local_property_cart_der[tid]=locy_this_local_property_cart_der;
  shzthis_block_local_property_cart_der[tid]=locz_this_local_property_cart_der;

  __syncthreads();

  // Reduction for the dot product 
  for (int s=tpb/2; s>0; s>>=1) // s=s/2'
    {
      if (tid < s)
	{
	  shxthis_block_local_property_cart_der[tid] +=shxthis_block_local_property_cart_der[tid + s];
	  shythis_block_local_property_cart_der[tid] +=shythis_block_local_property_cart_der[tid + s];
	  shzthis_block_local_property_cart_der[tid] +=shzthis_block_local_property_cart_der[tid + s];
	}
      __syncthreads();
    }

  if(tid == 0)
    {
      // No need for this to be atomic as they should be independent
      local_property_cart_der_d[ 3 * l_nn     ] = - shxthis_block_local_property_cart_der[0];
      local_property_cart_der_d[ 3 * l_nn + 1 ] = - shythis_block_local_property_cart_der[0];
      local_property_cart_der_d[ 3 * l_nn + 2 ] = - shzthis_block_local_property_cart_der[0];      
    }
 }



extern "C" void gpu_local_property_derivatives(int n_sites,
                                             double *Qss_d,int n_soap, int *l_index_d,
                                             double3 *soap_der_d, 
                                             double *local_property_cart_der_d, int n_pairs, hipStream_t *stream )
{
  dim3 nblocks(n_pairs,1);

  /*double *this_force_d; 
  hipMalloc((void**)&this_force_d,sizeof(double)*n_pairs*3);*/
  hipMemsetAsync(local_property_cart_der_d, 0, 3*n_pairs*sizeof(double), stream[0]);
     
  cuda_local_property_derivatives<<< nblocks, tpb, 0, stream[0]>>>(n_sites,
                                              Qss_d, n_soap, l_index_d, 
                                              soap_der_d,
                                              local_property_cart_der_d);

  /*gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() );*/

  return;
}



__global__ void cuda_get_soap_p(double *soap_d, double *sqrt_dot_p_d, double *multiplicity_array_d, 
                           hipDoubleComplex *cnk_d, bool *skip_soap_component_d,
                           int n_sites, int n_soap, int n_max, int l_max)
{
   int i_site = threadIdx.x+blockIdx.x*blockDim.x;
   int k_max=1+l_max*(l_max+1)/2+l_max;
   double my_sqrt_dot_p=0.0;
   if (i_site<n_sites){ 
    int counter=0;
    int counter2=0; 
    //int ssc_counter=0;
    for(int n=1;n<=n_max;n++){
      for(int np=n;np<=n_max;np++){
        for(int l=0;l<=l_max;l++){
          //if(!skip_soap_component_d[ssc_counter]){ //if( skip_soap_component(l, np, n) )cycle
          bool my_skip=skip_soap_component_d[l+(l_max+1)*(np-1+(n-1)*n_max)];
          if(!(my_skip)){ //if( skip_soap_component(l, np, n) )cycle
            
            counter++;
            double my_soap=0.0;//soap_d[counter-1+i_site*n_soap];
            for(int m=0;m<=l; m++){
              int k=1+l*(l+1)/2+m; //k = 1 + l*(l+1)/2 + m
              counter2++;
              hipDoubleComplex tmp_1_cnk_d=cnk_d[i_site+n_sites*((k-1)+(n-1)*k_max)];  //cnk_d[k-1+k_max*(n-1 +i_site*n_max)];
              hipDoubleComplex tmp_2_cnk_d=cnk_d[i_site+n_sites*((k-1)+(np-1)*k_max)]; //cnk_d[k-1+k_max*(np-1+i_site*n_max)];
              my_soap+=multiplicity_array_d[counter2-1]*(tmp_1_cnk_d.x*tmp_2_cnk_d.x+tmp_1_cnk_d.y*tmp_2_cnk_d.y); 
/*               if(isnan(my_soap)){
                printf("\n my_soap is nan %lf %lf %lf %lf %lf!!\n", multiplicity_array_d[counter2-1], tmp_1_cnk_d.x, tmp_1_cnk_d.y,tmp_2_cnk_d.x,tmp_2_cnk_d.y);
              } */
              //soap(counter, i) = soap(counter, i) + multiplicity * real(cnk(k, n, i) * conjg(cnk(k, np, i)))
            }
            soap_d[counter-1+i_site*n_soap]=my_soap;
            my_sqrt_dot_p+=my_soap*my_soap;
          }
        }
      }
    }
    my_sqrt_dot_p=sqrt(my_sqrt_dot_p);
    if(my_sqrt_dot_p<1.0e-5){
      my_sqrt_dot_p=1.0;
    }
    sqrt_dot_p_d[i_site]=my_sqrt_dot_p;
 }
}

extern "C" void gpu_get_sqrt_dot_p(double *sqrt_dot_d, double *soap_d, double *multiplicity_array_d, 
                                   hipDoubleComplex *cnk_d, bool *skip_soap_component_d, 
                                   int n_sites, int n_soap, int n_max, int l_max, hipStream_t *stream )
{
  dim3 nblocks=dim3((n_sites-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);
  cuda_get_soap_p<<< nblocks, nthreads,0 , stream[0]>>>(soap_d,sqrt_dot_d, multiplicity_array_d, cnk_d, skip_soap_component_d, 
                                         n_sites, n_soap, n_max, l_max);                                    
  return;
}



__global__ void cuda_get_soap_der_one(double *soap_rad_der_d, double *soap_azi_der_d, double *soap_pol_der_d, 
                                      double *multiplicity_array_d, 
                                      double *trans_soap_rad_der_d, double *trans_soap_azi_der_d,double *trans_soap_pol_der_d,
                                      hipDoubleComplex *cnk_d, 
                                      hipDoubleComplex *cnk_rad_der_d, hipDoubleComplex *cnk_azi_der_d, hipDoubleComplex *cnk_pol_der_d,
                                      int *k2_i_site_d, bool *skip_soap_component_d, 
                                      int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max)
{
  int k2 = threadIdx.x+blockIdx.x*blockDim.x;
  if (k2<n_atom_pairs){
    int i_site=k2_i_site_d[k2]-1;
    int counter=0;
    int counter2=0; 
    for(int n=1;n<=n_max;n++){
      for(int np=n;np<=n_max;np++){
        for(int l=0;l<=l_max;l++){
          if(!skip_soap_component_d[l+(l_max+1)*(np-1+(n-1)*n_max)]){ //if( skip_soap_component(l, np, n) )cycle // if it happens lots of time, do it in reverse
            counter++;
            double my_soap_rad_der=0; //trans_soap_rad_der_d[k2+(counter-1)*n_atom_pairs]; //soap_rad_der_d[counter-1+k2*n_soap];
            double my_soap_azi_der=0; //trans_soap_azi_der_d[k2+(counter-1)*n_atom_pairs]; //soap_azi_der_d[counter-1+k2*n_soap];
            double my_soap_pol_der=0; //trans_soap_pol_der_d[k2+(counter-1)*n_atom_pairs]; //soap_pol_der_d[counter-1+k2*n_soap];
            for(int m=0;m<=l; m++){
              int k=1+l*(l+1)/2+m; 
              counter2++;
              /* if(threadIdx.x==121 && blockIdx.x==154){
                printf("\n Pair  %d \n" , k2, i_site);              
              } */
              hipDoubleComplex tmp_1_cnk_d=cnk_d[i_site+n_sites*(k-1+(n-1)*k_max)]; //trans_cnk_d[i_site+n_sites*(k-1+(n-1)*k_max)];  //cnk_d[k-1+ k_max*(n-1 +i_site*n_max)];
              hipDoubleComplex tmp_2_cnk_d=cnk_d[i_site+n_sites*(k-1+(np-1)*k_max)]; //trans_cnk_d[i_site+n_sites*(k-1+(np-1)*k_max)]; //cnk_d[k-1+k_max*(np-1+i_site*n_max)];
              hipDoubleComplex tmp_1_cnk_rad_d=cnk_rad_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //trans_cnk_rad_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; // cnk_rad_der_d[k-1+k_max*(n-1 +k2*n_max)];
              hipDoubleComplex tmp_2_cnk_rad_d=cnk_rad_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //trans_cnk_rad_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; // cnk_rad_der_d[k-1+k_max*(np-1+k2*n_max)];
              hipDoubleComplex tmp_1_cnk_azi_d=cnk_azi_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //trans_cnk_azi_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //cnk_azi_der_d[k-1+k_max*(n-1 +k2*n_max)];
              hipDoubleComplex tmp_2_cnk_azi_d=cnk_azi_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //trans_cnk_azi_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //cnk_azi_der_d[k-1+k_max*(np-1+k2*n_max)];
              hipDoubleComplex tmp_1_cnk_pol_d=cnk_pol_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //trans_cnk_pol_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //cnk_pol_der_d[k-1+k_max*(n-1 +k2*n_max)];
              hipDoubleComplex tmp_2_cnk_pol_d=cnk_pol_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //trans_cnk_pol_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //cnk_pol_der_d[k-1+k_max*(np-1+k2*n_max)];
              my_soap_rad_der+=multiplicity_array_d[counter2-1]*(tmp_1_cnk_rad_d.x*tmp_2_cnk_d.x+tmp_1_cnk_rad_d.y*tmp_2_cnk_d.y+
                                                                 tmp_1_cnk_d.x*tmp_2_cnk_rad_d.x+tmp_1_cnk_d.y*tmp_2_cnk_rad_d.y);
              my_soap_azi_der+=multiplicity_array_d[counter2-1]*(tmp_1_cnk_azi_d.x*tmp_2_cnk_d.x+tmp_1_cnk_azi_d.y*tmp_2_cnk_d.y+
                                                                 tmp_1_cnk_d.x*tmp_2_cnk_azi_d.x+tmp_1_cnk_d.y*tmp_2_cnk_azi_d.y);
              my_soap_pol_der+=multiplicity_array_d[counter2-1]*(tmp_1_cnk_pol_d.x*tmp_2_cnk_d.x+tmp_1_cnk_pol_d.y*tmp_2_cnk_d.y+
                                                                 tmp_1_cnk_d.x*tmp_2_cnk_pol_d.x+tmp_1_cnk_d.y*tmp_2_cnk_pol_d.y);
            }   
            trans_soap_rad_der_d[k2+(counter-1)*n_atom_pairs]=my_soap_rad_der; //soap_rad_der_d[counter-1+k2*n_soap]=my_soap_rad_der;
            trans_soap_azi_der_d[k2+(counter-1)*n_atom_pairs]=my_soap_azi_der; //soap_azi_der_d[counter-1+k2*n_soap]=my_soap_azi_der;
            trans_soap_pol_der_d[k2+(counter-1)*n_atom_pairs]=my_soap_pol_der; //soap_pol_der_d[counter-1+k2*n_soap]=my_soap_pol_der;       
          }
        }
      }
    }
  }
}


 __global__ void cuda_get_soap_der_two_one(double *soap_d, double *sqrt_dot_p_d,
                                      double *soap_rad_der_d, double *soap_azi_der_d, double *soap_pol_der_d,
                                      double *trans_soap_rad_der_d, double *trans_soap_azi_der_d, double *trans_soap_pol_der_d,
                                      double *tdotoprod_der_rad, double *tdotoprod_der_azi, double *tdotoprod_der_pol,
                                      int *k2_i_site_d, 
                                      int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max)
{ 
  int k2=blockIdx.x;
  int tid=threadIdx.x;
  int i_site=k2_i_site_d[k2]-1;
  __shared__ double sh_soap_rad_der_dot[tpb];
  __shared__ double sh_soap_azi_der_dot[tpb];
  __shared__ double sh_soap_pol_der_dot[tpb];
  double this_dotprod_rad=0.0;double this_dotprod_azi=0.0;double this_dotprod_pol=0.0;
  
  for(int s=tid;s<n_soap;s=s+tpb){
    this_dotprod_rad+=soap_d[s+i_site*n_soap]*soap_rad_der_d[s+k2*n_soap];
    this_dotprod_azi+=soap_d[s+i_site*n_soap]*soap_azi_der_d[s+k2*n_soap];
    this_dotprod_pol+=soap_d[s+i_site*n_soap]*soap_pol_der_d[s+k2*n_soap];
  }
  sh_soap_rad_der_dot[tid]=this_dotprod_rad;
  sh_soap_azi_der_dot[tid]=this_dotprod_azi;
  sh_soap_pol_der_dot[tid]=this_dotprod_pol;
  __syncthreads();

  //reduction
  for (int s=tpb/2; s>0; s>>=1) // s=s/2
  {
    if (tid < s)
    {
      sh_soap_rad_der_dot[tid] +=sh_soap_rad_der_dot[tid + s];
      sh_soap_azi_der_dot[tid] +=sh_soap_azi_der_dot[tid + s];
      sh_soap_pol_der_dot[tid] +=sh_soap_pol_der_dot[tid + s];
    }
    __syncthreads();

  }
  for(int s=tid;s<n_soap;s=s+tpb){
    tdotoprod_der_rad[s+k2*n_soap]=sh_soap_rad_der_dot[0];
    tdotoprod_der_azi[s+k2*n_soap]=sh_soap_azi_der_dot[0];
    tdotoprod_der_pol[s+k2*n_soap]=sh_soap_pol_der_dot[0];
  }
} 

__global__ void cuda_get_soap_der_two_two(double *soap_d, double *sqrt_dot_p_d,
                                          double *soap_rad_der_d, double *soap_azi_der_d, double *soap_pol_der_d,
                                          double *tdotoprod_der_rad, double *tdotoprod_der_azi, double *tdotoprod_der_pol,
                                          int *k2_i_site_d, 
                                          int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max)
{ 
  int k2=blockIdx.x;
  int tid=threadIdx.x;
  int i_site=k2_i_site_d[k2]-1;
  double loc_sqrt_dot_p=sqrt_dot_p_d[i_site];
  for(int s=tid;s<n_soap;s=s+tpb){
    double my_soap=soap_d[s+i_site*n_soap];

    double my_soap_rad_der=soap_rad_der_d[s+k2*n_soap];
    double my_soap_azi_der=soap_azi_der_d[s+k2*n_soap];
    double my_soap_pol_der=soap_pol_der_d[s+k2*n_soap];

    double myprod_der_rad=tdotoprod_der_rad[s+k2*n_soap];
    double myprod_der_azi=tdotoprod_der_azi[s+k2*n_soap];
    double myprod_der_pol=tdotoprod_der_pol[s+k2*n_soap];


    soap_rad_der_d[s+k2*n_soap]=my_soap_rad_der/loc_sqrt_dot_p
                               -my_soap/(loc_sqrt_dot_p*loc_sqrt_dot_p*loc_sqrt_dot_p)*myprod_der_rad;
    soap_azi_der_d[s+k2*n_soap]=my_soap_azi_der/loc_sqrt_dot_p
                               -my_soap/(loc_sqrt_dot_p*loc_sqrt_dot_p*loc_sqrt_dot_p)*myprod_der_azi;
    soap_pol_der_d[s+k2*n_soap]=my_soap_pol_der/loc_sqrt_dot_p
                               -my_soap/(loc_sqrt_dot_p*loc_sqrt_dot_p*loc_sqrt_dot_p)*myprod_der_pol;
  }

}



__global__ void cuda_get_soap_der_thr_one(double3 *soap_cart_der_d,
                                          double *soap_rad_der_d, double *soap_azi_der_d, double *soap_pol_der_d,
                                          double *thetas, double *phis, double *rjs,
                                          int *k3_index, 
                                          int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max)
{ 
  int k2=blockIdx.x;
  int tid=threadIdx.x;
  int k3=k3_index[k2]-1;

  double my_theta=thetas[k2]; double my_phi=phis[k2]; double my_rj=rjs[k2];
  for(int s=tid;s<n_soap;s=s+tpb){  
    if(k3!=k2){
      double my_soap_rad_der=soap_rad_der_d[s+k2*n_soap];
      double my_soap_azi_der=soap_azi_der_d[s+k2*n_soap];
      double my_soap_pol_der=soap_pol_der_d[s+k2*n_soap];
      double3 my_soap_cart_der;
      my_soap_cart_der.x=sin(my_theta)*cos(my_phi)*my_soap_rad_der 
                        -cos(my_theta)*cos(my_phi)/my_rj*my_soap_pol_der
                        -sin(my_phi)/my_rj*my_soap_azi_der;
      my_soap_cart_der.y=sin(my_theta)*sin(my_phi)*my_soap_rad_der 
                        -cos(my_theta)*sin(my_phi)/my_rj*my_soap_pol_der
                        +cos(my_phi)/my_rj*my_soap_azi_der;
      my_soap_cart_der.z=cos(my_theta)*my_soap_rad_der 
                        +sin(my_theta)/my_rj*my_soap_pol_der;
      soap_cart_der_d[s+k2*n_soap]=my_soap_cart_der;
    }
  }
}

__global__ void cuda_get_soap_der_thr_two(double3 *soap_cart_der_d,
                                          double *soap_rad_der_d, double *soap_azi_der_d, double *soap_pol_der_d,
                                          double *thetas, double *phis, double *rjs,
                                          int *n_neigh_d, int *i_k2_start_d, int *k2_i_site_d, int *k3_index_d, 
                                          int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max, int maxneigh)
{ 
  int i_site=blockIdx.x;
  int tid=threadIdx.x;
  int my_start=i_k2_start_d[i_site]-1;
  int k3=my_start;
  int my_n_neigh=n_neigh_d[i_site];
  
  for(int s=tid;s<n_soap;s=s+tpb){
    double3 loc_sum;
    loc_sum.x=0,loc_sum.y=0,loc_sum.z=0;
    int k2=my_start+1;
    for(int j=1;j<my_n_neigh; j++){
      double3 my_soap_cart_der=soap_cart_der_d[s+k2*n_soap];
      loc_sum.x-=my_soap_cart_der.x;
      loc_sum.y-=my_soap_cart_der.y;
      loc_sum.z-=my_soap_cart_der.z;
      k2++;
    }
    soap_cart_der_d[s+k3*n_soap]=loc_sum;
  }
}


__global__ void naive_transpose_soap_rad_azi_pol(double *soap_rad_der_d,
                                            double *tran_soap_rad_der_d, 
                                            int n_soap, int n_atom_pairs)
{
  int i_g = threadIdx.x+blockIdx.x*blockDim.x;
  if(i_g<n_soap*n_atom_pairs){
    double loc_soap_rad=soap_rad_der_d[i_g];
    int k2=i_g/n_soap;
    int icount=i_g%n_soap;
    int new_i_g=k2+icount*n_atom_pairs;
    tran_soap_rad_der_d[new_i_g]=loc_soap_rad;
  }

}


__global__ void naive_transpose_cnk_arrays(hipDoubleComplex *C,
                                           hipDoubleComplex *tran_C, 
                                           int k_max, int n_max, int n_sites)
{
  // in Fortran is cnk( 1:k_max, 1:n_max, 1:n_sites) --> (1:n_sites,1:k_max, 1:n_max)
  //       cnk_rad_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
  int i_g = threadIdx.x+blockIdx.x*blockDim.x;
  if(i_g<k_max*n_max*n_sites){
    hipDoubleComplex loc_C=C[i_g];  // i_g=i_k+k_max*(i_n+i_site*n_max)
    int i_k=i_g%k_max;
    int i_z=i_g/k_max;
    int i_n=i_z%n_max;
    int i_site=i_z/n_max;
    int new_i_g=i_site+n_sites*(i_k+i_n*k_max);
    tran_C[new_i_g]=loc_C;
  }
}

extern "C" void gpu_get_soap_der(double *soap_d, double *sqrt_dot_d, double3 *soap_cart_der_d, 
                                 double *soap_rad_der_d, double *soap_azi_der_d, double *soap_pol_der_d, 
                                 double *thetas_d, double *phis_d, double *rjs_d, 
                                 double *multiplicity_array_d,
                                 hipDoubleComplex *cnk_d, 
                                 hipDoubleComplex *cnk_rad_der_d, hipDoubleComplex *cnk_azi_der_d, hipDoubleComplex *cnk_pol_der_d, 
                                 int *n_neigh_d, int *i_k2_start_d, int *k2_i_site_d, int *k3_index_d, bool *skip_soap_component_d, 
                                 int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max, int maxneigh, hipStream_t *stream )
{
  dim3 nblocks=dim3((n_atom_pairs-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  dim3 nblocks_get_soap_der_one=dim3((n_atom_pairs-1+tpb_get_soap_der_one)/tpb_get_soap_der_one,1,1);
  dim3 nthreads_get_soap_der_one=dim3(tpb_get_soap_der_one,1,1);
  //size_t mf, ma;
  //hipMemGetInfo(&mf, &ma);
  //printf("\n free: %zu total: %zu", mf, ma);
  double *tdotoprod_der_rad,*tdotoprod_der_azi,*tdotoprod_der_pol; 
  hipMallocAsync((void**)&tdotoprod_der_rad,sizeof(double)*n_atom_pairs*n_soap,stream[0]);
  hipMallocAsync((void**)&tdotoprod_der_azi,sizeof(double)*n_atom_pairs*n_soap,stream[0]);
  hipMallocAsync((void**)&tdotoprod_der_pol,sizeof(double)*n_atom_pairs*n_soap,stream[0]);
  

   double *trans_soap_rad_der_d, *trans_soap_azi_der_d, *trans_soap_pol_der_d;
  hipMallocAsync((void **)&trans_soap_rad_der_d, sizeof(double)*n_atom_pairs*n_soap,stream[0]);
  hipMallocAsync((void **)&trans_soap_azi_der_d, sizeof(double)*n_atom_pairs*n_soap,stream[0]);
  hipMallocAsync((void **)&trans_soap_pol_der_d, sizeof(double)*n_atom_pairs*n_soap,stream[0]);

                                            
  cuda_get_soap_der_one<<< nblocks_get_soap_der_one, nthreads_get_soap_der_one,0, stream[0]>>>(soap_rad_der_d,soap_azi_der_d, soap_pol_der_d, multiplicity_array_d, 
                                               trans_soap_rad_der_d, trans_soap_azi_der_d, trans_soap_pol_der_d, 
                                               cnk_d, cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d,
                                               k2_i_site_d, skip_soap_component_d, 
                                               n_sites,  n_atom_pairs, n_soap,  k_max, n_max, l_max);
                                           
  
  naive_transpose_soap_rad_azi_pol<<< (n_soap*n_atom_pairs+tpb-1)/tpb, tpb,0, stream[0]>>>(trans_soap_rad_der_d,
                                            soap_rad_der_d, 
                                            n_atom_pairs,n_soap);

  naive_transpose_soap_rad_azi_pol<<< (n_soap*n_atom_pairs+tpb-1)/tpb, tpb,0, stream[0]>>>(trans_soap_azi_der_d,
                                            soap_azi_der_d, 
                                            n_atom_pairs,n_soap);
                                            
  naive_transpose_soap_rad_azi_pol<<<(n_soap*n_atom_pairs+tpb-1)/tpb, tpb,0,  stream[0]>>>(trans_soap_pol_der_d,
                                            soap_pol_der_d, 
                                            n_atom_pairs,n_soap);
                                            
  cuda_get_soap_der_two_one<<<n_atom_pairs, nthreads,0, stream[0]>>>(soap_d,sqrt_dot_d, 
                                               soap_rad_der_d,soap_azi_der_d, soap_pol_der_d,
                                               trans_soap_rad_der_d, trans_soap_azi_der_d, trans_soap_pol_der_d,    
                                               tdotoprod_der_rad, tdotoprod_der_azi, tdotoprod_der_pol,                                            
                                               k2_i_site_d, 
                                               n_sites,  n_atom_pairs, n_soap,  k_max, n_max, l_max);
  

  cuda_get_soap_der_two_two<<<n_atom_pairs, nthreads,0, stream[0]>>>(soap_d, sqrt_dot_d,
                                               soap_rad_der_d,soap_azi_der_d, soap_pol_der_d,
                                               //trans_soap_rad_der_d, trans_soap_azi_der_d, trans_soap_pol_der_d,  
                                               tdotoprod_der_rad, tdotoprod_der_azi, tdotoprod_der_pol,                                               
                                               k2_i_site_d, 
                                               n_sites,  n_atom_pairs, n_soap,  k_max, n_max, l_max);
  
  cuda_get_soap_der_thr_one<<<n_atom_pairs, nthreads,0, stream[0]>>>(soap_cart_der_d,  
                                                         soap_rad_der_d,soap_azi_der_d, soap_pol_der_d, 
                                                         thetas_d, phis_d, rjs_d,
                                                         k3_index_d, 
                                                         n_sites,  n_atom_pairs, n_soap,  k_max, n_max, l_max);
  cuda_get_soap_der_thr_two<<<n_sites, nthreads,0, stream[0]>>>(soap_cart_der_d,  
                                                         soap_rad_der_d,soap_azi_der_d, soap_pol_der_d, 
                                                         thetas_d, phis_d, rjs_d,
                                                         n_neigh_d, i_k2_start_d, k2_i_site_d, k3_index_d, 
                                                         n_sites,  n_atom_pairs, n_soap,  k_max, n_max, l_max, maxneigh);                                                       
  //printf("\n YOLO \n");
  hipFreeAsync(tdotoprod_der_rad,   stream[0]);hipFreeAsync(tdotoprod_der_azi,   stream[0]);hipFreeAsync(tdotoprod_der_pol,   stream[0]);
  hipFreeAsync(trans_soap_rad_der_d,stream[0]);hipFreeAsync(trans_soap_azi_der_d,stream[0]);hipFreeAsync(trans_soap_pol_der_d,stream[0]);
  /* hipFreeAsync(trans_cnk_d,0); */
  /* hipFreeAsync(trans_cnk_rad_der_d,0);hipFreeAsync(trans_cnk_azi_der_d,0);hipFreeAsync(trans_cnk_pol_der_d,0); */
  
  // hipError_t code=hipDeviceSynchronize() ;
  // printf("\n %s \n", hipGetErrorString(code));
  // gpuErrchk( code );
  return;
}

__global__ void cuda_soap_normalize(double *soap_d, double *sqrt_dot_d, int n_soap, int n_sites)
{ 
  int i_site=blockIdx.x;
  int tid=threadIdx.x;
  double  my_sqrt_dot_p=sqrt_dot_d[i_site];
  for(int s=tid;s<n_soap;s=s+tpb){
    double my_soap_final=soap_d[s+i_site*n_soap]/my_sqrt_dot_p;
    soap_d[s+i_site*n_soap]=my_soap_final;
  }
}
extern "C" void gpu_soap_normalize(double *soap_d, double *sqrt_dot_d, int n_soap, int n_sites, hipStream_t *stream ){
  cuda_soap_normalize<<< n_sites, tpb,0,stream[0]>>> (soap_d, sqrt_dot_d, n_soap,n_sites);
}
/* 
__global__ void cuda_get_derivatives(double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d, double *radial_exp_coeff_der_d, 
                                   hipDoubleComplex *angular_exp_coeff_rad_der_d, hipDoubleComplex *angular_exp_coeff_azi_der_d, hipDoubleComplex *angular_exp_coeff_pol_der_d,
                                   hipDoubleComplex *cnk_rad_der_d, hipDoubleComplex *cnk_azi_der_d, hipDoubleComplex *cnk_pol_der_d,
                                   double *rjs_d, 
                                   double rcut_max,
                                   int n_atom_pairs, int n_sites, int n_soap, int k_max, int n_max, int l_max)
{
  int k2 = threadIdx.x+blockIdx.x*blockDim.x;
  double Pi=4.0*acos(-1.0);
  if(k2<n_atom_pairs){
    double my_rjs=rjs_d[k2];
    // printf(" %lf  %d\n", my_rjs, k2);
    if(my_rjs<rcut_max){
      for(int n=1;n<=n_max;n++){
        double my_radial_exp_c;
        my_radial_exp_c=radial_exp_coeff_d[n-1+k2*n_max];
        double my_radial_exp_c_der;
        my_radial_exp_c_der=radial_exp_coeff_der_d[n-1+k2*n_max];
        for(int  l=0;l<=l_max;l++){
          for(int m=0;m<=l;m++){
            int k=1+l*(l+1)/2+m;
            hipDoubleComplex my_cnk_rad_der;hipDoubleComplex my_cnk_azi_der;hipDoubleComplex my_cnk_pol_der;

            hipDoubleComplex my_ang_exp_c;
            my_ang_exp_c=angular_exp_coeff_d[k2+n_atom_pairs*(k-1)]; //angular_exp_coeff_d[k-1+k2*k_max];

            hipDoubleComplex my_ang_exp_c_rad_der;
            my_ang_exp_c_rad_der=angular_exp_coeff_rad_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_rad_der_d[k-1+k2*k_max];
            hipDoubleComplex my_ang_exp_c_azi_der;
            my_ang_exp_c_azi_der=angular_exp_coeff_azi_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_azi_der_d[k-1+k2*k_max];
            hipDoubleComplex my_ang_exp_c_pol_der;
            my_ang_exp_c_pol_der=angular_exp_coeff_pol_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_pol_der_d[k-1+k2*k_max];

            my_cnk_rad_der.x=Pi*(my_ang_exp_c.x*my_radial_exp_c_der+my_ang_exp_c_rad_der.x*my_radial_exp_c);
            my_cnk_rad_der.y=Pi*(my_ang_exp_c.y*my_radial_exp_c_der+my_ang_exp_c_rad_der.y*my_radial_exp_c);

            my_cnk_azi_der.x=Pi*(my_ang_exp_c_azi_der.x*my_radial_exp_c);
            my_cnk_azi_der.y=Pi*(my_ang_exp_c_azi_der.y*my_radial_exp_c);

            my_cnk_pol_der.x=Pi*(my_ang_exp_c_pol_der.x*my_radial_exp_c);
            my_cnk_pol_der.y=Pi*(my_ang_exp_c_pol_der.y*my_radial_exp_c);

            cnk_rad_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_rad_der;
            cnk_azi_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_azi_der;
            cnk_pol_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_pol_der;
          }
        }
      }
    }
  }
}
 */
/* 
__global__ void cuda_get_derivatives_new(double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d, double *radial_exp_coeff_der_d, 
                                   hipDoubleComplex *angular_exp_coeff_rad_der_d, hipDoubleComplex *angular_exp_coeff_azi_der_d, hipDoubleComplex *angular_exp_coeff_pol_der_d,
                                   hipDoubleComplex *cnk_rad_der_d, hipDoubleComplex *cnk_azi_der_d, hipDoubleComplex *cnk_pol_der_d,
                                   double *rjs_d, 
                                   double rcut_max,
                                   int n_atom_pairs, int n_sites, int n_soap, int k_max, int n_max, int l_max)
{
  int k2 =blockIdx.x;
  int k=threadIdx.x+1;
  double Pi=4.0*acos(-1.0);
  double my_rjs=rjs_d[k2];
  // printf(" %lf  %d\n", my_rjs, k2);
  if(my_rjs<rcut_max){
    
    for(int n=1;n<=n_max;n++){
      double my_radial_exp_c;
      my_radial_exp_c=radial_exp_coeff_d[n-1+k2*n_max];
      double my_radial_exp_c_der;
      my_radial_exp_c_der=radial_exp_coeff_der_d[n-1+k2*n_max];
      //int k=1+l*(l+1)/2+m;
      
      hipDoubleComplex my_cnk_rad_der;hipDoubleComplex my_cnk_azi_der;  hipDoubleComplex my_cnk_pol_der;
      
      hipDoubleComplex my_ang_exp_c; my_ang_exp_c=angular_exp_coeff_d[k2+n_atom_pairs*(k-1)]; //angular_exp_coeff_d[k-1+k2*k_max];
      
      hipDoubleComplex my_ang_exp_c_rad_der;  my_ang_exp_c_rad_der=angular_exp_coeff_rad_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_rad_der_d[k-1+k2*k_max];
      
      hipDoubleComplex my_ang_exp_c_azi_der;  my_ang_exp_c_azi_der=angular_exp_coeff_azi_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_azi_der_d[k-1+k2*k_max];
      
      hipDoubleComplex my_ang_exp_c_pol_der;  my_ang_exp_c_pol_der=angular_exp_coeff_pol_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_pol_der_d[k-1+k2*k_max];
      
      my_cnk_rad_der.x=Pi*(my_ang_exp_c.x*my_radial_exp_c_der+my_ang_exp_c_rad_der.x*my_radial_exp_c);
      my_cnk_rad_der.y=Pi*(my_ang_exp_c.y*my_radial_exp_c_der+my_ang_exp_c_rad_der.y*my_radial_exp_c);
      
      my_cnk_azi_der.x=Pi*(my_ang_exp_c_azi_der.x*my_radial_exp_c);
      my_cnk_azi_der.y=Pi*(my_ang_exp_c_azi_der.y*my_radial_exp_c);
      
      my_cnk_pol_der.x=Pi*(my_ang_exp_c_pol_der.x*my_radial_exp_c);
      my_cnk_pol_der.y=Pi*(my_ang_exp_c_pol_der.y*my_radial_exp_c);
      
      cnk_rad_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_rad_der;
      cnk_azi_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_azi_der;
      cnk_pol_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_pol_der;
    }
  }
}
 */

__global__ void cuda_get_derivatives_new_new(double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d, double *radial_exp_coeff_der_d, 
                                   hipDoubleComplex *angular_exp_coeff_rad_der_d, hipDoubleComplex *angular_exp_coeff_azi_der_d, hipDoubleComplex *angular_exp_coeff_pol_der_d,
                                   hipDoubleComplex *cnk_rad_der_d, hipDoubleComplex *cnk_azi_der_d, hipDoubleComplex *cnk_pol_der_d,
                                   double *rjs_d, 
                                   double rcut_max,
                                   int n_atom_pairs,int n_sites, int n_soap, int k_max, int n_max, int l_max)
{
  int k2 =threadIdx.x+blockDim.x*blockIdx.x;
  int n=blockIdx.y+1;
  int k=blockIdx.z+1;
  double Pi=4.0*acos(-1.0);
  if(k2<n_atom_pairs){
    double my_rjs=rjs_d[k2];
    if(my_rjs<rcut_max){
      double my_radial_exp_c;
      my_radial_exp_c=radial_exp_coeff_d[n-1+k2*n_max];
      double my_radial_exp_c_der;
      my_radial_exp_c_der=radial_exp_coeff_der_d[n-1+k2*n_max];
      //int k=1+l*(l+1)/2+m;
      hipDoubleComplex my_cnk_rad_der;hipDoubleComplex my_cnk_azi_der;hipDoubleComplex my_cnk_pol_der;
      hipDoubleComplex my_ang_exp_c; my_ang_exp_c=angular_exp_coeff_d[k2+n_atom_pairs*(k-1)]; //angular_exp_coeff_d[k-1+k2*k_max];
      hipDoubleComplex my_ang_exp_c_rad_der;my_ang_exp_c_rad_der=angular_exp_coeff_rad_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_rad_der_d[k-1+k2*k_max];
      hipDoubleComplex my_ang_exp_c_azi_der;my_ang_exp_c_azi_der=angular_exp_coeff_azi_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_azi_der_d[k-1+k2*k_max];
      hipDoubleComplex my_ang_exp_c_pol_der;my_ang_exp_c_pol_der=angular_exp_coeff_pol_der_d[k2+(k-1)*n_atom_pairs]; //angular_exp_coeff_pol_der_d[k-1+k2*k_max];
      
      my_cnk_rad_der.x=Pi*(my_ang_exp_c.x*my_radial_exp_c_der+my_ang_exp_c_rad_der.x*my_radial_exp_c);
      my_cnk_rad_der.y=Pi*(my_ang_exp_c.y*my_radial_exp_c_der+my_ang_exp_c_rad_der.y*my_radial_exp_c);
      
      my_cnk_azi_der.x=Pi*(my_ang_exp_c_azi_der.x*my_radial_exp_c);
      my_cnk_azi_der.y=Pi*(my_ang_exp_c_azi_der.y*my_radial_exp_c);
      
      my_cnk_pol_der.x=Pi*(my_ang_exp_c_pol_der.x*my_radial_exp_c);
      my_cnk_pol_der.y=Pi*(my_ang_exp_c_pol_der.y*my_radial_exp_c);
      
      //i_site+n_sites*(i_k+i_n*k_max);
      cnk_rad_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max)]=my_cnk_rad_der; //cnk_rad_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_rad_der;
      cnk_azi_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max)]=my_cnk_azi_der; //cnk_azi_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_azi_der;
      cnk_pol_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max)]=my_cnk_pol_der; //cnk_pol_der_d[k-1+k_max*(n-1+k2*n_max)]=my_cnk_pol_der;
    }
  }
}


extern "C" void gpu_get_derivatives(double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d, double *radial_exp_coeff_der_d, 
                                    hipDoubleComplex *angular_exp_coeff_rad_der_d, hipDoubleComplex *angular_exp_coeff_azi_der_d, hipDoubleComplex *angular_exp_coeff_pol_der_d,
                                    hipDoubleComplex *cnk_rad_der_d, hipDoubleComplex *cnk_azi_der_d, hipDoubleComplex *cnk_pol_der_d,
                                    double *rjs_d, double rcut_max,
                                    int n_atom_pairs, int n_sites, int n_soap, int k_max, int n_max, int l_max, hipStream_t *stream )
{
  /*dim3 nblocks=dim3((n_atom_pairs-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  hipEvent_t start, stop;
  hipEventCreate(&start);
  hipEventCreate(&stop);
  float milliseconds;
  hipEventRecord(start);
  for(int lll=1;lll<=1000;lll++){*//*    
  
  hipLaunchKernelGGL(cuda_get_derivatives, nblocks, nthreads, 0, 0, radial_exp_coeff_d, angular_exp_coeff_d, radial_exp_coeff_der_d,
                                              angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d,
                                              cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d,
                                              rjs_d, rcut_max,
                                              n_atom_pairs, n_soap,k_max, n_max, l_max);
  
  *//*}
  hipEventRecord(stop);
  hipEventSynchronize(stop);
   milliseconds = 0.0;
  hipEventElapsedTime(&milliseconds, start, stop);
  printf("\n Time of the first kernel in s %f\n", milliseconds/1000.0);

  

  hipEventRecord(start);
  for(int lll=1;lll<=1000;lll++){
  
  */
  
  cuda_get_derivatives_new_new<<<dim3((n_atom_pairs+tpbcnk-1)/tpbcnk,n_max,k_max), tpbcnk,0,stream[0]>>>(radial_exp_coeff_d, 
                                              angular_exp_coeff_d, radial_exp_coeff_der_d,
                                              angular_exp_coeff_rad_der_d, angular_exp_coeff_azi_der_d, angular_exp_coeff_pol_der_d,
                                              cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d,
                                              rjs_d, rcut_max,
                                              n_atom_pairs, n_sites, n_soap,k_max, n_max, l_max);
  
  /*}
  hipEventRecord(stop);
  hipEventSynchronize(stop);
   milliseconds = 0.0;
  hipEventElapsedTime(&milliseconds, start, stop);
  printf("\n Time of the second kernel in s %f\n", milliseconds/1000.0);

  exit(0);*/
                                              
}

__global__ void cuda_get_cnk_one(hipDoubleComplex *cnk_d, double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d,                             
                                 int *n_neigh_d, int  *k2_start_d, 
                                 int n_atom_pairs, int n_sites, int k_max, int n_max, int l_max)
{
  int i_site=threadIdx.x+blockIdx.x*blockDim.x;
  double pi=4.0*acos(-1.0);
  if(i_site<n_sites)
  {
    int k2=k2_start_d[i_site];
    int my_n_neigh=n_neigh_d[i_site];
    for(int j=1;j<=my_n_neigh; j++)
    {
      k2++;
      for(int n=1; n<=n_max; n++)
      {
        double loc_rad_exp_coeff=radial_exp_coeff_d[n-1+n_max*(k2-1)];
        for(int l=0; l<=l_max; l++)
        {
          for(int m=0; m<=l; m++)
          {
            int k=1+l*(l+1)/2+m; //k=1+
            hipDoubleComplex loc_cnk=cnk_d[i_site+n_sites*((k-1)+(n-1)*k_max)]; //cnk_d[k-1+k_max*(n-1+n_max*i_site)];
            hipDoubleComplex loc_ang_exp_coeff=angular_exp_coeff_d[k2-1+n_atom_pairs*(k-1)];  // angular_exp_coeff_d[k-1+k_max*(k2-1)];
            loc_cnk.x+=pi*loc_rad_exp_coeff*loc_ang_exp_coeff.x;
            loc_cnk.y+=pi*loc_rad_exp_coeff*loc_ang_exp_coeff.y;
            cnk_d[i_site+n_sites*((k-1)+(n-1)*k_max)]=loc_cnk; //cnk_d[k-1+k_max*(n-1+n_max*i_site)]=loc_cnk;
          }
        }
      }
    }
  }
}

__global__ void cuda_get_cnk_one_new_new(hipDoubleComplex *cnk_d, double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d,                             
                                 int *n_neigh_d, int  *k2_start_d, 
                                 int n_atom_pairs, int n_sites, int k_max, int n_max, int l_max)
{
  int i_site=threadIdx.x+blockIdx.x*blockDim.x;
  int n=blockIdx.y+1;
  int k=blockIdx.z+1;
  double pi=4.0*acos(-1.0);
  if(i_site<n_sites){
    int k2=k2_start_d[i_site];
    int my_n_neigh=n_neigh_d[i_site];
    hipDoubleComplex loc_cnk;
    loc_cnk.x=0.0; loc_cnk.y=0.0;
    /*if(k<=k_max){
       cnk_d[k-1+k_max*(n-1+n_max*i_site)]; // coalesced???
    }*/
      for(int j=1;j<=my_n_neigh; j++){
        k2++;
        double loc_rad_exp_coeff=radial_exp_coeff_d[n-1+n_max*(k2-1)]; //coalesced ???
        //int k=1+l*(l+1)/2+m;
        hipDoubleComplex loc_ang_exp_coeff=angular_exp_coeff_d[k2-1+n_atom_pairs*(k-1)];  // angular_exp_coeff_d[k-1+k_max*(k2-1)]; //coalesced ??
        loc_cnk.x+=pi*loc_rad_exp_coeff*loc_ang_exp_coeff.x;
        loc_cnk.y+=pi*loc_rad_exp_coeff*loc_ang_exp_coeff.y;
/*         if(isnan(loc_cnk.x)||isnan(loc_cnk.y)){
          printf("\n loc_cnk is nan %lf %lf %lf %lf %lf %lf", loc_cnk.x,loc_cnk.y, loc_ang_exp_coeff.x, loc_ang_exp_coeff.y,loc_rad_exp_coeff ,pi);
        } *//* 
        if(isnan(loc_ang_exp_coeff.x)||isnan(loc_ang_exp_coeff.y)){
          printf("\n loc_cnk is nan %lf %lf %lf %lf %lf %lf", loc_cnk.x,loc_cnk.y, loc_ang_exp_coeff.x, loc_ang_exp_coeff.y,loc_rad_exp_coeff,pi);
        } */
/*         if(isnan(loc_rad_exp_coeff)){
          printf("\n loc_rad_exp_coeff is nan %lf %lf %lf %lf %lf %lf", loc_cnk.x,loc_cnk.y, loc_ang_exp_coeff.x, loc_ang_exp_coeff.y,loc_rad_exp_coeff,pi);
        } */
      }
      if(k<=k_max){
        cnk_d[i_site+n_sites*((k-1)+(n-1)*k_max)]=loc_cnk; //cnk_d[k-1+k_max*(n-1+n_max*i_site)]=loc_cnk;
      }
  }
}

__global__ void cuda_get_cnk_two(hipDoubleComplex *cnk_d, double *atom_sigma_r, double *atom_sigma_t, double *rcut_hard, double *central_weight,  
                                 int *species, int *i_beg, int *i_end, int radial_enhancement, int *species_multiplicity,
                                 double *W, double *S,
                                 int n_sites, int k_max, int n_max, int size_species_1)
{
  int i_site=threadIdx.x+blockIdx.x*blockDim.x; //if (i_site >= n_sites) return;
  double pi=acos(-1.0);
  
  if(i_site<n_sites)
  {
    for (int k = 1; k <= species_multiplicity[i_site]; k++){
      int j = species[i_site*size_species_1+k-1]-1;
      double amplitude;
      double atom_sigma_r_j=atom_sigma_r[j];
      double atom_sigma_t_j=atom_sigma_t[j];
      double rcut_hard_j=rcut_hard[j];
      double central_weight_j=central_weight[j];
      if (radial_enhancement == 1){
        amplitude = sqrt(2.0/pi) * atom_sigma_r_j / rcut_hard_j;
      } else if (radial_enhancement == 2) {
        amplitude = (atom_sigma_r_j*atom_sigma_r_j) /(rcut_hard_j*rcut_hard_j);
      } else {
        amplitude = 1.0;
      }
      
      int i_beg_j = i_beg[j];
      int i_end_j = i_end[j];
      
      for (int n = i_beg_j; n <= i_end_j; n++) {
        double mmul_WS=0.0;
        for(int d=i_beg_j; d <= i_end_j; d++){
          mmul_WS+=W[n-1+(d-1)*n_max]*S[d-1+(i_end_j-1)*n_max];
          if(isnan(W[n-1+(d-1)*n_max])){
            printf("W is nan %lf\n", W[n-1+(d-1)*n_max]);
          }
          if(isnan(S[d-1+(i_end_j-1)*n_max])){
            printf("S is nan %lf\n",S[d-1+(i_end_j-1)*n_max]);
          } 
          /*if(n>n_max || d>n_max){
            printf("%d %d %d %d\n", i_site, l, d, n_max);
          }*/
        }
        hipDoubleComplex l_cnk=cnk_d[i_site+n_sites*((n-1)*k_max)]; //cnk_d[k_max*(n-1+n_max*i_site)];
        l_cnk.x +=amplitude * central_weight_j * sqrt(4.0*pi)*sqrt(sqrt(pi))*  
                  sqrt(atom_sigma_r_j/2.0)*
                 (rcut_hard_j*rcut_hard_j*rcut_hard_j)/(atom_sigma_t_j*atom_sigma_t_j)/
                  atom_sigma_r_j*mmul_WS;
        cnk_d[i_site+n_sites*((n-1)*k_max)]=l_cnk; //cnk_d[k_max*(n-1+n_max*i_site)]=l_cnk;
      }
    }
  }
}

extern "C" void gpu_get_cnk(double *radial_exp_coeff_d, hipDoubleComplex *angular_exp_coeff_d,
                            hipDoubleComplex *cnk_d, 
                            int *n_neigh_d, int  *k2_start_d,
                            int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max, int l_max,
                            int bintybint,
                            double *atom_sigma_r_d, double *atom_sigma_t_d, double *rcut_hard_d, 
                            double *central_weight_d,  int *species_d, int *i_beg_d, int *i_end_d, 
                            int radial_enhancement, int *species_multiplicity_d,
                            double *W_d, double *S_d, int size_species_1, hipStream_t *stream )
{
  //hipMemsetAsync(cnk_d,0, k_max*n_max*n_sites*sizeof(hipDoubleComplex),stream[0]);
  /*hipEvent_t start, stop;
  hipEventCreate(&start);
  hipEventCreate(&stop);
  float milliseconds;
  hipEventRecord(start);
  for(int lll=1;lll<=1000;lll++){*/
  
  /*dim3 nblocks=dim3((n_sites-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);
  hipLaunchKernelGGL(cuda_get_cnk_one, nblocks, nthreads, 0, 0, cnk_d, radial_exp_coeff_d, angular_exp_coeff_d,
                                         n_neigh_d, k2_start_d,
                                         n_sites, k_max, n_max, l_max);*/
  /*}
  hipEventRecord(stop);
  hipEventSynchronize(stop);
   milliseconds = 0.0;
  hipEventElapsedTime(&milliseconds, start, stop);
  printf("\n Time of the first kernel in s %f\n", milliseconds/1000.0);
  */
  
  dim3 nnth=dim3(tpbcnk, 1,1 ); // each block does the inner loops over l and m in total k_max per block 
  
  /*hipEventRecord(start);

  for(int lll=1;lll<=1000;lll++){*/
  cuda_get_cnk_one_new_new<<<dim3((n_sites+tpbcnk-1)/tpbcnk,n_max,k_max), nnth,0,stream[0]>>>(cnk_d, radial_exp_coeff_d, angular_exp_coeff_d,
                                          n_neigh_d, k2_start_d,
                                          n_atom_pairs, n_sites, k_max, n_max, l_max);                            
  /*}
  hipEventRecord(stop);
  hipEventSynchronize(stop);
  milliseconds = 0.0;
  hipEventElapsedTime(&milliseconds, start, stop);
  printf("\n Time of the second kernel in s %f\n", milliseconds/1000.0);
  exit(0);*/
  if(bintybint==1000){
    dim3 nblocks=dim3((n_sites-1+tpb)/tpb,1,1);
    dim3 nthreads=dim3(tpb,1,1);
    cuda_get_cnk_two<<<nblocks, nthreads,0,stream[0]>>>(cnk_d, atom_sigma_r_d, atom_sigma_t_d, rcut_hard_d, central_weight_d,  
                                           species_d, i_beg_d, i_end_d, radial_enhancement, species_multiplicity_d,
                                           W_d, S_d,
                                           n_sites, k_max, n_max, size_species_1);
  }
}

__global__ void cuda_get_plm_arrays_one(double *plm_array_global_d,int kmax, int lmax, double *thetas_d, int n_atom_pairs){
  int k_ij=threadIdx.x+blockIdx.x*blockDim.x;
  if(k_ij<n_atom_pairs){
    double x=cos(thetas_d[k_ij]);
    // compute the first 6 polynomials to initialize the recursion series
    plm_array_global_d[k_ij]=1.0;                             //plm_array_global_d[  k_ij*kmax]=1.0;
    plm_array_global_d[k_ij+1*n_atom_pairs]=x;                       //plm_array_global_d[1+k_ij*kmax]=x;               
    plm_array_global_d[k_ij+2*n_atom_pairs]=-sqrt(1.0-x*x);          // plm_array_global_d[2+k_ij*kmax]=-sqrt(1.0-x*x);
    plm_array_global_d[k_ij+3*n_atom_pairs]=1.5*x*x-0.5;             //plm_array_global_d[3+k_ij*kmax]=1.5*x*x-0.5;
    plm_array_global_d[k_ij+4*n_atom_pairs]=-3.0*x*sqrt(1.0-x*x);    //plm_array_global_d[4+k_ij*kmax]=-3.0*x*sqrt(1.0-x*x);
    plm_array_global_d[k_ij+5*n_atom_pairs]=3.0 -3.0*x*x;            //plm_array_global_d[5+k_ij*kmax]=3.0 -3.0*x*x;        
    
    for(int l=3;l<=lmax;l++){
      int k=0;
      for(int m=0;m<=l-2; m++){
        k=1+l*(l+1)/2+m;
        plm_array_global_d[k_ij+(k-1)*n_atom_pairs]=((2.0*l-1.0)*x*plm_array_global_d[k_ij+(k-l-1)*n_atom_pairs]-    //plm_array_global_d[k-1+k_ij*kmax]
                                                    (l-1.0+m)*plm_array_global_d[k_ij+(k-2*l+1-1)*n_atom_pairs])/(l-m);
      }
      k=k+1;
      plm_array_global_d[k_ij+(k-1)*n_atom_pairs]=x*(2.0*l-1.0)*plm_array_global_d[k_ij+(k-l-1)*n_atom_pairs];                //plm_array_global_d[k-1+k_ij*kmax]
      k=k+1;
      plm_array_global_d[k_ij+(k-1)*n_atom_pairs]=-(2.0*l-1.0)*sqrt(1.0-x*x)*plm_array_global_d[k_ij+(k-l-1-1)*n_atom_pairs]; //plm_array_global_d[k-1+k_ij*kmax]
    }
  }
}

extern "C" void  gpu_get_plm_array_global(double *plm_array_global_d, int n_atom_pairs, int kmax, 
                                     int lmax, double *thetas_d, hipStream_t *stream )
{
  dim3 nblocks=dim3((n_atom_pairs-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);
  cuda_get_plm_arrays_one<<< nblocks, nthreads,0,stream[0]>>>(plm_array_global_d,kmax,lmax, thetas_d, n_atom_pairs );
}


__global__ void cuda_get_exp_coeff_one(hipDoubleComplex *eimphi_global_d, double *rjs_d, double *phis_d,
                                      bool *mask_d, double *atom_sigma_in_d, double *atom_sigma_scaling_d,
                                      double rcut, int n_atom_pairs, int n_species, int lmax, int kmax,
                                      double *prefl_global_d, hipDoubleComplex *prefm_global_d, double *preflm_d, 
                                      double *plm_array_global_d, hipDoubleComplex *exp_coeff_d) //, double *fact_array_d)
{
  int k_ij=threadIdx.x+blockIdx.x*blockDim.x;
  double xcut = 1.0e-7;
  if(k_ij<n_atom_pairs){
    double rj=rjs_d[k_ij];
    double phi=phis_d[k_ij];
    if(rj<rcut){
      int i_sp=1;
      for(i_sp=1;i_sp<=n_species; i_sp++)
        if(mask_d[k_ij+(i_sp-1)*n_atom_pairs]) break;
      double atom_sigma=atom_sigma_in_d[i_sp-1] + atom_sigma_scaling_d[i_sp-1]*rj;
      double scaling=atom_sigma_scaling_d[i_sp-1];
      double rjbysigma=rj/atom_sigma;
      double amplitude=(rcut*rcut)/(atom_sigma*atom_sigma);
      double x=rjbysigma;
      double x2=x*x;
      double x4=x2*x2;
      double flm2, flm1,fl; 
      double coeff1 = 2.0*rj/atom_sigma*atom_sigma;
      double coeff2 = 1.0 - scaling*rj/atom_sigma;
      
      if(x>0){
        flm2=fabs((1.0-exp(-2.0*x2))/2.0/x2);
        flm1=fabs((x2-1.0+exp(-2.0*x2)*(x2+1.0))/2.0/x4);
      } else {
        flm2=1.0;
        flm1=0.0;
      }
      //  Complex exponential using Euler's formula and Chebyshev recursion
      double cosm2 = cos(phi);
      double cosphi2 = 2.0 * cosm2;
      double sinm2 = -sin(phi);
      double cosm1 = 1.0;
      double sinm1 = 0.0;
      hipDoubleComplex  loc_prefm;
      loc_prefm.x= 1.0;
      loc_prefm.y=0.0; // (1.d0, 0.d0)*
      prefm_global_d[k_ij] = loc_prefm;
      double ilexp=-500000;// no need for this, just makig sure the next lines were working
      
      double fact=1.0;
      int k=0;
      for(int l=0;l<=lmax; l++){
        if(l>0) fact=fact*(2.0*l+1.0);
        if(l==0){
          if(x< xcut) ilexp=1.0-x2;
          else ilexp=flm2;
        }
        else if(l==1){
          if(x2/1000.0<xcut) ilexp=(x2-x4)/fact; //fact_array_d[l-1];
          else ilexp=flm1;
        }
        else{
          if(pow(x2,l)/fact*l<xcut) fl=pow(x2,l)/fact;
          else fl=fabs(flm2-(2.0*l-1.0)/x2*flm1);
          flm2=flm1;
          flm1=fl;
          ilexp=fl;
        }
        if(l>0){
          double cos0 = cosphi2 * cosm1 - cosm2;
          double sin0 = cosphi2 * sinm1 - sinm2;
          cosm2 = cosm1;
          sinm2 = sinm1;
          cosm1 = cos0;
          sinm1 = sin0;
          loc_prefm.x=cos0;
          loc_prefm.y=-sin0;
          prefm_global_d[k_ij+l*n_atom_pairs] = loc_prefm;
        }
        prefl_global_d[k_ij+l*n_atom_pairs]=ilexp;
        for(int m=0;m<=l;m++){
          hipDoubleComplex loc_emphi;
          hipDoubleComplex tmp_prefm=prefm_global_d[k_ij+m*n_atom_pairs];
          loc_emphi.x=ilexp*tmp_prefm.x;
          loc_emphi.y=ilexp*tmp_prefm.y;
          eimphi_global_d[k_ij+k*n_atom_pairs]=loc_emphi;
          hipDoubleComplex loc_exp_coeff;
          loc_exp_coeff.x=amplitude*preflm_d[k]*plm_array_global_d[k_ij+k*n_atom_pairs]*loc_emphi.x;
          loc_exp_coeff.y=amplitude*preflm_d[k]*plm_array_global_d[k_ij+k*n_atom_pairs]*loc_emphi.y;
          exp_coeff_d[k_ij+k*n_atom_pairs]=loc_exp_coeff; 
          //exp_coeff_d[k+k_ij*kmax]=loc_exp_coeff;// naive transpose
          k++;
        }
      }
    }
  }
}


/*__global__ void cuda_get_fact_array(double *fact_array_d, int lmax)
{
  if(lmax>0){
    double fact=1.0;
    for(int l=1;l<=lmax; l++){
      fact=fact*(2.0*l+1.0);
      fact_array_d[l-1]=fact;
      //printf("%lf %d \n", fact, l);
    }
  }
}
*/



__global__ void cuda_get_exp_coeff_der_one(hipDoubleComplex *eimphi_global_d, double *rjs_d, double *phis_d,
                                      bool *mask_d, double *atom_sigma_in_d, double *atom_sigma_scaling_d,
                                      double rcut, int n_atom_pairs, int n_species, int lmax, int kmax,
                                      double *prefl_global_d, hipDoubleComplex *prefm_global_d, double *preflm_d, 
                                      double *prefl_global_der_d, 
                                      double *plm_array_global_d, double *plm_array_global_der_d, hipDoubleComplex *exp_coeff_d,
                                      hipDoubleComplex *eimphi_rad_der_global_d, hipDoubleComplex *eimphi_azi_der_global_d,
                                      double *plm_array_div_sin, double *plm_array_der_mul_sin, 
                                      hipDoubleComplex *exp_coeff_rad_der_d, hipDoubleComplex *exp_coeff_azi_der_d, hipDoubleComplex *exp_coeff_pol_der_d) 
{
  int k_ij=threadIdx.x+blockIdx.x*blockDim.x;

  if(k_ij<n_atom_pairs){
    double rj=rjs_d[k_ij];
    /*double phi=phis_d[k_ij];*/
    if(rj<rcut){
      int i_sp=1;
      for(i_sp=1;i_sp<=n_species; i_sp++){
        if(mask_d[k_ij+(i_sp-1)*n_atom_pairs]){
          break;
        }
      }
      double atom_sigma=atom_sigma_in_d[i_sp-1] + atom_sigma_scaling_d[i_sp-1]*rj;
      double scaling=atom_sigma_scaling_d[i_sp-1];
      double amplitude=(rcut*rcut)/(atom_sigma*atom_sigma);
      /*double rjbysigma=rj/atom_sigma;
      double x=rjbysigma;
      double x2=x*x;
      double x4=x2*x2;/
      double flm2, flm1,fl; */

      double coeff1 = 2.0*rj/(atom_sigma*atom_sigma);
      double coeff2 = 1.0 - scaling*rj/atom_sigma;
      //hipDoubleComplex  loc_prefm;
      //double ilexp=-500000;// no need for this, just makig sure the next lines were working
      

      int k=0;
      
      double ilexp_der;

      for(int l=0;l<=lmax; l++){
        if(l==0){
          ilexp_der=coeff1*(prefl_global_d[k_ij+n_atom_pairs]-prefl_global_d[k_ij]);
        }
        else{
          ilexp_der=(-coeff1-(2.0*l+2.0)/rj)*prefl_global_d[k_ij+l*n_atom_pairs]+coeff1*prefl_global_d[k_ij+(l-1)*n_atom_pairs];
        }
        if(rj<1.0e-5){
          ilexp_der=0.0;
        }
        ilexp_der*=coeff2;
        prefl_global_der_d[k_ij+l*n_atom_pairs]=ilexp_der;
        for(int m=0;m<=l;m++){
          hipDoubleComplex loc_exp_coeff=exp_coeff_d[k_ij+k*n_atom_pairs];
          double loc_preflm=preflm_d[k];
          hipDoubleComplex loc_emphi=eimphi_global_d[k_ij+k*n_atom_pairs];
          hipDoubleComplex tmp_prefm=prefm_global_d[k_ij+m*n_atom_pairs];
          hipDoubleComplex loc_emphi_rad_der;
          loc_emphi_rad_der.x=ilexp_der*tmp_prefm.x;
          loc_emphi_rad_der.y=ilexp_der*tmp_prefm.y;
          eimphi_rad_der_global_d[k_ij+k*n_atom_pairs]=loc_emphi_rad_der;

          hipDoubleComplex loc_emphi_azi_der;
          loc_emphi_azi_der.x=-loc_emphi.y;
          loc_emphi_azi_der.y= loc_emphi.x;
          eimphi_azi_der_global_d[k_ij+k*n_atom_pairs]=loc_emphi_azi_der;

          hipDoubleComplex loc_e_c_rad_der,loc_e_c_azi_der, loc_e_c_pol_der;
          loc_e_c_rad_der.x=amplitude*loc_preflm*plm_array_global_d[k_ij+k*n_atom_pairs]*loc_emphi_rad_der.x-
                            2.0* atom_sigma_scaling_d[i_sp-1]/atom_sigma*loc_exp_coeff.x;
          loc_e_c_rad_der.y=amplitude*loc_preflm*plm_array_global_d[k_ij+k*n_atom_pairs]*loc_emphi_rad_der.y-
                            2.0* atom_sigma_scaling_d[i_sp-1]/atom_sigma*loc_exp_coeff.y;
          
          loc_e_c_azi_der.x=amplitude*loc_preflm*plm_array_div_sin[k_ij+k*n_atom_pairs]*loc_emphi_azi_der.x;
          loc_e_c_azi_der.y=amplitude*loc_preflm*plm_array_div_sin[k_ij+k*n_atom_pairs]*loc_emphi_azi_der.y;     
          
          loc_e_c_pol_der.x=amplitude*loc_preflm*plm_array_der_mul_sin[k_ij+k*n_atom_pairs]*loc_emphi.x;
          loc_e_c_pol_der.y=amplitude*loc_preflm*plm_array_der_mul_sin[k_ij+k*n_atom_pairs]*loc_emphi.y;              

          // exp_coeff_rad_der_d[k+k_ij*kmax]=loc_e_c_rad_der;
          // exp_coeff_azi_der_d[k+k_ij*kmax]=loc_e_c_azi_der;
          // exp_coeff_pol_der_d[k+k_ij*kmax]=loc_e_c_pol_der;

          exp_coeff_rad_der_d[k_ij+k*n_atom_pairs]=loc_e_c_rad_der;
          exp_coeff_azi_der_d[k_ij+k*n_atom_pairs]=loc_e_c_azi_der;
          exp_coeff_pol_der_d[k_ij+k*n_atom_pairs]=loc_e_c_pol_der;

          k++;
        }
      }
    }
  }
}


__global__ void cuda_get_plm_arrays_der_one(double *plm_array_global_der_d,int kmax, int lmax, double *thetas_d, int n_atom_pairs,
                                            double *plm_array_div_sin, double *plm_array_der_mul_sin )
{
  int k_ij=threadIdx.x+blockIdx.x*blockDim.x;
  if(k_ij<n_atom_pairs){
    //double x=cos(thetas_d[k_ij]);
    double part1, part2;
    for(int l=0;l<=lmax;l++){
      for(int m=0; m<=l; m++){
        int k=1+l*(l+1)/2+m;
        int k_l_mp1=k+1;
        int k_l_mm1=k-1;
        int k_temp=-5;
        //       If m = 0 then we are asking for P_l^{-1}, which is not defined. We need
        //      to rewrite in terms of P_l^1:
        if(m==0){
          // P_0^1=0
          if(l==0){
            part1=0.0;
            // P_l^{-1} = - (l-1)!/(l+1)! * P_l^1
          }
          else{
            k_temp=1+l*(l+1)/2+1;
            part1= -0.5*plm_array_global_der_d[k_ij+(k_temp-1)*n_atom_pairs];
          }
        }
        else{
          part1=0.5*(l+m)*(l-m+1)*plm_array_global_der_d[k_ij+(k_l_mm1-1)*n_atom_pairs];
        }
        if(m==l){
          part2=0.0;
        }
        else{
          part2= -0.5*plm_array_global_der_d[k_ij+(k_l_mp1-1)*n_atom_pairs];
        }
        plm_array_der_mul_sin[k_ij+(k-1)*n_atom_pairs]=part1+part2;
      }
    }
    for(int l=0; l<=lmax;l++){
      for(int m=0; m<=l; m++){
        int k=1+l*(l+1)/2+m;
        if(m==0){
          plm_array_div_sin[k_ij+(k-1)*n_atom_pairs]=0.0;
        }
        else{
          int k_lp1_mp1 = 1 + (l+1)*(l+2)/2 + m + 1;
          int k_lp1_mm1 = 1 + (l+1)*(l+2)/2 + m - 1;
          part1=0.5*(l-m+1)*(l-m+2)*plm_array_global_der_d[k_ij+(k_lp1_mm1-1)*n_atom_pairs];
          part2=0.5*plm_array_global_der_d[k_ij+(k_lp1_mp1-1)*n_atom_pairs];
          plm_array_div_sin[k_ij+(k-1)*n_atom_pairs]=part1+part2;
        }
      }
    }
  }
}

extern "C" void  gpu_get_exp_coeff_array(hipDoubleComplex *eimphi_global_d, double *rjs_d,  double *phis_d,  double *thetas_d, 
                                             bool *mask_d, double *atom_sigma_in_d, double *atom_sigma_scaling_d, 
                                             double rcut, int n_atom_pairs, int n_species, int lmax, int kmax, 
                                             double *prefl_global_d, double *plm_array_global_d,double *plm_array_global_der_d, 
                                             double *prefl_global_der_d,
                                             double *preflm_d, hipDoubleComplex *exp_coeff_d, 
                                             bool c_do_derivatives, 
                                             hipDoubleComplex *eimphi_rad_der_global_d, hipDoubleComplex *eimphi_azi_der_global_d,
                                             double *plm_array_div_sin, double *plm_array_der_mul_sin, 
                                             hipDoubleComplex *exp_coeff_rad_der_d, hipDoubleComplex *exp_coeff_azi_der_d, hipDoubleComplex *exp_coeff_pol_der_d, hipStream_t *stream )
{
  dim3 nblocks=dim3((n_atom_pairs-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  hipDoubleComplex *prefm_global_d;
   
  gpuErrchk(hipMallocAsync(&prefm_global_d,  n_atom_pairs*(lmax+1)*sizeof(hipDoubleComplex) ,stream[0]));
  
  cuda_get_exp_coeff_one<<<nblocks, nthreads,0,stream[0]>>>(eimphi_global_d,rjs_d, phis_d,
                                                   mask_d, atom_sigma_in_d, atom_sigma_scaling_d,
                                                   rcut, n_atom_pairs, n_species, lmax, kmax,
                                                   prefl_global_d, prefm_global_d, preflm_d,
                                                   plm_array_global_d, exp_coeff_d); //, fact_array_d);


  
  if(c_do_derivatives){
   
  cuda_get_plm_arrays_der_one<<< nblocks, nthreads,0,stream[0]>>>(plm_array_global_der_d,kmax,lmax, thetas_d, n_atom_pairs,
                                                     plm_array_div_sin, plm_array_der_mul_sin );

  

  cuda_get_exp_coeff_der_one<<< nblocks, nthreads ,0,stream[0]>>>(eimphi_global_d, rjs_d, phis_d,
                                      mask_d, atom_sigma_in_d, atom_sigma_scaling_d,
                                      rcut, n_atom_pairs, n_species, lmax, kmax, 
                                      prefl_global_d, prefm_global_d, preflm_d, 
                                      prefl_global_der_d,
                                      plm_array_global_d, plm_array_global_der_d, exp_coeff_d,
                                      eimphi_rad_der_global_d, eimphi_azi_der_global_d,
                                      plm_array_div_sin, plm_array_der_mul_sin,
                                      exp_coeff_rad_der_d, exp_coeff_azi_der_d, exp_coeff_pol_der_d);
  
  
  }
  
  /*size_t free, total;
  hipMemGetInfo(& free, & total);
  counter++; 
  printf("\nFree memory %zu, from %zu in iteration %d\n", free/1024/1024, total/1024/1024, counter); */
  gpuErrchk(hipFreeAsync(prefm_global_d,stream[0]));
}


__global__
void check_nan(double *G, int Nt ){

 int id=threadIdx.x+threadIdx.x+blockIdx.x*blockDim.x;
 {
  if(isnan(G[id])&& id<Nt){
    printf("Is nan %lf at %d", G[id],id);
  }
 }

}



__global__
void cuda_global_scaling(double *radial_exp_coeff_d, 
                    int *i_beg_d, int *i_end_d, double *global_scaling_d,
                    int n_max, int n_atom_pairs, int n_species,
                    double *rcut_hard_d, int *k2_i_site_d, int *k2_start_d, int divide ){
                      
  int i_ij=threadIdx.x+blockIdx.x*blockDim.x;
  if(i_ij<n_atom_pairs){
    //int i_site=k2_i_site_d[i_ij];
    //int k2start=k2_start_d[i_site-1];
    //if(i_ij!=k2start)
    {
      int i_one=0;
      for(int i=0;i<n_species; i++){
        for(int ii=i_beg_d[i];ii<=i_end_d[i]; ii++){
          double loc_rad_exp_coeff=radial_exp_coeff_d[i_one+i_ij*n_max]*global_scaling_d[i]; //radial_exp_coeff_d[i_ij+i_one*size_radial_exp_coeff_two]*global_scaling_d[i];
          if(divide==0){
            loc_rad_exp_coeff*=sqrt(rcut_hard_d[i]);
          }
          if(divide==1){
            loc_rad_exp_coeff*=1.0/sqrt(rcut_hard_d[i]);
          } 
          radial_exp_coeff_d[i_one+i_ij*n_max]=loc_rad_exp_coeff; //radial_exp_coeff_d[i_ij+i_one*size_radial_exp_coeff_two]=loc_rad_exp_coeff;
          
          i_one++;
        }
      }  
    }
  }
}


__global__
void cuda_poly3gauss_one(double *radial_exp_coeff_d,
                    int *i_beg_d, int *i_end_d, double *global_scaling_d,
                    int n_max, int n_atom_pairs, int n_species,
                    double *rcut_hard_d, int *k2_i_site_d, int *k2_start_d){
                      
  int i_ij=threadIdx.x+blockIdx.x*blockDim.x;
  if(i_ij<n_atom_pairs){
    int i_site=k2_i_site_d[i_ij];
    int k2start=k2_start_d[i_site-1];
    if(i_ij!=k2start)
    {
      int i_one=0;
      for(int i=0;i<n_species; i++){
        for(int ii=i_beg_d[i];ii<=i_end_d[i]; ii++){
          double loc_rad_exp_coeff=0.0; //loc_rad_exp_coeff=radial_exp_coeff_d[i_one+i_ij*n_max]*global_scaling_d[i]; //radial_exp_coeff_d[i_ij+i_one*size_radial_exp_coeff_two]*global_scaling_d[i];
          
          loc_rad_exp_coeff*=sqrt(rcut_hard_d[i]);
          
          //radial_exp_coeff_d[i_one+i_ij*n_max]=loc_rad_exp_coeff; //radial_exp_coeff_d[i_ij+i_one*size_radial_exp_coeff_two]=loc_rad_exp_coeff;          
          i_one++;
        }
      }  
    }
  }
}

extern "C" void  gpu_get_radial_exp_coeff_poly3gauss(double *radial_exp_coeff_d, double *radial_exp_coeff_der_d, 
                                          int *i_beg_d, int *i_end_d, double *global_scaling_d,
                                          int size_radial_exp_coeff_one, int size_radial_exp_coeff_two, int n_species, 
                                          bool c_do_derivatives, int bintybint,
                                          double *rcut_hard_d,
                                          int *k2_i_site_d, int *k_2start_d,
                                          hipStream_t *stream ){
 
  dim3 nblocks=dim3((size_radial_exp_coeff_two-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1); 
  cuda_poly3gauss_one<<<nblocks, nthreads,0,stream[0]>>>(radial_exp_coeff_d, i_beg_d,i_end_d,global_scaling_d, 
                                          size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species,
                                          rcut_hard_d, k2_i_site_d, k_2start_d); 
  int divide;
  divide=0;
  cuda_global_scaling<<<nblocks, nthreads,0,stream[0]>>>(radial_exp_coeff_d, i_beg_d,i_end_d,global_scaling_d, 
                                          size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species,
                                          rcut_hard_d, k2_i_site_d, k_2start_d, divide);  
  /* gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() ); */
  if(c_do_derivatives){
    divide=1;
    cuda_global_scaling<<< nblocks, nthreads,0,stream[0] >>>(radial_exp_coeff_der_d, i_beg_d,i_end_d,global_scaling_d, 
                                size_radial_exp_coeff_one, size_radial_exp_coeff_two, n_species,
                                rcut_hard_d, k2_i_site_d, k_2start_d, divide);    
  } 
  /* gpuErrchk( hipPeekAtLastError() );
  gpuErrchk( hipDeviceSynchronize() );   */                          
}

__global__
void kernel_get_radial_poly3gauss(int n_atom_pairs, int n_species, bool *mask_d, double *rjs_d, double *rcut_hard_d, 
		                  int n_sites, int *n_neigh_d, int n_max, int n_temp, bool do_derivatives, double *exp_coeff_d, 
				  double *exp_coeff_der_d, double *rcut_soft_d, double *atom_sigma_d, double *exp_coeff_temp1_d,
				  double *exp_coeff_temp2_d,double *exp_coeff_der_temp_d, int *i_beg_d, int *i_end_d,
				  double *atom_sigma_scaling_d, int mode, int radial_enhancement, double *amplitude_scaling_d, 
				  int *alpha_max_d, double *nf_d, int n_temp_der, double *W_d){

  int n,d,k;
  double  ampli_tude, ampli_tude_der, atom_sigma_scaled, amplitude_scaling, C1, C2, W_exp, nf, atom_sigma_f, rj_f, sf2;
  double tmp1, tmp2, tmp3, tmp4, tmp5,tmp6,tmp7;

  int k_ij=threadIdx.x+blockIdx.x*blockDim.x;
  if(k_ij<n_atom_pairs){
    double rjs=rjs_d[k_ij];
    for(int i_sp=0;i_sp<n_species; i_sp++){
      if(mask_d[k_ij+i_sp*n_atom_pairs]){
        double rcut_hard_in = rcut_hard_d[i_sp];
        if(rjs<=rcut_hard_in){
	  int alpha_max = alpha_max_d[i_sp];
          int alpha_max_der = alpha_max;
          if (do_derivatives) alpha_max_der += 2;
          int j=0;
          int j1=0;
          for (n=0; n<n_sites; n++){
            j1 += n_neigh_d[n];
            if (k_ij<j1) break;
            j=j1;
          }
          if (k_ij>j) {
            double pi = acos(-1.0);
            double sq2 = sqrt(2.0);
	    double rcut_soft_in = rcut_soft_d[i_sp];
            double rcut_soft = rcut_soft_in/rcut_hard_in;
            double rcut_hard = 1.0;
	    double atom_sigma_in=atom_sigma_d[i_sp];
            double atom_sigma = atom_sigma_in/rcut_hard_in;
            double dr = 1.0 - rcut_soft_in/rcut_hard_in;
            double N_gauss = sqrt(2.0/atom_sigma) / pow(pi,0.25);
            double pref_f = 0.0;
            for (n=0; n<alpha_max_der; n++) {
              exp_coeff_temp1_d[k_ij*n_temp+n] = 0.0;
              exp_coeff_temp2_d[k_ij*n_temp+n] = 0.0;
            }
            if (do_derivatives)
               for (n=0; n<alpha_max; n++) exp_coeff_der_temp_d[k_ij*n_temp_der+n] = 0.0;
            double rj = rjs/rcut_hard_in;
	    double atom_sigma_scaling=atom_sigma_scaling_d[i_sp];
            double atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj;
            double s2 = pow(atom_sigma_scaled,2);
	    double amplitude_scaling = amplitude_scaling_d[i_sp];
	    tmp1 = 1.0 + rj * rj * (2.0 * rj - 3.0);
            tmp2 = atom_sigma_scaling / atom_sigma_scaled;
	    tmp3 = 6.0 / atom_sigma_scaled * rj * (rj - 1.0);
	    if (mode==mode_polynomial) {
              if( amplitude_scaling == 0.0 ){
                ampli_tude = 1.0 / atom_sigma_scaled;
                ampli_tude_der = - atom_sigma_scaling / s2;
              } else if( tmp1 <= 1.e-10 ){
                ampli_tude = 0.0;
                ampli_tude_der = 0.0;
              } else {
                if( amplitude_scaling == 1.0 ){
                  ampli_tude = 1.0 / atom_sigma_scaled * tmp1;
                  ampli_tude_der = tmp3 - tmp2 * ampli_tude;
                } else {
                  ampli_tude = 1.0 / atom_sigma_scaled * pow(tmp1,amplitude_scaling);
                  ampli_tude_der = amplitude_scaling * tmp3 * pow(tmp1,amplitude_scaling - 1.0) - tmp2 -ampli_tude;
                }
              }
            }
	    tmp3 = rj + sqrt(2.0/pi)*atom_sigma_scaled;
            tmp4 = sqrt(8.0/pi)*atom_sigma_scaled;
            tmp5 = rj*rj + s2 + tmp4*rj;
            if( radial_enhancement == 1 ){
              ampli_tude_der = ampli_tude * ( 1.0 + sqrt(2.0/pi)*atom_sigma_scaling ) + ampli_tude_der * tmp3;
              ampli_tude = ampli_tude * tmp3;
            } else if( radial_enhancement == 2 ){
              ampli_tude_der = ampli_tude*( 2.0*rj + 2.0*atom_sigma_scaled*atom_sigma_scaling + tmp4 + sqrt(8.0/pi)*rj*atom_sigma_scaling )
                               + ampli_tude_der*tmp5;
              ampli_tude = ampli_tude * tmp5;
            }
	    double I_n = 0.0;
            double N_n = 1.0;
            double N_np1 = N_a(rcut_hard, -2);
            double I_np1 = sqrt(pi/2.0) * atom_sigma_scaled * ( erf( (rcut_soft-rj)/sq2/atom_sigma_scaled ) - erf( (-rj)/sq2/atom_sigma_scaled ) ) / N_np1;
	    double I_np2, N_np2;
	    C1 = (rcut_hard_in == rcut_soft_in) ? 0.0 : s2 / dr * exp(-0.5 * pow(rcut_soft - rj,2) / s2);
            C2 = s2 / rcut_hard * exp(-0.5 * pow(rj,2) / s2);
            for (n = -1; n<=alpha_max_der-1;n++){
              C1 = C1 * dr;
              C2 = C2 * rcut_hard;
              N_np2 = N_a(rcut_hard, n);
              I_np2 = s2 * double(n+1) * N_n/ N_np2 * I_n - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 + C1 / N_np2  - C2 / N_np2;
              if(n > 0) exp_coeff_temp1_d[k_ij*n_temp+n-1] = I_np2;
              N_n = N_np1;
              N_np1 = N_np2;
              I_n = I_np1;
              I_np1 = I_np2;
	    }
	    if( do_derivatives ) {
              tmp1 = atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled;
              tmp2 = (rj - rcut_hard) / s2 * (tmp1 - 1.0);
              tmp3 = rcut_hard * ( 2.0 * tmp1 - 1.0 ) / s2;
              tmp4 = atom_sigma_scaling * rcut_hard * rcut_hard / pow(atom_sigma_scaled,3);
	      tmp5 = exp_coeff_temp1_d[k_ij*n_temp];
	      tmp6 = exp_coeff_temp1_d[k_ij*n_temp+1];
              for (n = 1; n<=alpha_max-1; n++){
	        tmp7 = exp_coeff_temp1_d[k_ij*n_temp+n+1];
                exp_coeff_der_temp_d[k_ij*n_temp_der+n-1] = tmp2 * tmp5 + tmp3 * N_a(rcut_hard, n+1) / N_a(rcut_hard, n) * tmp6 +
                                                            tmp4 * N_a(rcut_hard, n+2) / N_a(rcut_hard, n) * tmp7;
		tmp5 = tmp6;
		tmp6 = tmp7;
	      }
            }
            if (false || (rcut_soft - rj) < 4.0*atom_sigma_scaled) {
	      nf = nf_d[i_sp];
              tmp1 = dr * dr / nf / nf;
              atom_sigma_f = atom_sigma_scaled * dr / nf / sqrt(s2 + tmp1);
              rj_f = (s2 * rcut_soft + tmp1 * rj) / (s2 + tmp1);
              sf2 = pow(atom_sigma_f,2);
              pref_f = exp( -0.5 * pow(rcut_soft-rj,2) / ( s2 + tmp1) );
              I_n = 0.0;
              N_n = 1.0;
              N_np1 = N_a(rcut_hard, -2);
              I_np1 = sqrt(pi/2.0) * atom_sigma_f * ( erf( (rcut_hard-rj_f)/sq2/atom_sigma_f ) - erf( (rcut_soft-rj_f)/sq2/atom_sigma_f ) ) / N_np1;
              C2 = sf2 / dr * exp(-0.5 * pow(rcut_soft - rj_f,2) / sf2);
              for (n = -1; n<=alpha_max_der-1; n++){
                C2 *= dr;
                double N_np2 = N_a(rcut_hard, n);
                double I_np2 = sf2 * double(n+1) * N_n/ N_np2 * I_n - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1  - C2 / N_np2;
                if(n > 0) exp_coeff_temp2_d[k_ij*n_temp+n-1] += I_np2;
                N_n = N_np1;
                N_np1 = N_np2;
                I_n = I_np1;
                I_np1 = I_np2;
              }
	      if (do_derivatives) {
                double denom = s2 + tmp1;
                double der_pref_f = pref_f * ( (rcut_soft - rj) / denom + pow(rcut_soft - rj,2) / pow(denom,2)* atom_sigma_scaled * atom_sigma_scaling );
                double der_rjf_rj = (2.0*atom_sigma_scaled*rcut_soft*atom_sigma_scaling + tmp1) / denom - (s2*rcut_soft + tmp1 * rj) * 2.0 * 
			            atom_sigma_scaled * atom_sigma_scaling / pow(denom,2);
                double der_sjf_rj = atom_sigma_scaling * dr/nf / sqrt(denom) * (1.0 - pow(atom_sigma_scaled,2)/denom);
                tmp2 = (rj_f - rcut_hard) / sf2 * ( der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj );
                tmp3 = rcut_hard / sf2 * ( 2.0 * der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj );
                tmp4 = der_sjf_rj * rcut_hard * rcut_hard / pow(atom_sigma_f,3);
		tmp5 = exp_coeff_temp2_d[k_ij*n_temp];
		tmp6 = exp_coeff_temp2_d[k_ij*n_temp+1];
                for (n = 1; n <=alpha_max-1;n++){
	          tmp7 = exp_coeff_temp2_d[k_ij*n_temp+n+1];
                  exp_coeff_der_temp_d[k_ij*n_temp_der+n-1] += pref_f * ( tmp2 * tmp5 + tmp3 * N_a(rcut_hard, n+1) / N_a(rcut_hard, n) * tmp6 +
                                                              tmp4 * N_a(rcut_hard, n+2) / N_a(rcut_hard, n) * tmp7) + der_pref_f * tmp5;
		  tmp5 = tmp6;
		  tmp6 = tmp7;
		}
              }
            }
	    exp_coeff_temp1_d[k_ij*n_temp+alpha_max-1] = 0.0;
            exp_coeff_temp2_d[k_ij*n_temp+alpha_max-1] = 0.0;
            if (false || rj < 4.0*(atom_sigma+atom_sigma_scaled)) {
              double sigma_star = sqrt(pow(atom_sigma,2) + s2);
	      exp_coeff_temp1_d[k_ij*n_temp+alpha_max-1]= exp(- 0.5 * pow(rj,2) / pow(sigma_star,2) ) * sqrt(pi/2.0) * 
                                                          atom_sigma_scaled*atom_sigma / sigma_star * ( 1.0 + 
						          erf(atom_sigma/atom_sigma_scaled*rj/sq2/sigma_star) )* N_gauss;
              if (do_derivatives)
                exp_coeff_der_temp_d[k_ij*n_temp_der+alpha_max-1] = ( pow(rj,2) * atom_sigma_scaling / pow(atom_sigma_scaled,3) - 
				                                    rj/pow(sigma_star,2) + atom_sigma_scaling*pow(rj,2)*pow(atom_sigma,4)/
								    pow(atom_sigma_scaled,3)/pow(sigma_star,4) + atom_sigma_scaling*
								    pow(atom_sigma,2)/atom_sigma_scaled/pow(sigma_star,2) - 2.0*pow(rj,2)*
								    atom_sigma_scaling*pow(atom_sigma,2)/pow(atom_sigma_scaled,3)/
								    pow(sigma_star,2) ) * exp_coeff_temp1_d[k_ij*n_temp+alpha_max-1]+
                                                                    (1./s2 - 2.0*rj*atom_sigma_scaling/pow(atom_sigma_scaled,3)) * s2 * 
								    pow(atom_sigma,2) / pow(sigma_star,2) * sqrt(2.0/atom_sigma) / 
								    pow(pi,0.25) * exp(-0.5 * pow(rj,2) / pow(sigma_star,2) * (1.0 + 
								    pow(atom_sigma,2) / s2) ) + sqrt(2.0/atom_sigma) / pow(pi,0.25) * exp(-0.5*
								    pow( rj,2) / pow(sigma_star,2) * (1.0 + pow(atom_sigma,2) / s2) ) * 
								    atom_sigma_scaling / atom_sigma_scaled * rj*pow(atom_sigma,4)/pow(sigma_star,4);
            }
	    if (do_derivatives) {
              for (n=0; n<alpha_max; n++)
                exp_coeff_der_temp_d[k_ij*n_temp_der+n] = ampli_tude * exp_coeff_der_temp_d[k_ij*n_temp_der+n] +
                                                          ampli_tude_der * (exp_coeff_temp1_d[k_ij*n_temp+n] + pref_f * exp_coeff_temp2_d[k_ij*n_temp+n]);
	      for (d=i_beg_d[i_sp]; d<=i_end_d[i_sp]; d++){
	        W_exp = 0.0;
		k = 0;
		for (n=i_beg_d[i_sp]; n<=i_end_d[i_sp];n++){
	          W_exp += W_d[(n-1)*n_max+d-1]* exp_coeff_der_temp_d[k_ij*n_temp_der+k];
		  k +=1;
		}
		exp_coeff_der_d[k_ij*n_max+d-1]=W_exp;
              }
            }
  	    for (d=i_beg_d[i_sp]; d<=i_end_d[i_sp]; d++){
	      W_exp = 0.0;
	      k=0;
	      for (n=i_beg_d[i_sp]; n<=i_end_d[i_sp];n++){
	        W_exp += W_d[(n-1)*n_max+d-1]*(exp_coeff_temp1_d[k_ij*n_temp+k]+ pref_f*exp_coeff_temp2_d[k_ij*n_temp+k]);
		k+=1;
	      }
	      exp_coeff_d[k_ij*n_max+d-1]=ampli_tude*W_exp;
            }
	  }
	}
      }
    }
  }
}

extern "C" void  gpu_radial_poly3gauss(int n_atom_pairs, int n_species, bool *mask_d, double *rjs_d, double *rcut_hard_d, 
		                       int n_sites, int *n_neigh_d, int n_max, int n_temp, bool do_derivatives, double *exp_coeff_d,
                                       double *exp_coeff_der_d, double *rcut_soft_d, double *atom_sigma_d, double *exp_coeff_temp1_d,
                                       double *exp_coeff_temp2_d,double *exp_coeff_der_temp_d, int *i_beg, int *i_end, 
				       double *atom_sigma_scaling_d, int mode, int radial_enhancement, double *amplitude_scaling_d, 
				       int *alpha_max_d, double *nf_d, int n_temp_der, double *W_d,hipStream_t *stream ){

  dim3 nblocks=dim3((n_atom_pairs-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  kernel_get_radial_poly3gauss<<<nblocks,nthreads,0,stream[0] >>>(n_atom_pairs,n_species,mask_d,rjs_d,rcut_hard_d,n_sites,
	                                                          n_neigh_d, n_max, n_temp, do_derivatives, exp_coeff_d,
                                                                  exp_coeff_der_d, rcut_soft_d, atom_sigma_d, exp_coeff_temp1_d,
                                                                  exp_coeff_temp2_d, exp_coeff_der_temp_d, i_beg, i_end,
							          atom_sigma_scaling_d, mode, radial_enhancement, amplitude_scaling_d,
							          alpha_max_d, nf_d, n_temp_der, W_d);
}

__global__
void kernel_get_radial_poly3(int n_atom_pairs, int n_species, bool *mask_d, double *rjs_d, double *rcut_hard_d,
                                  int n_sites, int *n_neigh_d, int n_max, int n_temp, bool do_derivatives, double *exp_coeff_d,
                                  double *exp_coeff_der_d, double *rcut_soft_d, double *atom_sigma_d, double *exp_coeff_temp1_d,
                                  double *exp_coeff_temp2_d,double *exp_coeff_der_temp_d, int *i_beg_d, int *i_end_d,
                                  double *atom_sigma_scaling_d, int mode, int radial_enhancement, double *amplitude_scaling_d,
                                  int *alpha_max_d, double *nf_d, int n_temp_der, double *W_d, bool *do_central_d, double *central_weight_d){

  int n,d,k;
  double  ampli_tude, ampli_tude_der, atom_sigma_scaled, amplitude_scaling, C1, C2, W_exp, nf, atom_sigma_f, rj_f, sf2;
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  int k_ij=threadIdx.x+blockIdx.x*blockDim.x;
  if(k_ij<n_atom_pairs){
    double rjs=rjs_d[k_ij];
    for(int i_sp=0;i_sp<n_species; i_sp++){
      if(mask_d[k_ij+i_sp*n_atom_pairs]){
        double rcut_hard_in = rcut_hard_d[i_sp];
        if(rjs<=rcut_hard_in){
          int alpha_max = alpha_max_d[i_sp];
          int alpha_max_der = alpha_max;
          if (do_derivatives) alpha_max_der += 2;
          int j=0;
          int j1=0;
          for (n=0; n<n_sites; n++){
            j1 += n_neigh_d[n];
            if (k_ij<j1) break;
            j=j1;
          }
          if (k_ij>j || do_central_d[i_sp] ) {
            double pi = acos(-1.0);
            double sq2 = sqrt(2.0);
            double rcut_soft_in = rcut_soft_d[i_sp];
            double rcut_soft = rcut_soft_in/rcut_hard_in;
            double rcut_hard = 1.0;
            double atom_sigma_in=atom_sigma_d[i_sp];
            double atom_sigma = atom_sigma_in/rcut_hard_in;
            double dr = 1.0 - rcut_soft_in/rcut_hard_in;
            double pref_f = 0.0;
            for (n=0; n<alpha_max_der; n++) {
              exp_coeff_temp1_d[k_ij*n_temp+n] = 0.0;
              exp_coeff_temp2_d[k_ij*n_temp+n] = 0.0;
	    }
            if (do_derivatives)
               for (n=0; n<alpha_max; n++) exp_coeff_der_temp_d[k_ij*n_temp_der+n] = 0.0;
            double rj = rjs/rcut_hard_in;
            double atom_sigma_scaling=atom_sigma_scaling_d[i_sp];
            double atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj;
            double s2 = pow(atom_sigma_scaled,2);
            double amplitude_scaling = amplitude_scaling_d[i_sp];
	    tmp1 = 1.0 + rj * rj * (2.0 * rj - 3.0);
            tmp2 = atom_sigma_scaling / atom_sigma_scaled;
	    tmp3 = 6.0 / atom_sigma_scaled * rj * (rj - 1.0);
            if (mode==mode_polynomial) {
              if( amplitude_scaling == 0.0 ){
                ampli_tude = 1.0 / atom_sigma_scaled;
                ampli_tude_der = - atom_sigma_scaling / s2;
              } else if( tmp1 <= 1.e-10 ){
                ampli_tude = 0.0;
                ampli_tude_der = 0.0;
              } else {
                if( amplitude_scaling == 1.0 ){
                  ampli_tude = 1.0 / atom_sigma_scaled * tmp1;
                  ampli_tude_der = tmp3 - tmp2 * ampli_tude;
                } else {
                  ampli_tude = 1.0 / atom_sigma_scaled * pow(tmp1,amplitude_scaling);
                  ampli_tude_der = tmp3 * amplitude_scaling * pow(tmp1,amplitude_scaling - 1.0) - tmp2 * ampli_tude;
                }
              }
            }
	    if (k_ij==j) {
              ampli_tude = central_weight_d[i_sp] * ampli_tude;
              ampli_tude_der = central_weight_d[i_sp] * ampli_tude_der;
	    }
	    tmp3 = rj + sqrt(2.0/pi)*atom_sigma_scaled;
            tmp4 = sqrt(8.0/pi)*atom_sigma_scaled;
            tmp5 = rj*rj + s2 + tmp4*rj;
            if( radial_enhancement == 1 ){
              ampli_tude_der = ampli_tude * ( 1.0 + sqrt(2.0/pi)*atom_sigma_scaling ) + ampli_tude_der * tmp3;
              ampli_tude = ampli_tude * tmp3;
            } else if( radial_enhancement == 2 ){
              ampli_tude_der = ampli_tude*( 2.0*rj + 2.0*atom_sigma_scaled*atom_sigma_scaling + tmp4 + sqrt(8.0/pi)*rj*atom_sigma_scaling )
                               + ampli_tude_der*tmp5;
              ampli_tude = ampli_tude * tmp5;
            }
            double I_n = 0.0;
            double N_n = 1.0;
            double N_np1 = N_a(rcut_hard, -2);
            double I_np1 = sqrt(pi/2.0) * atom_sigma_scaled * ( erf( (rcut_soft-rj)/sq2/atom_sigma_scaled ) - erf( (-rj)/sq2/atom_sigma_scaled ) ) / N_np1;
            double I_np2, N_np2;
	    C1 = (rcut_hard_in == rcut_soft_in) ? 0.0 : s2 / dr * exp(-0.5 * pow(rcut_soft - rj,2) / s2);
            C2 = s2 / rcut_hard * exp(-0.5 * pow(rj,2) / s2);
            for (n = -1; n<=alpha_max_der;n++){
              C1 *= dr;
              C2 *= rcut_hard;
              N_np2 = N_a(rcut_hard, n);
              I_np2 = s2 * double(n+1) * N_n/ N_np2 * I_n - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 + C1 / N_np2  - C2 / N_np2;
              if(n > 0) exp_coeff_temp1_d[k_ij*n_temp+n-1] = I_np2;
              N_n = N_np1;
              N_np1 = N_np2;
              I_n = I_np1;
              I_np1 = I_np2;
            }
	    if( do_derivatives ) {
	      tmp1 = atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled;
              tmp2 = (rj - rcut_hard) / s2 * (tmp1 - 1.0);
              tmp3 = rcut_hard * ( 2.0 * tmp1 - 1.0 ) / s2;
              tmp4 = atom_sigma_scaling * rcut_hard * rcut_hard / pow(atom_sigma_scaled,3);
              tmp5 = exp_coeff_temp1_d[k_ij*n_temp];
              tmp6 = exp_coeff_temp1_d[k_ij*n_temp+1];
              for (n = 1; n<=alpha_max-1; n++) {
		tmp7 = exp_coeff_temp1_d[k_ij*n_temp+n+1];
                exp_coeff_der_temp_d[k_ij*n_temp_der+n-1] = tmp2 * tmp5 + tmp3 * N_a(rcut_hard, n+1) / N_a(rcut_hard, n) * tmp6 +
                                                            tmp4 * N_a(rcut_hard, n+2) / N_a(rcut_hard, n) * tmp7;
                tmp5 = tmp6;
                tmp6 = tmp7;
              }
            }
            if (false || (rcut_soft - rj) < 4.0*atom_sigma_scaled) {
	      nf = nf_d[i_sp];
              tmp1 = dr * dr / nf / nf;
              atom_sigma_f = atom_sigma_scaled * dr / nf / sqrt(s2 + tmp1);
              rj_f = (s2 * rcut_soft + tmp1 * rj) / (s2 + tmp1);
              sf2 = pow(atom_sigma_f,2);
              pref_f = exp( -0.5 * pow(rcut_soft-rj,2) / ( s2 + tmp1) );
              I_n = 0.0;
              N_n = 1.0;
              N_np1 = N_a(rcut_hard, -2);
              I_np1 = sqrt(pi/2.0) * atom_sigma_f * ( erf( (rcut_hard-rj_f)/sq2/atom_sigma_f ) - erf( (rcut_soft-rj_f)/sq2/atom_sigma_f ) ) / N_np1;
              C2 = sf2 / dr * exp(-0.5 * pow(rcut_soft - rj_f,2) / sf2);
              for (n = -1; n<=alpha_max_der-1; n++){
                C2 *= dr;
                double N_np2 = N_a(rcut_hard, n);
                double I_np2 = sf2 * double(n+1) * N_n/ N_np2 * I_n - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1  - C2 / N_np2;
                if(n > 0) exp_coeff_temp2_d[k_ij*n_temp+n-1] += I_np2;
                N_n = N_np1;
                N_np1 = N_np2;
                I_n = I_np1;
                I_np1 = I_np2;
              }
	      if (do_derivatives) {
                double denom = s2 + tmp1;
                double der_pref_f = pref_f * ( (rcut_soft - rj) / denom + pow(rcut_soft - rj,2) / pow(denom,2)* atom_sigma_scaled * atom_sigma_scaling );
                double der_rjf_rj = (2.0*atom_sigma_scaled*rcut_soft*atom_sigma_scaling + tmp1) / denom - (s2*rcut_soft + tmp1 * rj) * 2.0 *
                                    atom_sigma_scaled * atom_sigma_scaling / pow(denom,2);
                double der_sjf_rj = atom_sigma_scaling * dr/nf / sqrt(denom) * (1.0 - pow(atom_sigma_scaled,2)/denom);

		tmp2 = (rj_f - rcut_hard) / sf2 * ( der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj );
                tmp3 = rcut_hard / sf2 * ( 2.0 * der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj );
                tmp4 = der_sjf_rj * rcut_hard * rcut_hard / pow(atom_sigma_f,3);
                tmp5 = exp_coeff_temp2_d[k_ij*n_temp];
                tmp6 = exp_coeff_temp2_d[k_ij*n_temp+1];
                for (n = 1; n <=alpha_max-1;n++){
                  tmp7 = exp_coeff_temp2_d[k_ij*n_temp+n+1];
                  exp_coeff_der_temp_d[k_ij*n_temp_der+n-1] += pref_f * ( tmp2 * tmp5 + tmp3 * N_a(rcut_hard, n+1) / N_a(rcut_hard, n) * tmp6 +
                                                              tmp4 * N_a(rcut_hard, n+2) / N_a(rcut_hard, n) * tmp7) + der_pref_f * tmp5;
                  tmp5 = tmp6;
                  tmp6 = tmp7;
                }
              }
            }
            if (do_derivatives) {
              for (n=0; n<alpha_max; n++)
                exp_coeff_der_temp_d[k_ij*n_temp_der+n] = ampli_tude * exp_coeff_der_temp_d[k_ij*n_temp_der+n] +
                                                          ampli_tude_der * (exp_coeff_temp1_d[k_ij*n_temp+n] + pref_f * exp_coeff_temp2_d[k_ij*n_temp+n]);
              for (d=i_beg_d[i_sp]; d<=i_end_d[i_sp]; d++){
                W_exp = 0.0;
                k = 0;
                for (n=i_beg_d[i_sp]; n<=i_end_d[i_sp];n++){
                  W_exp += W_d[(n-1)*n_max+d-1]* exp_coeff_der_temp_d[k_ij*n_temp_der+k];
                  k +=1;
                }
                exp_coeff_der_d[k_ij*n_max+d-1]=W_exp;
              }
            }
            for (d=i_beg_d[i_sp]; d<=i_end_d[i_sp]; d++){
              W_exp = 0.0;
              k=0;
              for (n=i_beg_d[i_sp]; n<=i_end_d[i_sp];n++){
                W_exp += W_d[(n-1)*n_max+d-1]*(exp_coeff_temp1_d[k_ij*n_temp+k]+ pref_f*exp_coeff_temp2_d[k_ij*n_temp+k]);
                k+=1;
              }
              exp_coeff_d[k_ij*n_max+d-1]=ampli_tude*W_exp;
            }
          }
        }
      }
    }
  }
}

extern "C" void  gpu_radial_poly3(int n_atom_pairs, int n_species, bool *mask_d, double *rjs_d, double *rcut_hard_d,
                                  int n_sites, int *n_neigh_d, int n_max, int n_temp, bool do_derivatives, double *exp_coeff_d,
                                  double *exp_coeff_der_d, double *rcut_soft_d, double *atom_sigma_d, double *exp_coeff_temp1_d,
                                  double *exp_coeff_temp2_d,double *exp_coeff_der_temp_d, int *i_beg, int *i_end,
                                  double *atom_sigma_scaling_d, int mode, int radial_enhancement, double *amplitude_scaling_d,
                                  int *alpha_max_d, double *nf_d, int n_temp_der, double *W_d, bool *do_central_d, 
				  double *central_weight_d, hipStream_t *stream ){

  dim3 nblocks=dim3((n_atom_pairs-1+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  kernel_get_radial_poly3<<<nblocks,nthreads,0,stream[0] >>>(n_atom_pairs,n_species,mask_d,rjs_d,rcut_hard_d,n_sites,
                                                             n_neigh_d, n_max, n_temp, do_derivatives, exp_coeff_d,
                                                             exp_coeff_der_d, rcut_soft_d, atom_sigma_d, exp_coeff_temp1_d,
                                                             exp_coeff_temp2_d, exp_coeff_der_temp_d, i_beg, i_end,
                                                             atom_sigma_scaling_d, mode, radial_enhancement, amplitude_scaling_d,
                                                             alpha_max_d, nf_d, n_temp_der, W_d, do_central_d, central_weight_d);
}

__global__
void kernel_get_2b(int i_beg, int i_end, int n_sparse, double *energies_d, double e0, int *n_neigh_d, bool do_forces, double *forces_d,
	           double *virial_d, double *rjs_d, double rcut, int *species_d, int *neighbor_species_d, int sp1, int sp2,double buffer, 
		   double delta, double *cutoff_d, double *Qs_d, double sigma, double *alphas_d, double *xyz_d){

  int i_site=i_beg-1+threadIdx.x+blockIdx.x*blockDim.x;
  double forces_loc[3], virial_loc[9],energies_loc;
  int i,j,k,s,k1,k2;
  double rjs_k,fcut,dfcut,pi,sigma2,delta2,tmp;
  if(i_site<i_end){
    energies_loc = e0;
    if( do_forces ){
      for (i=0; i <3; i++) forces_loc[i] = 0.0;
      for (i=0; i <9; i++) virial_loc[i] = 0.0;
    }
    if( species_d[i_site] == sp1 ||  species_d[i_site] == sp2 ) {
      pi = acos(-1.0);
      k=0;
      for (i=i_beg-1; i<i_site; i++) k += n_neigh_d[i];
      for (j=2; j<=n_neigh_d[i_site]; j++) {
        k +=1;
        if( !((species_d[i_site]==sp1 && neighbor_species_d[k]==sp2) || (species_d[i_site]==sp2 && neighbor_species_d[k]==sp1)) ) continue;
        rjs_k = rjs_d[k];
        if( rjs_k < rcut ) {
	  fcut=( rjs_k<rcut-buffer) ? 1.0 : (cos(pi*(rjs_k-rcut+buffer)/buffer)+1.0)/2.0;
	  if( do_forces ) dfcut=( rjs_k<rcut-buffer) ? 0.0 : pi/2.0/buffer*sin(pi*(rjs_k-rcut+buffer)/buffer);
	  sigma2=sigma*sigma;
	  delta2=delta*delta;
          for (s = 0;s<n_sparse;s++){
            tmp = delta2*alphas_d[s]*cutoff_d[s]* exp(-0.5*pow(rjs_k-Qs_d[s],2)/sigma2);
            energies_loc += tmp*fcut;
            if( do_forces) {
              for (i=0;i<3;i++){
                forces_loc[i] = -2.0*tmp*xyz_d[3*k+i]/rjs_k*((rjs_k-Qs_d[s])/sigma2*fcut+dfcut);
                forces_d[3*i_site+i] += forces_loc[i];
	      }
              for (k2=0; k2<3; k2++)
                for (k1=0; k1<3; k1++)
                  virial_loc[3*k2+k1] += -0.5*(forces_loc[k1]*xyz_d[3*k+k2]+forces_loc[k2]*xyz_d[3*k+k1]);
            }
          }
        }
      }
    }
    energies_d[i_site] = energies_loc;
    for (k1=0;k1<9;k1++) atomicAdd(&virial_d[k1],0.5*virial_loc[k1]);
  }
}

extern "C" void  gpu_get_2b_forces_energies(int i_beg, int i_end, int n_sparse, double *energies_d, double e0, int *n_neigh_d, bool do_forces,
                                            double *forces_d, double *virial_d,double *rjs_d, double rcut, int *species_d,
                                            int *neighbor_species_d, int sp1, int sp2, double buffer, double delta, double *cutoff_d, 
					    double *Qs_d, double sigma, double *alphas_d, double *xyz_d, hipStream_t *stream ){

  dim3 nblocks=dim3((i_end-i_beg+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  kernel_get_2b<<<nblocks, nthreads,0,stream[0] >>>(i_beg, i_end, n_sparse, energies_d, e0, n_neigh_d, do_forces, forces_d, virial_d, 
		                                    rjs_d, rcut, species_d, neighbor_species_d, sp1, sp2, buffer, delta, cutoff_d, 
						    Qs_d, sigma, alphas_d, xyz_d);
  //temporary, to measure timings
  hipStreamSynchronize(stream[0]);
}

__global__
void kernel_get_core_pot(int i_beg, int i_end, bool do_forces, int *species_d, int sp1, int sp2, int *n_neigh_d, int *neighbor_species_d, 
		         double *rjs_d, int n_sparse, double *x_d, double *V_d, double *dVdx2_d, double yp1, double ypn, double *xyz_d, 
			 double *forces_d, double *virial_d, double *energies_d){

  int i_site=i_beg-1+threadIdx.x+blockIdx.x*blockDim.x;
  double forces_loc[3], virial_loc[9],energies_loc;
  int i,j,k,k1,k2;
  double rjs_k, rcut, Vint, d_Vint;

  if(i_site<i_end){
    energies_loc = 0.0;
    if( do_forces ){
      for (i=0; i <3; i++) forces_loc[i] = 0.0;
      for (i=0; i <9; i++) virial_loc[i] = 0.0;
    }
    if( species_d[i_site] == sp1 ||  species_d[i_site] == sp2 ) {
      k=0;
      for (i=i_beg-1; i<i_site; i++) k += n_neigh_d[i];
      for (j=1; j < n_neigh_d[i_site]; j++) {
        k +=1;
	if( !((species_d[i_site]==sp1 && neighbor_species_d[k]==sp2) || (species_d[i_site]==sp2 && neighbor_species_d[k]==sp1)) ) continue;
        rjs_k = rjs_d[k];
        rcut=x_d[0];
        for (i=1; i<n_sparse; i++)
	  if (rcut<=x_d[i]) rcut=x_d[i];

	//	printf("GPU_CORE_POT: rjs_k = %lf, k = %d, rcut = %lf, i_site = %d \n",rjs_k, k, rcut,  i_site);    	  	
        if( rjs_k < rcut ) {
	  energies_loc += 0.5 * gpu_spline(n_sparse, rjs_k, x_d, V_d, dVdx2_d, rcut, yp1, ypn);

          if( do_forces) {
	    d_Vint=gpu_spline_der(n_sparse, rjs_k, x_d, V_d, dVdx2_d, rcut, yp1, ypn);
            for (i=0; i<3; i++) {
              forces_loc[i] = d_Vint * xyz_d[3*k+i] / rjs_k;
              forces_d[3*i_site+i] += forces_loc[i];
            }
            for (k2=0; k2<3; k2++)
              for (k1=0; k1<3; k1++)
                virial_loc[3*k2+k1] += - 0.5 * (forces_loc[k1]*xyz_d[3*k+k2] + forces_loc[k2]*xyz_d[3*k+k1]);
          }
        }
      }
    }
    energies_d[i_site] = energies_loc;
    for (k1=0;k1<9;k1++) atomicAdd(&virial_d[k1],0.5*virial_loc[k1]);
  }
}

extern "C" void  gpu_get_core_pot_energy_and_forces(int i_beg, int i_end, bool do_forces, int *species_d, int sp1, int sp2, int *n_neigh_d, 
		                                    int *neighbor_species_d, double *rjs_d, int n_sparse, double *x_d, double *V_d, 
						    double *dVdx2_d, double yp1, double ypn, double *xyz_d, double *forces_d, 
						    double *virial_d, double *energies_d, hipStream_t *stream){

  dim3 nblocks=dim3((i_end-i_beg+tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);

  kernel_get_core_pot<<<nblocks, nthreads,0,stream[0] >>>(i_beg, i_end, do_forces, species_d, sp1, sp2, n_neigh_d, neighbor_species_d,
                                                          rjs_d, n_sparse, x_d, V_d, dVdx2_d, yp1, ypn, xyz_d, forces_d, virial_d, 
							  energies_d);
}

extern "C" void gpu_device_sync()
{
  gpuErrchk( hipDeviceSynchronize() );
}
/* 
extern "C" void cuda_malloc_double(double **a_d, int Np)
{
   gpuErrchk(hipMalloc( (void **) a_d, sizeof(double) * Np ));
   return;
} */

/* extern "C" void cuda_malloc_double_complex(hipDoubleComplex **a_d, int Np)
{

  gpuErrchk(hipMalloc((void **)  a_d, sizeof(hipDoubleComplex) * Np));
  //gpuErrchk(hipMallocAsync( a_d, sizeof(hipDoubleComplex) * Np,0));
   return;
} */

/* extern "C" void cuda_malloc_int(int **a_d, int Np)
{
   // Allocate memory on GPU
   gpuErrchk(hipMalloc( (void **) a_d, sizeof(int) * Np ));
   return;
}
 */

/* extern "C" void cuda_malloc_bool(void **a_d, int Np)
{
   // Allocate memory on GPU
   gpuErrchk(hipMalloc( (void **) a_d, sizeof(bool) * Np ));
   return;
}
 */
