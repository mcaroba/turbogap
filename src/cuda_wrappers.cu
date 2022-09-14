 // wrappers file
// compile with:
// rm cuda_wrappers.o; nvcc -lcublas -lcurand -arch=sm_70 src/cuda_wrappers.cu -c;
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <curand.h>
#include <assert.h>
#include <cuComplex.h>

#define tpb 128

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
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__global__ void vect_dble(double *a, int N)
{
   int idx = threadIdx.x+blockIdx.x*gridDim.x;
   if (idx<N)printf(" %lf \n", a[idx]);
}

extern "C" void cuda_malloc_double(void **a_d, int Np)
{
   // Allocate memory on GPU
   //double **loc_a_d=(double **) a_d;
   //printf(" malloc double \n" );
   //if ( cudaSuccess != cudaMalloc( loc_a_d, sizeof(double) * Np )); exit(0);
   //gpuErrchk(cudaMalloc( a_d, sizeof(double) * Np ));
   //printf("Error in malloc double \n" );
   /*if ( cudaSuccess != cudaMalloc( a_d, sizeof(double) * Np ));
   printf("Error in malloc double \n" );*/
   cudaMallocAsync( a_d, sizeof(double) * Np,0);
   return;
}

extern "C" void cuda_malloc_double_complex(void **a_d, int Np)
{
   // Allocate memory on GPU
   //double **loc_a_d=(double **) a_d;
   //printf(" malloc double \n" );
   //if ( cudaSuccess != cudaMalloc( loc_a_d, sizeof(double) * Np )); exit(0);
   //gpuErrchk(cudaMalloc( a_d, sizeof(double) * Np ));
   //printf("Error in malloc double \n" );
   /*if ( cudaSuccess != cudaMalloc( a_d, sizeof(double) * Np ));
   printf("Error in malloc double \n" );*/
   cudaMallocAsync( a_d, sizeof(cuDoubleComplex) * Np,0);
   return;
}

extern "C" void cuda_malloc_int(int **a_d, int Np)
{
   // Allocate memory on GPU
   //gpuErrchk(cudaMalloc( a_d, sizeof(int) * Np ));
   cudaMallocAsync( a_d, sizeof(int) * Np, 0 );
   return;
}


extern "C" void cuda_free(double **a_d)
{
   cudaFree(*a_d);
   //printf("GPU memory freed \n");
   return;
}

extern "C" void GPU_fill_rand(double *A, int N, int ccc) {
	// Create a pseudo-random number generator
	curandGenerator_t prng;
	curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);

	// Set the seed for the random number generator using the system clock
	curandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock() + (unsigned long long)  ccc * N);

	// Fill the array with random numbers on the device
	curandGenerateUniformDouble(prng, A,N);
  //vect_dble<<<(N+128-1)/128,128>>>(A,N);
  //cudaDeviceSynchronize();
  printf("\n Filled \n");
}

extern "C" void cuda_cpy_double_htod(double *a, double *a_d, int N)
{
   //cudaMemcpy(a_d, a, sizeof(double) * N, cudaMemcpyHostToDevice);
   cudaMemcpyAsync(a_d, a, sizeof(double) * N, cudaMemcpyHostToDevice);


   return;
}

extern "C" void cuda_cpy_double_complex_htod(double *a, double *a_d, int N)
{
   //cudaMemcpy(a_d, a, sizeof(double) * N, cudaMemcpyHostToDevice);
   cudaMemcpyAsync(a_d, a, sizeof(cuDoubleComplex) * N, cudaMemcpyHostToDevice);


   return;
}



extern "C" void cuda_cpy_int_htod(int *a, int *a_d, int N)
{

   //cudaMemcpy(a_d, a, sizeof(int) * N, cudaMemcpyHostToDevice );
   cudaMemcpyAsync(a_d, a, sizeof(int) * N, cudaMemcpyHostToDevice );
   return;
}


extern "C" void cuda_cpy_double_dtoh(double *a_d, double *a ,int N)
{
  //cudaMemcpyAsync( a, a_d, sizeof(double) * N, cudaMemcpyDeviceToHost );
  cudaMemcpy( a, a_d, sizeof(double) * N, cudaMemcpyDeviceToHost );
  //gpuErrchk(cudaMemcpy( a, a_d, sizeof(double) * N, cudaMemcpyDeviceToHost ));
   //printf("\nTest cpy D to H \n");
   
   return;
}
extern "C" void cuda_cpy_double_dtod(double *b_d, double *a_d,int N)
{
   cudaMemcpyAsync( a_d, b_d, sizeof(double) * N, cudaMemcpyDeviceToDevice );
   //cudaMemcpy( a_d, b_d, sizeof(double) * N, cudaMemcpyDeviceToDevice );

   return;
}


extern "C" void create_cublas_handle(cublasHandle_t *handle)
{
 	 cublasCreate(handle);
   //printf("\n cublas handle created \n");

   return;
}

extern "C" void destroy_cublas_handle(cublasHandle_t *handle)
{
 	 // Destroy the handle
   cublasDestroy(*handle);
   //printf("\n cublas handle destroyed \n");

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

extern "C" void gpu_kernels_pow(double *a,double *b, double zeta, int size)
{
  int ntpb=256;
  int nblocks=(size+ntpb-1)/ntpb;
  gpu_pow<<<nblocks,ntpb>>>(a,b,zeta, size);
  //gpuErrchk( cudaPeekAtLastError() );
  //gpuErrchk( cudaDeviceSynchronize() );
  return;
}

extern "C" void gpu_blas_mmul_t_n(cublasHandle_t handle, const double *Qs_d, const double *soap_d, double *kernels_d, const int n_sparse, const int n_soap, const int n_sites)
//                                                           const double *A,     const double *B,         double *C,       const int nAx,
// const int nAy,      const int nBy,double *b, double zeta, int N)
{
// (cublasHandle_t handle, const double *Qs_d, const double *soap_d, double *kernels_d, const int n_sparse, const int n_soap, const int n_sites,double *b, double zeta, int N)
	const double alf = 1;
	const double bet = 0;

// soap(n_soap,n_sites)
// Qs(1:n_soap, 1:n_sparse)
// kernels(1:n_sites, 1:n_sparse)
// call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, kernels, n_sites)

	// Do the actual multiplication
  //printf("\n Pongo Longo \n");
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_sites, n_sparse, n_soap, &alf, soap_d, n_soap, Qs_d, n_soap, &bet, kernels_d, n_sites);

  return;
}

extern "C" void gpu_blas_mvmul_n(cublasHandle_t handle,  double *kernels_copy_d, const double *alphas_d, double *energies_d, const int n_sites, const int n_sparse)
{

	const double alf = 1;
	const double bet = 0;
	const double *alpha = &alf;
	const double *beta = &bet;

	// Do the actual multiplication
  cublasDgemv(handle, CUBLAS_OP_N, n_sites,n_sparse, alpha, kernels_copy_d, n_sites, alphas_d, 1, beta, energies_d, 1);
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

extern "C" void gpu_axpc(double *a, double dccc, double e0, int size)
{
  int ntpb=256;
  int nblocks=(size+ntpb-1)/ntpb;
  gpu_simpleaxpc<<<nblocks,ntpb>>>(a,dccc,e0, size);
  /*gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );*/
  return;

}

extern "C" void wrappers_all(double *soap, double *kernels, double *kernels_copy, double *Qs, double *energies, double delta, double zeta, double e0, int n_sites, int n_soap, int n_sparse, int size_kernels, int size_soap, int size_Qs, int size_alphas, int  size_energies)
{
  int ntpb=256;
  int nblocks=(size_kernels+ntpb-1)/ntpb;
  // Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasCreate(&handle);
  double *kernels_d, *kernels_copy_d, *soap_d, *Qs_d, *energies_d;
  cudaMalloc( &kernels_d, sizeof(double) * size_kernels );
  cudaMalloc( &kernels_copy_d, sizeof(double) * size_kernels );
  cudaMalloc( &soap_d, sizeof(double) * size_soap );
  cudaMalloc( &Qs_d, sizeof(double) * size_Qs );
  cudaMalloc( &energies_d, sizeof(double)*size_energies);


  const double alf = 1;
  const double bet = 0;

  cudaMemcpy(kernels_d, kernels, sizeof(double) * size_kernels, cudaMemcpyHostToDevice );
  cudaMemcpy(soap_d, soap, sizeof(double) * size_soap, cudaMemcpyHostToDevice );
  cudaMemcpy(Qs_d, Qs, sizeof(double) * size_Qs, cudaMemcpyHostToDevice );
  // Do the actual multiplication

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_sites, n_sparse, n_soap, &alf, soap_d, n_soap, Qs_d, n_soap, &bet, kernels_d, n_sites);
//cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N,  nBy, nAx, nAy, alpha, B, nAy, A, nAy, beta, C, nBy);
    //printf("\n cublasDgemm \n");
  // gpu_blas_mmul_t_n(cubhandle,     A,     B,      C,         nAx,      nAy,       nBy,             bb, zeta, N)
  // gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites, kernels_copy_d, zeta, size_kernels)

  cudaMemcpy( kernels, kernels_d, sizeof(double) * size_kernels, cudaMemcpyDeviceToHost );
  gpu_pow<<<dim3(nblocks,1,1),dim3(ntpb,1,1)>>>(kernels_d,kernels_copy_d, zeta, size_kernels);
  cudaMemcpy( kernels_copy, kernels_copy_d, sizeof(double) * size_kernels, cudaMemcpyDeviceToHost );
	// Destroy the handle
	cublasDestroy(handle);
  cudaFree(kernels_d);
  cudaFree(kernels_copy_d);
  cudaFree(soap_d);
  cudaFree(Qs_d);
  cudaFree(energies_d);
  //printf("\n %d %d %d %d %d %d %d %d  \n", n_sites, n_soap, n_sparse, size_kernels,  size_soap,  size_Qs,  size_alphas,  size_energies);
  //printf("\n %d %d %d\n", nblocks,ntpb, size_kernels);
  //exit(0);
 return;
}

extern "C" void cuda_set_device( int my_rank)
{

  int  num_gpus=0;
  gpuErrchk(cudaGetDeviceCount(&num_gpus));
  gpuErrchk(cudaSetDevice(my_rank%num_gpus));
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

extern "C" void cuda_matvect_kernels(double *kernels_d, double *alphas_d,int  n_sites, int n_sparse)
{
  int  ntpb=256;
  int nblocks=(n_sites*n_sparse+ntpb-1)/ntpb;
  matvect_kernels<<<nblocks,ntpb>>>(kernels_d,alphas_d,n_sites,n_sparse);
  /*gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );*/
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

extern "C" void cuda_matvect_qs(double *qs_d,double *qs_copy_d, double *alphas_d,int  n_soap, int n_sparse)
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
  matvect_qs<<<nblocks,ntpb>>>(qs_d,qs_copy_d,alphas_d,n_soap,n_sparse);
  /*gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );*/
  return;
}


// gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, n_soap, n_sites, cdelta)
extern "C" void gpu_blas_mmul_n_t(cublasHandle_t handle, const double *kernels_der_d, const double *Qs_copy_d, double *Qss_d, const int n_sparse, const int n_soap, const int n_sites, double cdelta)
{

	const double alf = cdelta;
	const double bet = 0;
	const double *alpha = &alf;
	const double *beta = &bet;
// soap(n_soap,n_sites)
// Qs(1:n_soap, 1:n_sparse)
// kernels(1:n_sites, 1:n_sparse)
// call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, kernels, n_sites)
// cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, n_sites, n_sparse, n_soap, alpha, soap_d, n_soap, Qs_d, n_soap, beta, kernels_d, n_sites);

// allocate( kernels_der(1:n_sites, 1:n_sparse)
// allocate( Qs_copy(1:n_soap, 1:n_sparse) ))
// allocate( Qss(1:n_sites, 1:n_soap) )
// call dgemm("n", "t", n_sites, n_soap, n_sparse, cdelta, kernels_der, n_sites, Qs_copy, n_soap, 0.d0, Qss, n_sites)
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, n_sites,  n_soap, n_sparse, alpha,  kernels_der_d, n_sites,Qs_copy_d, n_soap, beta, Qss_d, n_sites);
}

/*
!!$OMP parallel do private(i,j,l,j2,this_Qss)
!      l = 0
      do i = 1, n_sites
        this_Qss = Qss(i,1:n_soap)
!        do j = 1, n_neigh(i)
!          l = l + 1
         do l = neighbors_beg(i), neighbors_end(i)
          j2 = mod(neighbors_list(l)-1, n_sites0) + 1
          do k = 1, 3
            this_force(k) = dot_product(this_Qss, soap_der(k,:,l))
            forces(k, j2) = forces(k, j2) + this_force(k)
          end do
!         This is a many body potential, so there's no factor of 1/2 here
!          virial = virial + dot_product( this_force(1:3), xyz(1:3,l) )
          do k1 = 1, 3
            do k2 =1, 3
              virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,l) + this_force(k2)*xyz(k1,l))
            end do
          end do
        end do
      end do
!!$OMP end parallel do
*/
__global__ void cuda_soap_forces_virial(int *n_neigh_d,int n_sites,
                                        double *Qss_d,int n_soap,
                                        int *neighbors_beg_d,
                                        double3 *soap_der_d,
                                        double3 *xyz_d, double *virial_d,
                                        int *neighbors_list_d,int n_sites0, double *forces_d)
{
  int i_site=blockIdx.x;
  int i_nn=blockIdx.y;
  int n_nn_i_site=n_neigh_d[i_site];
  int tid=threadIdx.x;

  if(i_nn<n_nn_i_site)  // only the blocks with index smaller than the number of neighbors are participating in the calcultions
  {
    int l_nn=neighbors_beg_d[i_site]-1+i_nn;


    __shared__ double shxthis_block_force[tpb];
    __shared__ double shythis_block_force[tpb];
    __shared__ double shzthis_block_force[tpb];

    shxthis_block_force[tid]=0;
    shythis_block_force[tid]=0;
    shzthis_block_force[tid]=0;
    double locx_this_force=0;
    double locy_this_force=0;
    double locz_this_force=0;


    for(int ii=tid; ii< n_soap;ii=ii+tpb)
    {
        int i_Qss=i_site+ii*n_sites;
        double loc_this_Qss=Qss_d[i_Qss];// this rea<d  seems OK
        int in_soap_der=(l_nn*n_soap+ii);
        double3 loc_soap_der=soap_der_d[in_soap_der];
        locx_this_force+=loc_this_Qss*loc_soap_der.x;
        locy_this_force+=loc_this_Qss*loc_soap_der.y;
        locz_this_force+=loc_this_Qss*loc_soap_der.z;
      }

    shxthis_block_force[tid]=locx_this_force;
    shythis_block_force[tid]=locy_this_force;
    shzthis_block_force[tid]=locz_this_force;

    __syncthreads();

    //reduction
    for (int s=tpb/2; s>0; s>>=1) // s=s/2
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
      int j2=(neighbors_list_d[l_nn]-1) % (n_sites0);
      if(j2>= n_sites0)
      {printf("j2 the error! \n");}

      atomicAdd(&forces_d[j2*3]  , shxthis_block_force[0]);
      atomicAdd(&forces_d[j2*3+1], shythis_block_force[0]);
      atomicAdd(&forces_d[j2*3+2], shzthis_block_force[0]);

      // now the virial
      double this_force[3];
      this_force[0]=shxthis_block_force[0];
      this_force[1]=shythis_block_force[0];
      this_force[2]=shzthis_block_force[0];

      double3 tmp_xyz;
      tmp_xyz=xyz_d[l_nn];
      double this_xyz[3];
      this_xyz[0]=tmp_xyz.x;
      this_xyz[1]=tmp_xyz.y;
      this_xyz[2]=tmp_xyz.z;

      for(int k1=0;k1<3;k1++)
      {
        for(int k2=0;k2<3;k2++)
        {
          double loc_viri=0.5*(this_force[k1]*this_xyz[k2]+this_force[k2]*this_xyz[k1]);
          atomicAdd(&virial_d[k2+3*k1], loc_viri);
        }
      }
    }
  }
}

extern "C" void gpu_final_soap_forces_virial(int *n_neigh_d, int n_sites, int maxnn,
                                             double *Qss_d,int n_soap, int *neighbors_beg_d, double3 *soap_der_d,
                                             double3 *xyz_d, double *virial_d,
                                             int *neighbors_list_d,int n_sites0, double *forces_d)
{

     dim3 nblocks(n_sites,maxnn,1);
     cudaMemsetAsync(forces_d,0, 3*n_sites0*sizeof(double));
     cudaMemsetAsync(virial_d,0, 9*sizeof(double));
     cuda_soap_forces_virial<<<nblocks,tpb>>>(n_neigh_d,n_sites,
                                              Qss_d,n_soap, neighbors_beg_d,
                                              soap_der_d, xyz_d, virial_d,
                                              neighbors_list_d, n_sites0, forces_d);

     /*gpuErrchk( cudaPeekAtLastError() );
     gpuErrchk( cudaDeviceSynchronize() );*/

     return;
}


extern "C" void gpu_soap_energies_forces_virial(int *n_neigh_d, int n_sites, int maxnn,
                                             double *Qss_d,int n_soap, int *neighbors_beg_d,
                                             double3 *soap_der_d,
                                             double3 *xyz_d, double *virial_d,
                                             int *neighbors_list_d,int n_sites0, double *forces_d,
                                             cublasHandle_t handle, double *kernels_der_d, double *Qs_copy_d,
                                             const int n_sparse, double cdelta_force,
                                             double *alphas_d,
                                             double *kernels_d, double mzetam, int size_kernels,
                                             int do_forces,
                                             double *energies_d, double cdelta_ene, double e0, int size_energies,
                                             double *Qs_d, int size_Qs,
                                             double  *kernels_copy_d,
                                             double zeta,
                                             double *soap_d )
{
  gpu_blas_mmul_t_n(handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites);
  gpu_kernels_pow( kernels_d, kernels_copy_d,zeta, size_kernels);
  gpu_blas_mvmul_n(handle, kernels_copy_d, alphas_d, energies_d, n_sites, n_sparse);
  gpu_axpc( energies_d,cdelta_ene,e0, size_energies);
  if(do_forces==1)
  {
    cuda_cpy_double_dtod(Qs_d,   Qs_copy_d ,size_Qs);

    gpu_kernels_pow(kernels_d, kernels_der_d,mzetam, size_kernels);
    if(n_sites<n_soap)
    {
      cuda_matvect_kernels(kernels_der_d, alphas_d, n_sites, n_sparse);
    }
    else
    {
     cuda_matvect_kernels(Qs_copy_d, alphas_d, n_soap, n_sparse);
    }
       gpu_blas_mmul_n_t(handle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse,n_soap, n_sites, cdelta_force);

       gpu_final_soap_forces_virial(n_neigh_d, n_sites, maxnn,
                                     Qss_d, n_soap,neighbors_beg_d,
                                     soap_der_d, xyz_d, virial_d,
                                     neighbors_list_d, n_sites0, forces_d);
       /*dim3 nblocks(n_sites,maxnn,1);*/
       /*gpuErrchk(cudaMemset(forces_d,0, 3*n_sites0*sizeof(double)));
       gpuErrchk(cudaMemset(virial_d,0, 9*sizeof(double)));*/

       /*cudaMemsetAsync(forces_d,0, 3*n_sites0*sizeof(double));
       cudaMemsetAsync(virial_d,0, 9*sizeof(double));

       cuda_soap_forces_virial<<<nblocks,tpb>>>(n_neigh_d,n_sites,
                                            Qss_d,n_soap,neighbors_beg_d,
                                            soap_der_d,
                                            xyz_d, virial_d,
                                            neighbors_list_d, n_sites0, forces_d);*/
  }

     /*gpuErrchk( cudaPeekAtLastError() );
     gpuErrchk( cudaDeviceSynchronize() );*/
     return;
}



__global__ void cuda_get_soap_p(double *soap_d, double *multiplicity_array_d, 
                           cuDoubleComplex *cnk_d, int *skip_soap_component_d,
                           int n_sites, int n_soap, int n_max, int l_max)
{
   int i_site = threadIdx.x+blockIdx.x*blockDim.x;
   int k_max=1+l_max*(l_max+1)/2+l_max;
   if (i_site<n_sites){ 
    int counter=0;
    int counter2=0; 
    for(int n=0;n<n_max;n++){
      for(int np=n;np<n_max;np++){
        for(int l=0;l<=l_max;l++){
          if(skip_soap_component_d[l+(l_max+1)*(np+n*n_max)]!=0){
            counter++;
            for(int m=0;m<=l; m++){
              int k=1+l*(l+1)/2+m;
              counter2++;
              cuDoubleComplex tmp_cnk_d=cnk_d[k+k_max*(n+i_site*n_max)];
              soap_d[counter-1+i_site*n_soap]+=multiplicity_array_d[counter2-1]*tmp_cnk_d.x*tmp_cnk_d.y;
            }

          }
        }
      }
    }
 }
}

__global__ void cuda_get_sqrt_dot_p(double *soap_d, double *sqrt_dot_p_d,  
                                int n_sites, int n_soap)
{
  int i_site=blockIdx.x;
  int tid=threadIdx.x;
  __shared__ double sh_sqrt_dot_p_d[tpb];
  double this_dotprod=0.0;
  for(int s=0;s<n_soap;s=s+tpb){
    if(s<n_soap){
      this_dotprod+=soap_d[tid+s*tpb+i_site*n_soap]*soap_d[tid+s*tpb+i_site*n_soap];
    } 
  }
  sh_sqrt_dot_p_d[tid]=this_dotprod;
  __syncthreads();

  //reduction
  for (int s=tpb/2; s>0; s>>=1) // s=s/2
  {
    if (tid < s)
    {
      sh_sqrt_dot_p_d[tid] +=sh_sqrt_dot_p_d[tid + s];
    }
    __syncthreads();

  }
  //  at this point this_force is computed
  if(tid==0){
    double final_dotprod=1.0;
    if(sh_sqrt_dot_p_d[0]>=1.0e-5){
      final_dotprod=sh_sqrt_dot_p_d[0];
    }
    sqrt_dot_p_d[i_site]=final_dotprod;
  }
}

extern "C" void gpu_get_sqrt_dot_p(double *sqrt_dot_d, double *soap_d, double *multiplicity_array_d, 
                                   cuDoubleComplex *cnk_d, int *skip_soap_component_d, 
                                   int n_sites, int n_soap, int n_max, int l_max)
{
  dim3 nblocks=dim3((n_sites+1-tpb)/tpb,1,1);
  dim3 nthreads=dim3(tpb,1,1);
  cuda_get_soap_p<<<nblocks, nthreads>>>(soap_d, multiplicity_array_d, cnk_d, skip_soap_component_d, 
                                         n_sites, n_soap, n_max, l_max);
  cuda_get_sqrt_dot_p<<<n_sites,tpb>>>(soap_d,sqrt_dot_d, n_sites, n_soap);                                      
  return;
}