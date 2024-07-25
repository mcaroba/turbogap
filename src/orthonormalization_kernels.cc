#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <hip/hip_runtime.h>
#include <hip/hip_runtime.h>
#include <hipblas.h>
#include <hipsolver.h>
#include <chrono>
#include <numbers>
#include <cstring>

//first version of kernels, no stream no gpu pointers in/out


#include <limits>	
template<bool is_gauss>
__global__ void init_S (int alpha_max, double* S)
{
        //fortran/c row column major should not matter since is specular across diagonal? it matters for coalesced but if i flip here i/j should not change result. however thats not true for following pieces of kernel so i should already setup for flipped row/col
        //
        //i -> row -> far addresses -> divisio
        //j -> col -> closest addresses -> module
//        printf("test %p\n",S);
        auto tid = blockIdx.x * blockDim.x + threadIdx.x;
        auto i = tid / alpha_max;
        auto j = tid % alpha_max;
        int end_value;
        if constexpr (is_gauss)
          end_value = alpha_max -1;
        else
          end_value = alpha_max;
        if (i < end_value  && j < end_value )
        {

                if (i == j)
                {
                        S[tid] = 1.0;
                }
                else
                {
			//to get same number i need to add 1 to both i and j because fortran loop start from 1 and c starts from 0
			auto value1 = static_cast<double>(5+2*(i+1)) * static_cast<double>(5+2*(j+1));
			auto value2 = static_cast<double>(5+i+j+2);
			auto v1 = sqrt(value1);
			auto result= v1/value2;
			S[tid] =result;
                }
        }
	//initialize leftover to 0
	else if (tid < alpha_max * alpha_max)
	{
		S[tid] = 0.0;
	}
}


//this copy is only for the turbogap get_orthonormalization_matrix_poly3gauss, it copies the missing triangle of the symmetric matrix, inserts the diagonal of "1"s and than copies the resulting S matrix into W, provided that enough threads are used.
__global__ void orthonormalization_copy_matrix(double* S, double* edge_S, double* W, const int alpha_max)
{
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  auto i = tid / alpha_max;
  auto j = tid % alpha_max;
  if(i == alpha_max -1 && j == alpha_max -1)
  {
    S[tid]= 1.0;
  }
  else if(i == alpha_max -1 )
  {
    S[tid] = edge_S[j];
  }

  else if ( j == alpha_max-1)
  {
    S[tid] = edge_S[i];
  }
  
  W[tid] = S[tid];
}


inline double N_a (double rcut, int a)
{

    const int b = 2*a + 5;
    return sqrt( rcut / static_cast<double>(b) );

}

//creates a diagonal matrix with as values the square roots of the values provided by diagS parameter, provided that enough threads are used.
__global__ void sqr_diag_mat(double* S, double* diagS, const int alpha_max)
{
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  auto i = tid / alpha_max;
  auto j = tid % alpha_max;
  if(i == j)
    S[tid] = sqrt(diagS[i]);
  else 
    S[tid] = 0;

}

//copies the missing triangular matrix on the other half, to create a symmetric matrix, provided that enough threads are used.
__global__ void mirror_diag_mat(double* S, const int alpha_max)
{
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  auto i = tid / alpha_max;
  auto j = tid % alpha_max;
  if(j > i)
    S[tid] = S[j*alpha_max +i];

}

//creates an Id matrix of order alpha_max, provided that enough threads are used.
__global__ void init_id_mat(double* S, const int alpha_max)
{
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  auto i = tid / alpha_max;
  auto j = tid % alpha_max;
  if(j == i)
    S[tid] = 1;
  else
    S[tid] = 0;
}
void sqrt_mat(double* W, const int alpha_max, hipStream_t& s, hipblasHandle_t& handle )
{
//things to get W^0.5: 
//step 1: dgesvd
  double* U;
  double* VT;
  double* svd;
  double * work;
  int* devinfo;
  hipMalloc((void**) &U, alpha_max*alpha_max*sizeof(double));
  hipMalloc((void**) &VT, alpha_max*alpha_max*sizeof(double));
  hipMalloc((void**) &svd, alpha_max*sizeof(double));
  hipMalloc((void**) &devinfo, 1*sizeof(int));
  int lworksize; 
  { 
     auto status = hipsolverDgesvd_bufferSize(handle, 'A', 'A', alpha_max, alpha_max, &lworksize);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR in hipsolverDgesvd_bufferSize, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipMalloc((void**) &work, lworksize*sizeof(double));
  

  { 
    auto status =  hipsolverDnDgesvd( handle, 'A', 'A', alpha_max, alpha_max, W, alpha_max, svd, U, alpha_max, VT, alpha_max, work,  lworksize, nullptr, devinfo);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipFree(work);
//step2: diag mat with sqroots of svd diagonal
  auto total_threads= alpha_max*alpha_max;
  auto thread_per_block = 64;
  auto num_blocks = total_threads/thread_per_block >= 1 ? total_threads/thread_per_block : 1 ;
  sqr_diag_mat<<<num_blocks,thread_per_block,0,s>>>(W,svd,alpha_max);
//step3: matmuls
  double* tmpmat;
  hipMalloc((void**) &tmpmat, alpha_max*alpha_max*sizeof(double));
  {
    //U*W
    auto alpha = 1.0;
    auto beta = 0.0;
    auto status = hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, alpha_max, alpha_max, alpha_max, &alpha, U, alpha_max, W, alpha_max, &beta, tmpmat, alpha_max);
    if (status != HIPBLAS_STATUS_SUCCESS) {
      std::cout << "ERROR IN hipblasDgemm, status is" << status << std::endl;
      throw std::runtime_error("");
    }
  }
  {
    //W*VT
    auto alpha = 1.0;
    auto beta = 0.0;
    auto status = hipblasDgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, alpha_max, alpha_max, alpha_max, &alpha, tmpmat, alpha_max, VT, alpha_max, &beta, W, alpha_max);
    if (status != HIPBLAS_STATUS_SUCCESS) {
      std::cout << "ERROR IN hipblasDgemm, status is" << status << std::endl;
      throw std::runtime_error("");
    }
  }
  hipFree(tmpmat);

}


void invert_mat_other(double* W, const int alpha_max, hipStream_t& s, hipblasHandle_t& handle )
{
  double * work;
  int lworksize;
  int* devinfo;
// printf("starting the invert mat other function \n");
  hipMalloc((void**) &devinfo, 1*sizeof(int));
  std::cout<<std::flush;

  {
    auto status =  hipsolverDpotrf_bufferSize ( handle, HIPSOLVER_FILL_MODE_UPPER, alpha_max, W, alpha_max,&lworksize);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }

  hipMalloc((void**) &work, lworksize*sizeof(double));
  {
    auto status =  hipsolverDpotrf( handle, HIPSOLVER_FILL_MODE_UPPER, alpha_max, W, alpha_max,work,lworksize, devinfo);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipFree(work);

 
  {
    auto status =  hipsolverDpotri_bufferSize ( handle, HIPSOLVER_FILL_MODE_UPPER, alpha_max, W, alpha_max,&lworksize);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipMalloc((void**) &work, lworksize*sizeof(double));
  { 
    auto status =  hipsolverDpotri( handle, HIPSOLVER_FILL_MODE_UPPER, alpha_max, W, alpha_max,work,lworksize, devinfo);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipFree(work);
  auto total_threads= alpha_max*alpha_max;
  auto thread_per_block = 64;
  auto num_blocks = total_threads/thread_per_block >= 1 ? total_threads/thread_per_block : 1 ;
  mirror_diag_mat<<<num_blocks,thread_per_block,0,s>>>(W,alpha_max);

}
void invert_mat(double* &W, const int alpha_max, hipStream_t& s, hipblasHandle_t& handle )
{
  double * work;
  double * IdMat;
  int * Ipivot;
  int lworksize;
  int* devinfo;
  hipMalloc((void**) &IdMat, alpha_max*alpha_max*sizeof(double));
  hipMalloc((void**) &Ipivot, alpha_max*sizeof(int));
  auto total_threads= alpha_max*alpha_max;
  auto thread_per_block = 64;
  auto num_blocks = total_threads/thread_per_block >= 1 ? total_threads/thread_per_block : 1 ;
  init_id_mat<<<num_blocks,thread_per_block,0,s>>>(IdMat,alpha_max);

//hipsolverStatus_t hipsolverDgetrf_bufferSize(hipsolverHandle_t handle, int m, int n, double *A, int lda, int *lwork)

  {
    auto status =  hipsolverDgetrf_bufferSize( handle, alpha_max, alpha_max, W, alpha_max,&lworksize);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipMalloc((void**) &work, lworksize*sizeof(double));
//hipsolverStatus_t hipsolverDgetrf(hipsolverHandle_t handle, int m, int n, double *A, int lda, double *work, int lwork, int *devIpiv, int *devInfo)  

  {
    auto status =  hipsolverDgetrf( handle, alpha_max, alpha_max, W, alpha_max,work,lworksize,Ipivot ,devinfo);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipFree(work);
//hipsolverStatus_t hipsolverDgetrs_bufferSize(hipsolverHandle_t handle, hipsolverOperation_t trans, int n, int nrhs, double *A, int lda, int *devIpiv, double *B, int ldb, int *lwork)

  {
    auto status =  hipsolverDgetrs_bufferSize( handle, HIPSOLVER_OP_N, alpha_max, alpha_max, W,alpha_max ,Ipivot, IdMat, alpha_max,&lworksize);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
    hipMalloc((void**) &work, lworksize*sizeof(double));
  } 
  {
//hipsolverStatus_t hipsolverDgetrs(hipsolverHandle_t handle, hipsolverOperation_t trans, int n, int nrhs, double *A, int lda, int *devIpiv, double *B, int ldb, double *work, int lwork, int *devInfo)
    auto status =  hipsolverDgetrs( handle,HIPSOLVER_OP_N,alpha_max, alpha_max, W, alpha_max,Ipivot,IdMat,alpha_max,work, lworksize ,devinfo);
    if(status != HIPSOLVER_STATUS_SUCCESS)
    {
      printf ("\n\n ERROR, status is: %d \n\n", status);
      throw std::runtime_error("");
    }
  }
  hipFree(work);
  hipFree(Ipivot);
  //note that after hipsolverDgetrs, which solve system AX=B, results are in B. in this case B was the identity matrix so our system is AA^-1 = I, so the results are stored in the IdMat value on the gpu.
  //I don't like this but if i want to hide complexity is kinda needed. N.B. matrix W is passed as ref to pointer for this reason. and I don't really want to do a memcpy here, since "old" W is no more needed
  hipFree(W);
  W = IdMat;
}



//specific for the orthonormals: it will copy on the diagonal of the supermatrix	
__global__ void orthonormalization_copy_to_global_mats(double* src, double* dest, int alpha_max, int d_beg, int d_colsize)
{
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  auto i = tid / alpha_max;
  auto j = tid % alpha_max;
  if (tid < alpha_max*alpha_max)
  {
    dest[i*d_colsize+j+d_colsize*d_beg+d_beg] = src[tid];
  }

}



































// actual kernels in cpp

void get_orthonormalization_matrix_poly3gauss_cc(const int alpha_max, const double atom_sigma_in, const double rcut_hard_in, double *S, double* W, hipblasHandle_t handle, hipStream_t s){

  //local data(?)
  double * edge_S;
  double * edge_SH;

  hipMalloc((void**) &edge_S, alpha_max*sizeof(double));
  hipHostMalloc((void**)&edge_SH, alpha_max*sizeof(double));
//  hipDeviceSynchronize();
//  hipError_t err = hipGetLastError();
//  if (err != hipSuccess) {
//    printf("CUDA error: %s\n", hipGetErrorString(err));
//  fflush(stdout); 
//} 
  
  auto total_threads= alpha_max*alpha_max;
  auto thread_per_block = 64;
  auto num_blocks = total_threads/thread_per_block >= 1 ? total_threads/thread_per_block : 1 ;
//  printf("begin of orthonormalization poly3gauss: before init s\n");
//  fflush(stdout); 
  init_S<true><<<num_blocks,thread_per_block,0,s>>>(alpha_max,S);

  hipDeviceSynchronize();
//  printf("begin of orthonormalization poly3gauss: after init s \n");
//  fflush(stdout); 
//  err = hipGetLastError();
//  if (err != hipSuccess) {
//    printf("CUDA error in init_S: %s\n", hipGetErrorString(err));
//  fflush(stdout); 
//} 
  {
    const auto atom_sigma =  atom_sigma_in/rcut_hard_in;
    const auto rcut_hard = 1.0;
    const auto s2 = atom_sigma * atom_sigma;
    const auto sq2 = sqrt(2.0);
    auto pi = std::numbers::pi_v<double>;  
    auto I_n = 0.0;
    auto N_n = 1.0;
    auto N_np1 = N_a(rcut_hard, -2);
    auto I_np1 = sqrt(pi/2.0) * atom_sigma * erf( rcut_hard/sq2/atom_sigma ) / N_np1;
    auto C2 = s2 / rcut_hard;
 
 
    //loop here is completely sequential, doesnt really make sense to gpu this... cant even diagonalize and shfl. guess best option is to go for parallel cpu execution and then copy results in gpu matrix
    for (auto i=-1; i< alpha_max; ++i)
    {
      C2 = C2 * rcut_hard;
      auto N_np2 = N_a(rcut_hard, i);
      auto I_np2 = s2 * static_cast<double>(i+1) * N_n/ N_np2 * I_n + N_np1 * rcut_hard / N_np2 * I_np1 - C2 / N_np2; 
      if(i > 0)
      {
              //Include the normalization factor of the Gaussian
              edge_SH[i-1] = I_np2 * sq2 / sqrt(atom_sigma) / std::pow(pi,0.25);
      }
      N_n = N_np1;
      N_np1 = N_np2;
      I_n = I_np1;
      I_np1 = I_np2;
    }
  }

  hipMemcpy(edge_S,edge_SH,alpha_max*sizeof(double),hipMemcpyHostToDevice);
//  hipDeviceSynchronize();
//  err = hipGetLastError();
//  if (err != hipSuccess) {
//    printf("CUDA error in memcpy: %s\n", hipGetErrorString(err));
//  fflush(stdout); 
//  } 
  //finalize S matrix, and initialize the W matrix at the same time
  orthonormalization_copy_matrix<<<num_blocks,thread_per_block,0,s>>>(S,edge_S,W,alpha_max);
//  hipDeviceSynchronize();
//  err = hipGetLastError();
//  if (err != hipSuccess) {
//    printf("CUDA error in copy matrix: %s\n", hipGetErrorString(err));
//  fflush(stdout); 
//  } 

  sqrt_mat(W, alpha_max, s, handle);
//  hipDeviceSynchronize();
//  err = hipGetLastError();
//  if (err != hipSuccess) {
//    printf("CUDA error in sqrt mat: %s\n", hipGetErrorString(err));
//  fflush(stdout); 
//  } 
  //invert_mat(W, alpha_max, s, handle);
  invert_mat_other(W, alpha_max, s, handle);
//  hipDeviceSynchronize();
//  err = hipGetLastError();
//  if (err != hipSuccess) {
//    printf("CUDA error in invert mat: %s\n", hipGetErrorString(err));
//  fflush(stdout); 
//  } 



//
//  printf (" after invert matrix, W \n\n");
//  hipMemcpy(HW,W,alpha_max*alpha_max*sizeof(double),hipMemcpyDeviceToHost);
//  for (int i=0; i<alpha_max*alpha_max; i++)
//    printf("%.17g ",HW[i]);
// 
//  printf("\n");
//  std::memcpy(W_inout, HW, alpha_max*alpha_max*sizeof(double));
//
//  printf (" after invert matrix, S \n\n");
//  hipMemcpy(HS,S,alpha_max*alpha_max*sizeof(double),hipMemcpyDeviceToHost);
//  for (int i=0; i<alpha_max*alpha_max; i++)
//    printf("%.17g ",HS[i]);
// 
//  printf("\n");
//  std::memcpy(S_inout, HS, alpha_max*alpha_max*sizeof(double));

}


void get_orthonormalization_matrix_poly3_cc(const int alpha_max, double* S, double* W, hipblasHandle_t handle, hipStream_t s){
  

  auto total_threads= alpha_max*alpha_max;
  auto thread_per_block = 64;
  auto num_blocks = total_threads/thread_per_block >= 1 ? total_threads/thread_per_block : 1 ;
  init_S<false><<<num_blocks,thread_per_block,0,s>>>(alpha_max,S);
  //cpy S to W
  hipMemcpy(W,S,alpha_max*alpha_max*sizeof(double),hipMemcpyDeviceToDevice);
  
  sqrt_mat(W, alpha_max, s, handle);
  //invert_mat(W, alpha_max, s, handle);
  invert_mat_other(W, alpha_max, s, handle);




//  printf (" after invert matrix, W \n\n");
//  hipMemcpy(HW,W,alpha_max*alpha_max*sizeof(double),hipMemcpyDeviceToHost);
//  for (int i=0; i<alpha_max*alpha_max; i++)
//    printf("%.17g ",HW[i]);
// 
//  printf("\n");
//  std::memcpy(W_inout, HW, alpha_max*alpha_max*sizeof(double));
//  printf (" after invert matrix, S \n\n");
//  hipMemcpy(HS,S,alpha_max*alpha_max*sizeof(double),hipMemcpyDeviceToHost);
//  for (int i=0; i<alpha_max*alpha_max; i++)
//    printf("%.17g ",HS[i]);
// 
//  printf("\n");
//  std::memcpy(S_inout, HS, alpha_max*alpha_max*sizeof(double));
//
}

void orthonormalization_copy_to_global_matrix_cc(double* src, double* dest, const int src_rowsize, const int dest_start, const int dest_rowsize,hipStream_t s ) 
{
  auto total_threads= src_rowsize*src_rowsize;
  auto thread_per_block = 64;
  auto num_blocks = total_threads/thread_per_block >= 1 ? total_threads/thread_per_block : 1 ;
  orthonormalization_copy_to_global_mats<<<num_blocks,thread_per_block,0,s>>>(src,dest,src_rowsize,dest_start,dest_rowsize);
}




//exposed C APIs 
//stream and handle are passed as "ptr_c" type in fortran so they need to be dereferenced.
extern "C"
{
  void get_orthonormalization_matrix_poly3(const int alpha_max, double* S, double* W,hipblasHandle_t* handle, hipStream_t *s){
    get_orthonormalization_matrix_poly3_cc(alpha_max, S, W,*handle,*s);
  }
  void get_orthonormalization_matrix_poly3gauss(const int alpha_max, const double atom_sigma_in, const double rcut_hard_in, double *S, double* W,hipblasHandle_t* handle, hipStream_t* s){
//    printf("calling ortho from extern c, pointers are W%p and S%p\n",W,S);
//    fflush(stdout); 

    get_orthonormalization_matrix_poly3gauss_cc(alpha_max, atom_sigma_in, rcut_hard_in, S, W,*handle,*s);
  }

  void copy_to_global_matrix(double* src, double* dest, const int src_rowsize, const int dest_start, const int dest_rowsize,hipStream_t* s ) {
//    printf("copy to matrix begin, values are %d, %d, %d \n",src_rowsize,dest_start,dest_rowsize);
//    fflush(stdout); 
    orthonormalization_copy_to_global_matrix_cc(src, dest, src_rowsize, dest_start,  dest_rowsize, *s ); 
//    printf("copy to matrix end\n");
//    fflush(stdout); 
  }
}
