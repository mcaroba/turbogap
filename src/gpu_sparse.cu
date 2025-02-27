#include <cuda_runtime_api.h> // cudaMalloc, cudaMemcpy, etc.
#include <cusparse.h>         // cusparseSpMV
#include <stdio.h>            // printf
#include <stdlib.h>           // EXIT_FAILURE

#define CHECK_CUDA(func)                                                       \
{                                                                              \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
        printf("CUDA API failed at line %d with error: %s (%d)\n",             \
               __LINE__, cudaGetErrorString(status), status);                  \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(status), status);              \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

__global__ void print_something(
                int A_num_rows,
                int A_num_cols,
                int A_nnz,
                double* dA_values,
                int* dA_columns,
                int* dA_csrOffsets,
                int* ipointer_d,
                double* dY,
                double* dX){

        for (int i=0;i<9;i++){
        printf("\n dA in C %lf ", dA_values[i]);
}

}


extern "C" void gpu_sparse_matrix_mul_kernel( 
		int A_num_rows,
	       	int A_num_cols, 
		int A_nnz,
		double* dA_values, 
		int* dA_columns, 
		int* dA_csrOffsets, 
		int* ipointer_d, 
		double* dY, 
		double* dX, 
		cudaStream_t gpu_stream, cusparseHandle_t handle ){

   	cusparseSpMatDescr_t matA;
   	cusparseDnVecDescr_t vecX, vecY;
   	void*                dBuffer = NULL;
   	size_t               bufferSize = 0;
   	// Perform sparse matrix-vector multiplication using cuSPARSE
   	const double alpha = 1.0;  // Scalar multiplier for A*x
   	const double beta = 0.0;   // Scalar multiplier for y
   	// const int  A_num_cols =  A_num_rows;

	cusparseStatus_t  status ;
	// Create sparse matrix A in CSR format
	status =    cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_nnz,
			dA_csrOffsets, dA_columns, dA_values,
			CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
			CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) ;

//    double     hA_values[]     = { -1.0f, -2.0f, 3.0f, 4.0f, 5.0f,
//                                  6.0f, -7.0f, 8.0f, 9.0f };
//cudaMemcpy(hA_values, dA_values, A_nnz * sizeof(double),
//                           cudaMemcpyDeviceToHost);

//        for (int i=0;i<9;i++){
//        printf("\n hA in C %lf ", hA_values[i]);
//}

//	printf("Something");
//        print_something<<<1,1>>>(A_num_rows,
//                A_num_cols, 
//                A_nnz,
//                dA_values,
//                dA_columns,
//                dA_csrOffsets,
//                ipointer_d,
//                dY,
//                dX);

//	cudaDeviceSynchronize();
//	exit(0);

	// Create dense vector X
 	status =  cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_64F)  ;
	// Create dense vector y
 	status =  cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_64F)  ;
	// allocate an external buffer if needed

//        cudaDeviceSynchronize();
//        exit(0);



 	status =  cusparseSpMV_bufferSize(
			handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			&alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
			CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize)  ;

	//  	CHECK_CUSPARSE( cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_nnz,
	//  				dA_csrOffsets, dA_columns, dA_values,
	//  				CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
	//  				CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )
	//      // Create dense vector X
	//      CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_64F) )
	//      // Create dense vector y
	//      CHECK_CUSPARSE( cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_64F) )
	//      // allocate an external buffer if needed
	//      CHECK_CUSPARSE( cusparseSpMV_bufferSize(
	//                                   handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
	//                                   &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
	//                                   CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) )
	//      CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) )
	//  
	//    // execute preprocess (optional)
	//    CHECK_CUSPARSE( cusparseSpMV_preprocess(
	//                                 handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
	//                                 &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
	//                                 CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) )

	// execute SpMV
        cudaMalloc(&dBuffer, bufferSize);

 	status =  cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			&alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
			CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);

	// destroy matrix/vector descriptors
	status =  cusparseDestroySpMat(matA)  ;
	status =  cusparseDestroyDnVec(vecX)  ;
	status =  cusparseDestroyDnVec(vecY)  ;
	
	cudaFree(&dBuffer) ;
}



int main_blah(void) {
    // Host problem definition
    const int A_num_rows      = 4;
    const int A_num_cols      = 4;
    const int A_nnz           = 9;
    int       hA_csrOffsets[] = { 0, 3, 4, 7, 9 };
    int       hA_columns[]    = { 0, 2, 3, 1, 0, 2, 3, 1, 3 };
    double     hA_values[]     = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f,
                                  6.0f, 7.0f, 8.0f, 9.0f };
    double     hX[]            = { 1.0f, 2.0f, 3.0f, 4.0f };
    double     hY[]            = { 0.0f, 0.0f, 0.0f, 0.0f };
    double     hY_result[]     = { 19.0f, 8.0f, 51.0f, 52.0f };
    double     alpha           = 1.0f;
    double     beta            = 0.0f;
    //--------------------------------------------------------------------------
    // Device memory management
    int   *dA_csrOffsets, *dA_columns;
    double *dA_values, *dX, *dY;
    CHECK_CUDA( cudaMalloc((void**) &dA_csrOffsets,
                           (A_num_rows + 1) * sizeof(int)) )
    CHECK_CUDA( cudaMalloc((void**) &dA_columns, A_nnz * sizeof(int))        )
    CHECK_CUDA( cudaMalloc((void**) &dA_values,  A_nnz * sizeof(double))      )
    CHECK_CUDA( cudaMalloc((void**) &dX,         A_num_cols * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &dY,         A_num_rows * sizeof(double)) )

    CHECK_CUDA( cudaMemcpy(dA_csrOffsets, hA_csrOffsets,
                           (A_num_rows + 1) * sizeof(int),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(dA_columns, hA_columns, A_nnz * sizeof(int),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(dA_values, hA_values, A_nnz * sizeof(double),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(dX, hX, A_num_cols * sizeof(double),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(dY, hY, A_num_rows * sizeof(double),
                           cudaMemcpyHostToDevice) )
    //--------------------------------------------------------------------------
    // CUSPARSE APIs
    cusparseHandle_t     handle = NULL;
    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vecX, vecY;
    void*                dBuffer    = NULL;
    size_t               bufferSize = 0;
    CHECK_CUSPARSE( cusparseCreate(&handle) )
    // Create sparse matrix A in CSR format
    CHECK_CUSPARSE( cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_nnz,
                                      dA_csrOffsets, dA_columns, dA_values,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )
    // Create dense vector X
    CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_64F) )
    // Create dense vector y
    CHECK_CUSPARSE( cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_64F) )
    // allocate an external buffer if needed
    CHECK_CUSPARSE( cusparseSpMV_bufferSize(
                                 handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                                 CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) )
    CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) )

//     // execute preprocess (optional)
//     CHECK_CUSPARSE( cusparseSpMV_preprocess(
//                                  handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
//                                  &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
//                                  CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) )

    // execute SpMV
    CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                                 CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) )

    // destroy matrix/vector descriptors
    CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(vecX) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(vecY) )
    CHECK_CUSPARSE( cusparseDestroy(handle) )
    //--------------------------------------------------------------------------
    // device result check
    CHECK_CUDA( cudaMemcpy(hY, dY, A_num_rows * sizeof(double),
                           cudaMemcpyDeviceToHost) )
    int correct = 1;
    for (int i = 0; i < A_num_rows; i++) {
        if (hY[i] != hY_result[i]) { // direct floating point comparison is not
            correct = 0;             // reliable
            break;
        }
    }
    if (correct)
        printf("spmv_csr_example test PASSED\n");
    else
        printf("spmv_csr_example test FAILED: wrong result\n");
    //--------------------------------------------------------------------------
    // device memory deallocation
    CHECK_CUDA( cudaFree(dBuffer) )
    CHECK_CUDA( cudaFree(dA_csrOffsets) )
    CHECK_CUDA( cudaFree(dA_columns) )
    CHECK_CUDA( cudaFree(dA_values) )
    CHECK_CUDA( cudaFree(dX) )
    CHECK_CUDA( cudaFree(dY) )
    return EXIT_SUCCESS;
}




extern "C" void create_cusparse_handle(cusparseHandle_t *handle,cudaStream_t *stream )
{
 	  cusparseCreate(handle);
    //cudaStreamCreate(stream);
    //cusparseSetStream(*handle, *stream);
//    hipsolverCreate(handle);
//    hipsolverSetStream(*handle,*stream);
    /*printf("\n cusparse handle created \n");
    exit(0);*/

   return;
}

extern "C" void destroy_cusparse_handle(cusparseHandle_t *handle,cudaStream_t *stream )
{
 	 // Destroy the handle
   cusparseDestroy(*handle);
   //printf("\n cusparse handle destroyed. \n The End? \n");
   return;
}

