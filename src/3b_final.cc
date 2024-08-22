#include <hip/hip_runtime.h> 
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <hipblas/hipblas.h>
//#include <hipsolver.h>
//#include <hiprand/hiprand.h>
#include <assert.h>
#include <hip/hip_complex.h>
#include <numbers>
#include <iostream>
#include <iomanip>
#include <bit>
#include <chrono>

//stupid thing needed to have 0 starting exponential printing, debug only to compare with fortran exponential notation.

class Double {
public:
    Double(double x): value(x) {}
    const double value;
};

std::ostream & operator<< (std::ostream & stream, const Double & x) {
    // So that the log does not scream
    if (x.value == 0.) {
        stream << 0.;
        return stream;
    }

    int exponent = floor(log10(std::abs(x.value)));
    double base = x.value / pow(10, exponent);

    // Transform here
    //base /= 10;
    //exponent += 1;
    if (exponent < 0)
    stream << base << 'E' << std::setw(3) << std::setfill('0') << std::internal<<exponent; // Change the format as needed
    else 
    stream << base << "E+" << std::setw(2) << std::setfill('0') <<exponent; // Change the format as needed

    return stream;
}

enum kern_type {exp_k, pol_k};

void hip_check_error(hipError_t err)
{

  std::cout<< "reported error is " <<	hipGetErrorString(err) <<std::flush << std::endl;

}



// pref(1:n_sparse) = 2.d0 * delta**2 * alphas(1:n_sparse) * cutoff(1:n_sparse)
__global__ void eval_pref(double* __restrict__ pref, const double* __restrict__ alpha, const double delta, const double* __restrict__ cutoff, const int size){
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  // for (auto i = 0; i<= size/blockDim.x;++i)
  // {
  //   auto idx = tid + i*blockDim.x;
  //   if (idx < size)
  //     pref[idx] = 2.0 * ( delta * delta ) * alpha[idx] * cutoff[idx];
  // }
  
  auto idx = tid; // + i*blockDim.x;
  if (idx < size)
  {
      pref[idx] = 2.0 * ( delta * delta ) * alpha[idx] * cutoff[idx];
  }
}

__global__ void add_init_energy(double* __restrict__ arr, const double value, const int size, const int i_beg)
{
  auto tid = i_beg -1 + blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < size && tid > i_beg-2)
      arr[tid] += value;

}


__device__ double warp_red(double data) {

   double res = data;
   for (int i =warpSize/2; i!=0; i=i>>1) {
      res += __shfl_down_sync(0xffffffff,res, i,warpSize);
   }
   return res;
}

// __device__ double warp_red(double data) {

//    double res = data;
//    for (int i = 32; i!=0; i=i>>1) {
//       res += __shfl_down(res, i, 64);
//    }
//    return res;
// }

void setup_kappas(int* kappas_array, const int n_sites, const int* n_neigh)
{
  kappas_array[0]=0;
  for (int i=1; i<n_sites;++i)
  {
    kappas_array[i] = kappas_array[i-1]+n_neigh[i-1];
//    printf("for i %d, k starts at %d \n",i,kappas_array[i]);
//    fflush(stdout);
  }
//  printf("all kappas done \n");
//  fflush(stdout);
  
}
#if 1
//idea behind templating this is that i want to minimize ifs inside the kernel, even if they would be "branch taken" or "not taken" for the whole warp. this (should) also allow to allocate registers only when the variables are really needed (attention was given to scope of variables), thus allowing for having the needed things into registers and limiting "swap" between register and global memory
template< bool do_forces, kern_type type>
__device__ void eval_energies_d(
    const double rcut,
    const double buffer,
    const int n_sparse,
    const double* __restrict__ xyz, 
    const double* __restrict__ xyz_red, 
    const double* __restrict__ r, 
    const double* __restrict__ sigma,
    const double* __restrict__ qs, 
    const double* __restrict__ pref,
    double* energies, 
    double* forces, 
    double* virial, 
    double* fcut, 
    double* dfcut,
    const int i, 
    const int i3, 
    const int j3, 
    const int k3,
    const int n_sites0
    )
{
  double q[3];
  //const auto pi = std::numbers::pi_v<double>; 
  const auto pi =3.14159265358979323846264338327950288419716939; // std::numbers::pi_v<double>; 
  int tid= threadIdx.x;
  //setup some common, read only values into shmemory 
  if (tid < 6)
  {
    if (r[tid/3] < rcut - buffer)
    {
      if(tid%3 == 0) fcut[tid/ 3] = 1.0;  //need to assign only once, no bidimensional
      dfcut[tid] = 0.0;
    }
    else
    { 
      if(tid%3 == 0) fcut[tid/ 3] = ( cos( pi*(r[tid/3] - rcut + buffer) / buffer ) + 1.0 ) / 2.0; //need to assign only once, no bidimensional
      dfcut[tid] = sin( pi*(r[tid/3] - rcut + buffer) / buffer ) / 2.0 * pi / buffer * xyz_red[tid]; 
    }
  
  }
  __syncthreads();
  if(tid < 3)
  {
    if(tid%3 == 0) fcut[2] = fcut[0] * fcut[1];  // fcut = fcut 12*fcut 13 has been put in same array.
  }
  __syncthreads();
  if(tid < 3)
  {
    dfcut[6+tid] = fcut[0] * dfcut[3+tid] + fcut[1] * dfcut[tid];
  }
  __syncthreads();

  //setup some variables used more times in the loop. those could also be shmem tbh, never tested the difference. possible optimization here if needed.
  q[0] = r[0]+r[1];
  q[1] = (r[0]-r[1])*(r[0]-r[1]) ;
  q[2] = r[2];

  auto sigma2_0 = sigma[0]*sigma[0];
  auto sigma2_1 = sigma[1]*sigma[1];
  auto sigma2_2 = sigma[2]*sigma[2];
  int offset=tid-warpSize;
  int nitera=(n_sparse+warpSize-1)/warpSize;
  //for(auto offset=tid; offset < n_sparse; offset+=warpSize)
  for(int  iter=0; iter< nitera; iter++)
  {
    offset+=warpSize;
    double local_r=1.0;
    if(offset<n_sparse)
    {
      local_r       =  pow(q[0] -qs[offset],2)            / sigma2_0;
      local_r      +=  pow(q[1] -qs[offset+n_sparse],2)   / sigma2_1;
      local_r      +=  pow(q[2] -qs[offset+2*n_sparse],2) / sigma2_2;
    }
      local_r      =  sqrt(local_r);
  
    if constexpr(type == pol_k)
    {
      constexpr int j_int = 3/2 + 1 + 1;
      constexpr double j = j_int;
      double my_energy=0.0;
      double energy_to_reduce;
      double covpp = 0.0;
      if(offset<n_sparse){
        //auto covpp = local_r >= 1 ? 0.0 : local_r <= 0 ? 1.0 : pow((1.0 - local_r),(j_int+1)) * ((j+1.0)*local_r + 1.0);
        covpp = local_r >= 1 ? 0.0 : local_r <= 0 ? 1.0 : pow((1.0 - local_r),(j_int+1)) * ((j+1.0)*local_r + 1.0);
        my_energy = pref[offset] * covpp ; 
      }
      energy_to_reduce = my_energy;
  
      //reduction of energy from warp to single value, syncthreads should never happen (left for sanity reasons)
      __syncthreads();
      double tmp = warp_red(energy_to_reduce);
      __syncthreads();
  
      if(tid == 0)
      {
        atomicAdd(&(energies[i]) , fcut[2]*tmp);
      }
  
      if constexpr (do_forces)
      {
        double drdq[3];
        drdq[0] = 0.0 ; drdq[1] = 0.0;  drdq[2] = 0.0;
        if(offset<n_sparse){
          drdq[0] = ( q[0] - qs[offset]           ) / sigma2_0 ;
          drdq[1] = ( q[1] - qs[offset+n_sparse]  ) / sigma2_1 ;
          drdq[2] = ( q[2] - qs[offset+2*n_sparse]) / sigma2_2 ;
        }
        double kernel_der =0.0;
        if(offset<n_sparse){
          auto cov_pp_der = local_r <= 0.0 || local_r >= 1.0 ? 0.0 : (( -(j+1.0) * pow((1.0 - local_r),j_int) * ((j+1.0)*local_r + 1.0) + pow((1.0 - local_r),(j_int+1)) * (j + 1.0) ) / local_r);
          kernel_der = pref[offset] * cov_pp_der ;
        }
        //split loop of drdx/forces into 3 loops, one for each force since drdx is only needed inside. can help reducing scope of variables? which means keeping things in registers and bla bla 
  #pragma unroll
        for (int l = 0; l < 3; ++l)
        {   
        // For atom 1
        // drdx1(1:n_sparse, l) = drdq(1:n_sparse, 1) * ( xyz12_red(l) + xyz13_red(l) ) + drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * ( xyz12_red(l) - xyz13_red(l) )
          auto drdx1 = drdq[0] * ( xyz_red[l] + xyz_red[l+3] ) + drdq[1] * 2.0 * (r[0]-r[1]) * ( xyz_red[l] - xyz_red[3+l] );
	  //force1(l) = - sum( kernel(1:n_sparse) * dfcut(l) - kernel_der(1:n_sparse) * fcut*drdx1(1:n_sparse, l) )
          auto force1 = my_energy * dfcut[6+l] - kernel_der * fcut[2] *drdx1;
          __syncthreads();
  	  auto tmp = warp_red(force1);
          __syncthreads();
          if(tid == 0)
          {
  	    atomicAdd(&(forces[i3*3+l]) , -tmp);
          }
        }//atom1, pol_k
  
  #pragma unroll
        for (int l = 0; l < 3; ++l)
        {   
          // For atom 2
          // drdx2(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz12_red(l) - drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * xyz12_red(l) + drdq(1:n_sparse, 3) * xyz23_red(l)
          auto drdx2 = - drdq[0] * xyz_red[l  ] - drdq[1] * 2.0 * (r[0]-r[1]) * xyz_red[l  ] + drdq[2] * xyz_red[l+6];
          // force2(l) = - sum( - kernel(1:n_sparse) * fcut13 * dfcut12(l) - kernel_der(1:n_sparse) * fcut*drdx2(1:n_sparse, l) )
          auto force2 = - my_energy * fcut[1] * dfcut[l] - kernel_der * fcut[2] *drdx2;
          __syncthreads();
  	  auto tmp = warp_red(force2);  
          __syncthreads();
          if(tid == 0)
          {
    	    atomicAdd(&(forces[j3*3+l] ), -tmp);
            #pragma unroll
  	    for (auto k=0; k<3; ++k)
  	    {
  	      virial[l*3+k] += 0.5 * -tmp/*l*/ * xyz[k];
  	      virial[k*3+l] += 0.5 * -tmp/*l*/ * xyz[k]; //this is not really the original formula, it should be k*3+l = force(k)*xyz(l) but i am exploiting the fact that virials are a symmetric matrix to speed up computations by flipping the indexes here and avoid to store forces(1:3) into variables -> less register usage/memory accesses -> higher efficiency
  	    }
            }
        }//atom2, pol_k
        #pragma unroll
        for (int l = 0; l < 3; ++l)
        {   
        //! For atom 3
        //  drdx3(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz13_red(l) - drdq(1:n_sparse, 2) * 2.d0 * (r13-r12) * xyz13_red(l) - drdq(1:n_sparse, 3) * xyz23_red(l)
          auto drdx3 = - drdq[0] * xyz_red[l+3] - drdq[1] * 2.0 * (r[1]-r[0]) * xyz_red[l+3] - drdq[2] * xyz_red[6+l];
  	__syncthreads();
	//                  force3(l) = - sum( - kernel(1:n_sparse) * fcut12 * dfcut13(l) - kernel_der(1:n_sparse) * fcut*drdx3(1:n_sparse, l) )
          auto force3 = - my_energy * fcut[0] * dfcut[3+l] - kernel_der * fcut[2] *drdx3;
          __syncthreads();
  	auto tmp = warp_red(force3);
          __syncthreads();
          if(tid == 0)
          {
  	    atomicAdd(&(forces[k3*3+l]) ,-tmp);
            #pragma unroll
	    for (auto k=0; k<3; ++k)
	    {
	      virial[l*3+k] += 0.5 * -tmp/*l*/ * xyz[3+k];
	      virial[k*3+l] += 0.5 * -tmp/*l*/ * xyz[3+k];//see above comment
	    }
          }
        }//atom3, pol_k
      }//do forces
    }//pol_k

    if constexpr(type == exp_k)
    { 
      printf("\n This part has a bug related to the __syncthreads(). \n Look above to see how to fix it\n");
      const auto my_energy = pref[offset] * exp(-0.5 * local_r * local_r) ;
      auto energy_to_reduce = my_energy;
      //reduction of energy from warp to single value, syncthreads should never happen (left for sanity reasons)
      __syncthreads();
      double tmp = warp_red(energy_to_reduce);
      __syncthreads();
  
      if(tid == 0)
      {
        atomicAdd(&(energies[i]) , fcut[2]*tmp);
      }
      if constexpr (do_forces)
      {
        double drdq[3];
        drdq[0] = ( q[0] - qs[offset]           ) / sigma2_0 ;
        drdq[1] = ( q[1] - qs[offset+n_sparse]  ) / sigma2_1 ;
        drdq[2] = ( q[2] - qs[offset+2*n_sparse]) / sigma2_2 ;
        #pragma unroll
        for (int l = 0; l < 3; ++l)
        {   
        // For atom 1
        // drdx1(1:n_sparse, l) = drdq(1:n_sparse, 1) * ( xyz12_red(l) + xyz13_red(l) ) + drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * ( xyz12_red(l) - xyz13_red(l) )
          auto drdx1 = drdq[0] * ( xyz_red[l] + xyz_red[l+3] ) + drdq[1] * 2.0 * (r[0]-r[1]) * ( xyz_red[l] - xyz_red[3+l] );
	  //force1(l) = - sum( kernel(1:n_sparse) * (dfcut(l) + fcut*drdx1(1:n_sparse, l)) )
	  auto force1 = my_energy * (dfcut[6+l] + fcut[2] *drdx1);
          __syncthreads();
  	  auto tmp = warp_red(force1);
          __syncthreads();
          if(tid == 0)
          {
  	    atomicAdd(&(forces[i3*3+l]) , -tmp);
          }
        }//atom1, exp_k
  
  #pragma unroll
        for (int l = 0; l < 3; ++l)
        {   
          // For atom 2
          // drdx2(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz12_red(l) - drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * xyz12_red(l) + drdq(1:n_sparse, 3) * xyz23_red(l)
          auto drdx2 = - drdq[0] * xyz_red[l  ] - drdq[1] * 2.0 * (r[0]-r[1]) * xyz_red[l  ] + drdq[2] * xyz_red[l+6];
	  // force2(l) = - sum( kernel(1:n_sparse) * (-fcut13 * dfcut12(l) + fcut*drdx2(1:n_sparse, l)) )
          auto force2 = my_energy *(- fcut[1] * dfcut[l] + fcut[2] *drdx2);
          __syncthreads();
  	  auto tmp = warp_red(force2);  
          __syncthreads();
          if(tid == 0)
          {
    	    atomicAdd(&(forces[j3*3+l] ), -tmp);
            #pragma unroll
  	    for (auto k=0; k<3; ++k)
  	    {
  	      virial[l*3+k] += 0.5 * -tmp/*l*/ * xyz[k];
  	      virial[k*3+l] += 0.5 * -tmp/*l*/ * xyz[k]; //this is not really the original formula, it should be k*3+l = force(k)*xyz(l) but i am exploiting the fact that virials are a symmetric matrix to speed up computations by flipping the indexes here and avoid to store forces(1:3) into variables -> less register usage/memory accesses -> higher efficiency
  	    }
            }
        }//atom2, exp_k
        #pragma unroll
        for (int l = 0; l < 3; ++l)
        {   
        //! For atom 3
        //  drdx3(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz13_red(l) - drdq(1:n_sparse, 2) * 2.d0 * (r13-r12) * xyz13_red(l) - drdq(1:n_sparse, 3) * xyz23_red(l)
          auto drdx3 = - drdq[0] * xyz_red[l+3] - drdq[1] * 2.0 * (r[1]-r[0]) * xyz_red[l+3] - drdq[2] * xyz_red[6+l];
  	__syncthreads();
	// force3(l) = - sum( - kernel(1:n_sparse) * fcut12 * dfcut13(l) - kernel_der(1:n_sparse) * fcut*drdx3(1:n_sparse, l) )
	// force3(l) = - sum( kernel(1:n_sparse) * (-fcut12 * dfcut13(l) + fcut*drdx3(1:n_sparse, l)) )
          auto force3 = my_energy * ( -fcut[0] * dfcut[3+l] + fcut[2] *drdx3);
          __syncthreads();
  	  auto tmp = warp_red(force3);
          __syncthreads();
          if(tid == 0)
          {
  	    atomicAdd(&(forces[k3*3+l]) ,-tmp);
            #pragma unroll
	    for (auto k=0; k<3; ++k)
	    {
	      virial[l*3+k] += 0.5 * -tmp/*l*/ * xyz[3+k];
	      virial[k*3+l] += 0.5 * -tmp/*l*/ * xyz[3+k];//see above comment
	    }
          }
        }//atom3, exp_k
      }//do_forces
    }//exp_k
  }//loop on n_sparse
}


#endif 

//grid is different "i" iterations (i.e. sites), block is inner loops.
template< bool do_forces, kern_type type>
__global__ void kernel_2nd_try(

    const int n_sparse, 
    const int n_sites, 
    const int n_atom_pairs, 
    const int n_sites0, 
    const int sp0,
    const int sp1,
    const int sp2,
    const double* alpha,
    const double delta,
    const double* cutoff,
    //host arrays
    const double* rjs,
    const double* xyz,
    const int* n_neigh,
    const int* species,
    const int* neighbors_list,
    const int* neighbor_species ,
    const int* kappas_array,
    const double rcut,
    const double buffer,
    //const double e0   ,
    //device array      
    const double* sigma,
    const double* qs,
    const double* pref_d,
    //output arrays
    double* energies,
    double* forces,
    double* virial,
    const int i_beg
    )
{
//    if (threadIdx.x == 0 && blockIdx.x == 0)
//  {
//    printf("################### inside kernel parameters are: nsparse %d, nsites %d, natompair %d, nsite0 %d \n",n_sparse,n_sites,n_atom_pairs,n_sites0);
//    printf("sp0 %d, sp1 %d, sp2 %d, alpha[0] %f, delta %f, cutoff[0] %f, rjs[0] %f, xyz[0] %f, n_neigh[0] %d,",sp0, sp1, sp2, alpha[0], delta, cutoff[0], rjs[0], xyz[0], n_neigh[0]);
//    printf("neighbors_list[0] %d, neighbor_species[0] %d, kappas_array_d[0] %d, rcut %f, buffer %f, ",neighbors_list[0], neighbor_species[0], kappas_array[0], rcut, buffer);
//    printf("sigma[0] %f, qs[0] %f, pref_d[0] %f, energies[0] %f, forces[0] %f, virial[0] %f \n",sigma[0], qs[0], pref_d[0], energies[0], forces[0], virial[0]);
//  }

  __shared__ double xyz_shmem[9];
  __shared__ double xyz_red_shmem[9];
  __shared__ double r_shmem[3];
  __shared__ double fcut[3];
  __shared__ double dfcut[9];
  __shared__ double virial_loc[9];
  
  int tid = threadIdx.x;
  int counter_first_cont =0;
  int counter_second_cont =0;
  int counter_launch =0;


  const int i = blockIdx.x +i_beg -1; //i_beg-1 + blockIdx.x
  if(species[i] != sp0)
  {
    printf("block %d has wrong species, returning \n",i);
    return;
  }
  if(tid==0)
  {
    #pragma unroll
    for(auto tmp=0; tmp<9; ++tmp)
    {
      virial_loc[tmp]=0;
    }
  }
  __syncthreads();
  //int k = i_beg != 1 ?  kappas_array[i-i_beg+1] : kappas_array[i] ;
  
   int k =  kappas_array[i];
  // int nblk = blockIdx.x;
  auto i3 = ((neighbors_list[k]-1)% n_sites0) ;//+ 1;
  //if(tid==0)
  //  printf ("hello from i %d, my k is %d, i3 is: %d \n", i, k, i3);

  //  if(tid==0) {printf("starting kernel, got k and i3\n");}
  for (int j = 1; j < n_neigh[i]; ++j)
  {
    k++;
    if(threadIdx.x == 0)
    {
      r_shmem[0] = rjs[k];
    }
    __syncthreads();
    if((neighbor_species[k] != sp1 && neighbor_species[k] != sp2) || r_shmem[0] > rcut )	  
    {
    // if(tid==0)printf("NOT launching kernel, loop variables are: i i%d j: %d and n_neigh[i] is: %d \n ",i,j,n_neigh[i]);
      counter_first_cont++;
      continue;
    }
    auto j3 =  ((neighbors_list[k]-1)%n_sites0);//+1;
    if (threadIdx.x < 3)
    {
      xyz_shmem[threadIdx.x] = xyz[k*3+threadIdx.x];
    }
    __syncthreads();
    if (threadIdx.x < 3)
    {
      xyz_red_shmem[threadIdx.x] = xyz_shmem[threadIdx.x]/r_shmem[0];
    }
    auto k2 = k;    
    __syncthreads();
    for(int j2=j+1; j2 < n_neigh[i]; ++j2)
    {
      k2++;
      if(threadIdx.x == 0)
      {
	r_shmem[1] = rjs[k2];
      }
      __syncthreads();

      if(!( ((neighbor_species[k] == sp1 && neighbor_species[k2] == sp2) ||
          (neighbor_species[k] == sp2 && neighbor_species[k2] == sp1) ) && r_shmem[1] < rcut ))
      {
        // if(tid==0)printf("NOT launching kernel, loop variables are: i %d j: %d and j2 %d and n_neigh[i] is: %d \n ",i,j,j2,n_neigh[i]);
        counter_second_cont++;
        continue;
      }
      auto k3 = ((neighbors_list[k2]-1)% n_sites0);// +1;
      if (threadIdx.x < 3)
      {
        xyz_shmem[3+threadIdx.x] = xyz[k2*3+threadIdx.x];}
    __syncthreads();
    if (threadIdx.x < 3)
    {
        xyz_red_shmem[3+threadIdx.x] = xyz_shmem[3+threadIdx.x]/r_shmem[1];}
    __syncthreads();
    if (threadIdx.x < 3)
    {
        xyz_shmem[6+threadIdx.x] = xyz_shmem[3+threadIdx.x] - xyz_shmem[threadIdx.x];
      }
      __syncthreads();
      if (threadIdx.x == 0)
      {
      r_shmem[2] = sqrt( xyz_shmem[6]* xyz_shmem[6]+ xyz_shmem[7]* xyz_shmem[7]+ xyz_shmem[8]* xyz_shmem[8]);
      }
      __syncthreads();
      if (threadIdx.x == 0)
      {
      xyz_red_shmem[6] = xyz_shmem[6]/r_shmem[2];
      xyz_red_shmem[7] = xyz_shmem[7]/r_shmem[2];
      xyz_red_shmem[8] = xyz_shmem[8]/r_shmem[2];
      }
      __syncthreads();
    /*
      if (tid%64 == 0) printf (
	  "launching kernel, loop variables are: i %d j: %d j2: %d and n_neigh[i] is: %d i3 is %d, j3 is %d , k3 is %d with 1 r: %g, xyz: %g,%g,%g, xyzred: %g,%g,%g with 2 r: %g, xyz: %g,%g,%g, xyzred: %g,%g,%g with 3 r: %g, xyz: %g,%g,%g, xyzred: %g,%g,%g\n"
	  ,i,j,j2,n_neigh[i],i3,j3,k3
	  ,r_shmem[0],xyz_shmem[0],xyz_shmem[1],xyz_shmem[2], xyz_red_shmem[0],xyz_red_shmem[1],xyz_red_shmem[2]
	  ,r_shmem[1],xyz_shmem[3],xyz_shmem[4],xyz_shmem[5], xyz_red_shmem[3],xyz_red_shmem[4],xyz_red_shmem[5]
	  ,r_shmem[2],xyz_shmem[6],xyz_shmem[7],xyz_shmem[8], xyz_red_shmem[6],xyz_red_shmem[7],xyz_red_shmem[8]
	  );*/
      //since is in the kernel i want to make the switch an "if constexpr" wall, and bring the switch outside the gpu
      counter_launch++;
      if constexpr (type == exp_k){
	if constexpr (do_forces == true){
          eval_energies_d<true,exp_k>(rcut, buffer, n_sparse, xyz_shmem, xyz_red_shmem, r_shmem, sigma, qs, pref_d, energies, forces, virial_loc, fcut, dfcut, i, i3, j3, k3, n_sites0) ;
	}else{
          eval_energies_d<false,exp_k>(rcut, buffer, n_sparse, xyz_shmem, xyz_red_shmem, r_shmem, sigma, qs, pref_d, energies, forces, virial_loc, fcut, dfcut, i, i3, j3, k3, n_sites0) ;
	}
      }


      if constexpr (type == pol_k){
	if constexpr (do_forces == true){
          eval_energies_d<true,pol_k>(rcut, buffer, n_sparse, xyz_shmem, xyz_red_shmem, r_shmem, sigma, qs, pref_d, energies, forces, virial_loc, fcut, dfcut, i, i3, j3, k3, n_sites0) ;
	}else{
          eval_energies_d<false,pol_k>(rcut, buffer, n_sparse, xyz_shmem, xyz_red_shmem, r_shmem, sigma, qs, pref_d, energies, forces,virial_loc, fcut, dfcut, i, i3, j3, k3, n_sites0) ;
	}
      }
    }//j2 loop
//   printf("tid %d from block: %d, i have skipped first loop %d times, second %d times and launched kernel %d times\n",tid,i,
//   counter_first_cont, counter_second_cont, counter_launch );
  }//j loop

  if constexpr (do_forces == true)
  {
    if (tid == 0)
    {
      #pragma unroll
      for(auto tmp=0; tmp<9; ++tmp)
      {
        atomicAdd(&virial[tmp], virial_loc[tmp]);
      }
    }
  }//doforces

}


#if 0
    __global__ void sanity_check(const double* rjs, const int* neighbors_list, const int n_atom_pairs){


      printf("rjs is:\n");
      for (int i=0; i < n_atom_pairs; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",rjs[i]);
      }
      printf("neighbors list is:\n");
      for (int i=0; i < n_atom_pairs; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%d ",neighbors_list[i]);
      }
    }
#endif

//    __global__ void sanity_check(const double* rjs, const int* neighbors_list, const int n_atom_pairs){
__global__ void sanity_check(
    const int n_sparse, 
    const int n_sites, 
    const int n_atom_pairs, 
    const int n_sites0, 
    const int sp0,
    const int sp1,
    const int sp2,
    const double* alpha,
    const double delta,
    const double* cutoff,
    //host arrays
    const double* rjs,
    const double* xyz,
    const int* n_neigh,
    const int* species,
    const int* neighbors_list,
    const int* neighbor_species ,
    const int* kappas_array,
    const double rcut,
    const double buffer,
    //const double e0   ,
    //device array      
    const double* sigma,
    const double* qs,
    const double* pref_d,
    //output arrays
    double* energies,
    double* forces,
    double* virial
){

      printf("alpha is:\n");
      for (int i=0; i < n_sparse; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",alpha[i]);
      }
      printf("\ncutoff is:\n");
      for (int i=0; i < n_sparse; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",cutoff[i]);
      }
      printf("\nrjs is:\n");
      for (int i=0; i < n_atom_pairs; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",rjs[i]);
      }
      printf("\nxyz is:\n");
      for (int i=0; i <3*n_atom_pairs; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",xyz[i]);
      }
      printf("\nn_neigh is:\n");
      for (int i=0; i < n_sites; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%d ",n_neigh[i]);
      }
      printf("\nspecies is:\n");
      for (int i=0; i < n_sites; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%d ",species[i]);
      }
      printf("\nneighbors list is:\n");
      for (int i=0; i < n_atom_pairs; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%d ",neighbors_list[i]);
      }
      printf("\nneighbors species is:\n");
      for (int i=0; i < n_atom_pairs; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%d ",neighbor_species[i]);
      }
      printf("\nkappas is:\n");
      for (int i=0; i < n_sites; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%d ",kappas_array[i]);
      }
      printf("\nsigma is:\n");
      for (int i=0; i < 3; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",sigma[i]);
      }
      printf("\nqs is:\n");
      for (int i=0; i < 3*n_sparse; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",qs[i]);
      }
      printf("\npref_d is:\n");
      for (int i=0; i < n_sparse; ++i)
      {
	if (i%10 == 0) printf("\n");
	printf("%.16f ",pref_d[i]);
      }

    }

void setup_3bresult_arrays(void** energies,void** forces,void** virial,const hipStream_t s, const bool do_forces, const int n_sites, const int n_sites0)
{
  printf("inside setup: nsites is %d, nsites0 is %d\n",n_sites,n_sites0);
  if( do_forces )
  {
    //3 seems to be hardcoded in allocate(forces) in original turbogap code.
    hip_check_error( hipMallocAsync(forces, n_sites0 *3*sizeof(double), s));
    hip_check_error( hipMemsetAsync(forces, 0, n_sites0 *3*sizeof(double), s));
    hip_check_error( hipMallocAsync(virial, 9*sizeof(double), s));
    hip_check_error( hipMemsetAsync(virial, 0, 9*sizeof(double), s));      
  }
  hip_check_error(hipMallocAsync(energies, n_sites* n_sites *sizeof(double), s));
  hip_check_error(hipMemsetAsync(energies, 0, n_sites*n_sites*sizeof(double), s));      
  printf("allocated energies at address %p \n",energies);
  {
    //using scope to override "global" threads, since this values are only used for initialization.
//    auto threads = n_sites  < 1024 ? n_sites : 1024;
//    auto grids = n_sites < 1024 ? 1 : n_sites/1024 +1;
//    init_energy<<<grids,threads,0,s>>>(energies, e0, n_sites*n_sites);

    // double * energies_h;
    // hipHostMalloc((void**)&energies_h, n_sites*sizeof(double));
    // hipMemcpyAsync(energies_h, energies, n_sites*sizeof(double),hipMemcpyDeviceToHost ,s);
    // hipDeviceSynchronize();
  }
}

void cleanup_3bresult_arrays(void** energy_d, void** forces_d, void** virials_d, void** energy_h, void** forces_h, void** virials_h,const hipStream_t s, const bool do_forces, const int n_sites, const int n_sites0)
{

  if (do_forces)
  {
    hip_check_error( hipMemcpyAsync(*forces_h, *forces_d, n_sites0 *3*sizeof(double),hipMemcpyDeviceToHost ,s));
    hip_check_error( hipFreeAsync(forces_d, s));
    hip_check_error( hipMemcpyAsync(*virials_h, *virials_d,9*sizeof(double),hipMemcpyDeviceToHost ,s));
    hip_check_error( hipFreeAsync(virials_d, s));
  }
  hip_check_error( hipMemcpyAsync(*energy_h, *energy_d,n_sites0*sizeof(double),hipMemcpyDeviceToHost ,s));
  hip_check_error( hipFreeAsync(energy_d, s));
}

__global__ void print_kappas (int* kappas_array_d, const int n_sites)
{
  auto tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid< n_sites)
    printf("tid %d has kappa %d\n",tid, kappas_array_d[tid]);
}


void create_kappas (int* kappas_array_d, int* n_neigh_host,const hipStream_t s, const int n_sites )
{
  int* kappas_array_h; 
  hip_check_error(hipHostMalloc((void**) &kappas_array_h, n_sites*sizeof(int)));
   
  setup_kappas(kappas_array_h,n_sites,n_neigh_host);

  hip_check_error(hipMemcpyAsync(kappas_array_d, kappas_array_h, n_sites*sizeof(int), hipMemcpyHostToDevice, s));
  hip_check_error(hipStreamSynchronize(s));
  hip_check_error(hipHostFree(kappas_array_h));
//    auto threads = n_sites  < 1024 ? n_sites : 1024;
//    auto grids = n_sites < 1024 ? 1 : n_sites/1024 +1;
//    print_kappas<<<grids,threads,0,s>>>(kappas_array_d, n_sites);

//  hip_check_error(hipStreamSynchronize(s));
}



void gpu_3b_cc_2nd_try(

const int n_sparse    , 
const int n_sites     , 
const int n_atom_pairs, 
const int n_sites0    , 
const int sp0         ,
const int sp1         ,
const int sp2         ,
const double* alpha   ,
const double delta    ,
const double e0   ,
const double* cutoff  ,
const hipStream_t s   ,
//host arrays
const double* rjs           ,
const double* xyz           ,
const int* n_neigh          ,
const int* species          ,
const int* neighbors_list   ,
const int* neighbor_species ,
const bool do_forces,
const double rcut     ,
const double buffer   ,
//device array      
const double* sigma   ,
const double* qs   ,
const kern_type kernel_type,
double* energies,
double* forces, 
double* virial,
const int i_beg,
const int i_end,
const int* kappas_array_d
)
{
  //const auto pi = std::numbers::pi_v<double>; 
  const auto pi =3.14159265358979323846264338327950288419716939; // std::numbers::pi_v<double>; 
  double* pref_d;
  double* pref_h; 
  int device_id;
  hip_check_error(hipGetDevice (&device_id));
  hipDeviceProp_t props;
  hip_check_error(hipGetDeviceProperties(&props, device_id));
  auto threads = props.warpSize; 
//  std::cout<<"running with "<<threads<<" threads"<<std::endl<<std::flush;

//  hip_check_error(hipHostMalloc((void**)&pref_h, n_sparse*sizeof(double)));
  
  hip_check_error( hipDeviceSynchronize());
  hip_check_error(hipMallocAsync((void**)&pref_d, n_sparse*sizeof(double), s));
  // eval_pref<<<1,threads,0,s>>>(pref_d, alpha, delta, cutoff,n_sparse);
  eval_pref<<<(n_sparse+256-1)/256,256,0,s>>>(pref_d, alpha, delta, cutoff,n_sparse);
//  std::cout << std::setprecision(16) << std::fixed;
//    printf("################### after eval pref: C interface #############\n");
//    fflush(stdout);

  hip_check_error( hipDeviceSynchronize());
//dbg
//hip_check_error(hipMemcpyAsync(pref_h, pref_d, n_sparse*sizeof(double), hipMemcpyDeviceToHost, s));
//hip_check_error(hipDeviceSynchronize());
//for (int i = 0; i < n_sparse; ++i)
//  printf("%.17e ",pref_h[i]);
//printf("\n");
//end dbg

 {
    //using scope to override "global" threads, since this values are only used for initialization.
//    printf("################### C interface before init energy #############\n");
//    fflush(stdout);
    auto threads = n_sites  < 1024 ? n_sites : 1024;
    auto grids = n_sites < 1024 ? 1 : n_sites/1024 +1;
//    printf(" grids is: %d, threads is %d, n_sites is %d, ptr to energies is %p \n",grids, threads, n_sites, energies);
//    fflush(stdout);
    printf(" \n Energy cont %lf \n", e0);
    fflush(stdout);
    if (e0 != 0)
    {
      add_init_energy<<<grids,threads,0,s>>>(energies, e0, i_end, i_beg);
    }
 }
  hip_check_error( hipDeviceSynchronize());
//    printf("################### C interface after init energy #############\n");
//    fflush(stdout);


//  int* kappas_array_d;
//  int* kappas_array_h; 
//  hip_check_error(hipHostMalloc((void**) &kappas_array_h, n_sites*sizeof(int)));
//  hip_check_error(hipMallocAsync((void**)&kappas_array_d, n_sites*sizeof(int), s));
//   
//  setup_kappas(kappas_array_h,n_sites,n_neigh_host);
//
//  hip_check_error(hipMemcpyAsync(kappas_array_d, kappas_array_h, n_sites*sizeof(int), hipMemcpyHostToDevice, s));
//  create_kappas (kappas_array_d, n_neigh_host, n_sites ); 
  auto grids=i_end - i_beg+1;//n_sites;
//    printf("################### C interface: grids are %d, threads are %d  \n parameters are: nsparse %d, nsites %d, natompair %d, nsite0 %d \n",grids,threads,n_sparse,n_sites,n_atom_pairs,n_sites0);
//    printf("sp0 %d, sp1 %d, sp2 %d, alpha[0] %f, delta %f, cutoff[0] %f, rjs[0] %f, xyz[0] %f, n_neigh[0] %d,",sp0, sp1, sp2, alpha[0], delta, cutoff[0], rjs[0], xyz[0], n_neigh[0]);
//    printf("species[0] %d, neighbors_list[0] %d, neighbor_species[0] %d, kappas_array_d[0] %d, rcut %f, buffer %f, ",species[0], neighbors_list[0], neighbor_species[0], kappas_array_d[0], rcut, buffer);
//    printf("sigma[0] %f, qs[0] %f, pref_d[0] %f, energies[0] %f, forces[0] %f, virial[0] %f, i_beg %d \n",sigma[0], qs[0], pref_d[0], energies[0], forces[0], virial[0], i_beg);
//    fflush(stdout);
  
//   hip_check_error( hipStreamSynchronize(s));
//   hip_check_error( hipDeviceSynchronize());
//   auto err = 	hipPeekAtLastError();
//   std::cout<< "reported error is " <<	hipGetErrorString(err) <<std::flush << std::endl;
//   printf("launching kernel\n");
  printf("starting kernel with %d\n", threads);
  fflush(stdout);

  auto clock_start = std::chrono::system_clock::now();
  switch(kernel_type)
  {
    case exp_k:
    {
      do_forces ?
        kernel_2nd_try<true,exp_k> <<<grids,threads,0,s>>>(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, cutoff, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, kappas_array_d, rcut, buffer, /*e0,*/ sigma, qs, pref_d, energies, forces, virial, i_beg) :
        kernel_2nd_try<false,exp_k><<<grids,threads,0,s>>>(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, cutoff, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, kappas_array_d, rcut, buffer, /*e0,*/ sigma, qs, pref_d, energies, forces, virial, i_beg) ;
      break;
    }
    case pol_k:
    {
      do_forces ?
        kernel_2nd_try<true,pol_k> <<<grids,threads,0,s>>>(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, cutoff, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, kappas_array_d, rcut, buffer, /*e0,*/ sigma, qs, pref_d, energies, forces, virial, i_beg) :
        kernel_2nd_try<false,pol_k><<<grids,threads,0,s>>>(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, cutoff, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, kappas_array_d, rcut, buffer, /*e0,*/ sigma, qs, pref_d, energies, forces, virial, i_beg) ;
      break;
    }
  }
  // sanity_check<<<1,1>>>(rjs,neighbors_list,n_atom_pairs);
  // sanity_check<<<1,1>>>(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, cutoff, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, kappas_array_d, rcut, buffer, /*e0,*/ sigma, qs, pref_d, energies, forces, virial);

  hip_check_error( hipDeviceSynchronize());
//  printf("outta loops\n");
//  fflush(stdout);
  auto clock_now = std::chrono::system_clock::now();
  float currentTime = float(std::chrono::duration_cast <std::chrono::microseconds> (clock_now - clock_start).count());
//  std::cout << "Elapsed Time: " << currentTime /1000000 << " S \n";
//  fflush(stdout);
//  if (do_forces)
//  {
//    double * forces_h;
//    hip_check_error( hipHostMalloc((void**)&forces_h, n_sites0 *3*sizeof(double)));
//    hip_check_error( hipMemcpyAsync(forces_h, forces, n_sites0 *3*sizeof(double),hipMemcpyDeviceToHost ,s));
//    double * virials_h;
//    hip_check_error( hipHostMalloc((void**)&virials_h, 9*sizeof(double)));
//    hip_check_error( hipMemcpyAsync(virials_h, virial,9*sizeof(double),hipMemcpyDeviceToHost ,s));
//  }
//  double * energies_h;
//  hip_check_error( hipHostMalloc((void**)&energies_h, n_sites*sizeof(double)));
//  hip_check_error( hipMemcpyAsync(energies_h, energies,n_sites*sizeof(double),hipMemcpyDeviceToHost ,s));
//  hip_check_error( hipDeviceSynchronize());
////  printf("forces are: ");
////  for (int i=0; i<n_sites0*3; ++i)
////  {
////    if ((i%10) == 0) printf("\n");
////    std::cout<< Double(forces_h[i])<<",";
////  }
//  printf("\n energies are: ");
//  for (int i=0; i<n_sites; ++i)
//  {
//    if ((i%10) == 0) printf("\n");
//    std::cout<< Double(energies_h[i])<<",";
//  }
//  printf("\n");

//  printf("\n virials are: \n");
//  for (int i=0; i<9; ++i)
//  {
//    std::cout<< Double(virials_h[i])<<",";
//  }
//  printf("\n");
  hip_check_error( hipFreeAsync(pref_d,s) );

}
extern "C"
{
  void gpu_3b(
const int n_sparse    , 
const int n_sites     , 
const int n_atom_pairs, 
const int n_sites0    , 
const int sp0         ,
const int sp1         ,
const int sp2         ,
const double* alpha   ,
const double delta    ,
const double e0       ,
const double* cutoff  ,
const hipStream_t* s  ,
//host arrays
const double* rjs           ,
const double* xyz           ,
const int* n_neigh          ,
const int* species          ,
const int* neighbors_list   ,
const int* neighbor_species ,
const bool do_forces,
const double rcut,
const double buffer,
const double* sigma,
const double* qs,
const char* str,
const int i_beg,
const int i_end, //not sure is really needed though. need to check with nsites.
double* energy,
double* forces,
double* virials,
const int* kappas_array_d
)
  {

    kern_type k;
    if(!strcmp(str,"pol"))
    {
      k = pol_k;
//      printf("calling pol_k\n");
    }
    else if(!strcmp(str,"exp"))
    {
      k = exp_k;
//      printf("calling exp_k\n");
    }
    else 
    {
      printf("kernel type is not recognized!\n");
    }
    gpu_3b_cc_2nd_try(n_sparse, n_sites, n_atom_pairs, n_sites0, sp0, sp1, sp2, alpha, delta, e0, cutoff, *s, rjs, xyz, n_neigh, species, neighbors_list, neighbor_species, do_forces, rcut, buffer, sigma, qs, k, energy, forces, virials, i_beg, i_end,kappas_array_d);
  }


void create_kappas_cwrap (int* kappas_array_d, int* n_neigh_host,
    const hipStream_t* s  ,
    const int n_sites )
{
  create_kappas(kappas_array_d, n_neigh_host, *s, n_sites);
}

//allocate and 0 init device arrays
void setup_3bresult_arrays_cwrap(void** energy_d, void** forces_d, void** virials_d, 
    const hipStream_t* s  ,
    const bool do_forces, const int n_sites, const int n_sites0)
{
//    printf("################### begin C interface: setup 3b array #############\n");
  setup_3bresult_arrays(energy_d, forces_d, virials_d, *s, do_forces, n_sites,n_sites0);
  hip_check_error( hipDeviceSynchronize());
}

//copy back device to host, deallocate device
void cleanup_3bresult_arrays_cwrap(void** energy_d, void** forces_d, void** virials_d,
    void** energy_h, void** forces_h, void** virials_h, 
    const hipStream_t* s  ,
    const bool do_forces, const int n_sites, const int n_sites0)
{
  hip_check_error( hipDeviceSynchronize());
//    printf("################### begin C interface: cleanup 3b array #############\n");
  cleanup_3bresult_arrays(energy_d, forces_d, virials_d, energy_h, forces_h, virials_h, *s, do_forces, n_sites,n_sites0);
}
}
