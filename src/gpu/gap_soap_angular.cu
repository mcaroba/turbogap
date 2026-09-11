// The soap_turbo angular expansion: the associated Legendre arrays, the
// exp(i m phi) coefficients and their derivatives, and the cnk contraction of
// the radial and angular parts.
#include "gpu_common.h"
#include "gap_gpu.h"

#define tpb 64
#define tpbcnk 64 // this is because k_max is 45???

__global__ void cuda_get_cnk_one(hipDoubleComplex* cnk_d, double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d,
                                 int* n_neigh_d, int* k2_start_d, int n_atom_pairs, int n_sites, int k_max, int n_max, int l_max) {
  int i_site = threadIdx.x + blockIdx.x * blockDim.x;
  double pi = 4.0 * acos(-1.0);
  if (i_site < n_sites) {
    int k2 = k2_start_d[i_site];
    int my_n_neigh = n_neigh_d[i_site];
    for (int j = 1; j <= my_n_neigh; j++) {
      k2++;
      for (int n = 1; n <= n_max; n++) {
        double loc_rad_exp_coeff = radial_exp_coeff_d[n - 1 + n_max * (k2 - 1)];
        for (int l = 0; l <= l_max; l++) {
          for (int m = 0; m <= l; m++) {
            int k = 1 + l * (l + 1) / 2 + m;                                                  //k=1+
            hipDoubleComplex loc_cnk = cnk_d[i_site + n_sites * ((k - 1) + (n - 1) * k_max)]; //cnk_d[k-1+k_max*(n-1+n_max*i_site)];
            hipDoubleComplex loc_ang_exp_coeff =
                angular_exp_coeff_d[k2 - 1 + n_atom_pairs * (k - 1)]; // angular_exp_coeff_d[k-1+k_max*(k2-1)];
            loc_cnk.x += pi * loc_rad_exp_coeff * loc_ang_exp_coeff.x;
            loc_cnk.y += pi * loc_rad_exp_coeff * loc_ang_exp_coeff.y;
            cnk_d[i_site + n_sites * ((k - 1) + (n - 1) * k_max)] = loc_cnk; //cnk_d[k-1+k_max*(n-1+n_max*i_site)]=loc_cnk;
          }
        }
      }
    }
  }
}

__global__ void cuda_get_cnk_one_new_new(hipDoubleComplex* cnk_d, double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d,
                                         int* n_neigh_d, int* k2_start_d, int n_atom_pairs, int n_sites, int k_max, int n_max,
                                         int l_max) {
  int i_site = threadIdx.x + blockIdx.x * blockDim.x;
  int n = blockIdx.y + 1;
  int k = blockIdx.z + 1;
  double pi = 4.0 * acos(-1.0);
  if (i_site < n_sites) {
    int k2 = k2_start_d[i_site];
    int my_n_neigh = n_neigh_d[i_site];
    hipDoubleComplex loc_cnk;
    loc_cnk.x = 0.0;
    loc_cnk.y = 0.0;
    /*if(k<=k_max){
	cnk_d[k-1+k_max*(n-1+n_max*i_site)]; // coalesced???
	}*/
    for (int j = 1; j <= my_n_neigh; j++) {
      k2++;
      double loc_rad_exp_coeff = radial_exp_coeff_d[n - 1 + n_max * (k2 - 1)]; //coalesced ???
      //int k=1+l*(l+1)/2+m;
      hipDoubleComplex loc_ang_exp_coeff =
          angular_exp_coeff_d[k2 - 1 + n_atom_pairs * (k - 1)]; // angular_exp_coeff_d[k-1+k_max*(k2-1)]; //coalesced ??
      loc_cnk.x += pi * loc_rad_exp_coeff * loc_ang_exp_coeff.x;
      loc_cnk.y += pi * loc_rad_exp_coeff * loc_ang_exp_coeff.y;
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
    if (k <= k_max) {
      cnk_d[i_site + n_sites * ((k - 1) + (n - 1) * k_max)] = loc_cnk; //cnk_d[k-1+k_max*(n-1+n_max*i_site)]=loc_cnk;
    }
  }
}

__global__ void cuda_get_cnk_two(hipDoubleComplex* cnk_d, double* atom_sigma_r, double* atom_sigma_t, double* rcut_hard,
                                 double* central_weight, int* species, int* i_beg, int* i_end, int radial_enhancement,
                                 int* species_multiplicity, double* W, double* S, int n_sites, int k_max, int n_max,
                                 int size_species_1) {
  int i_site = threadIdx.x + blockIdx.x * blockDim.x; //if (i_site >= n_sites) return;
  double pi = acos(-1.0);

  if (i_site < n_sites) {
    for (int k = 1; k <= species_multiplicity[i_site]; k++) {
      int j = species[i_site * size_species_1 + k - 1] - 1;
      double amplitude;
      double atom_sigma_r_j = atom_sigma_r[j];
      double atom_sigma_t_j = atom_sigma_t[j];
      double rcut_hard_j = rcut_hard[j];
      double central_weight_j = central_weight[j];
      if (radial_enhancement == 1) {
        amplitude = sqrt(2.0 / pi) * atom_sigma_r_j / rcut_hard_j;
      } else if (radial_enhancement == 2) {
        amplitude = (atom_sigma_r_j * atom_sigma_r_j) / (rcut_hard_j * rcut_hard_j);
      } else {
        amplitude = 1.0;
      }

      int i_beg_j = i_beg[j];
      int i_end_j = i_end[j];

      for (int n = i_beg_j; n <= i_end_j; n++) {
        double mmul_WS = 0.0;
        for (int d = i_beg_j; d <= i_end_j; d++) {
          mmul_WS += W[n - 1 + (d - 1) * n_max] * S[d - 1 + (i_end_j - 1) * n_max];
          if (isnan(W[n - 1 + (d - 1) * n_max])) {
            printf("W is nan %lf\n", W[n - 1 + (d - 1) * n_max]);
          }
          if (isnan(S[d - 1 + (i_end_j - 1) * n_max])) {
            printf("S is nan %lf\n", S[d - 1 + (i_end_j - 1) * n_max]);
          }
          /*if(n>n_max || d>n_max){
	      printf("%d %d %d %d\n", i_site, l, d, n_max);
	      }*/
        }
        hipDoubleComplex l_cnk = cnk_d[i_site + n_sites * ((n - 1) * k_max)]; //cnk_d[k_max*(n-1+n_max*i_site)];
        l_cnk.x += amplitude * central_weight_j * sqrt(4.0 * pi) * sqrt(sqrt(pi)) * sqrt(atom_sigma_r_j / 2.0) *
                   (rcut_hard_j * rcut_hard_j * rcut_hard_j) / (atom_sigma_t_j * atom_sigma_t_j) / atom_sigma_r_j * mmul_WS;
        cnk_d[i_site + n_sites * ((n - 1) * k_max)] = l_cnk; //cnk_d[k_max*(n-1+n_max*i_site)]=l_cnk;
      }
    }
  }
}

extern "C" void gpu_get_cnk(double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d, hipDoubleComplex* cnk_d,
                            int* n_neigh_d, int* k2_start_d, int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max,
                            int l_max, int bintybint, double* atom_sigma_r_d, double* atom_sigma_t_d, double* rcut_hard_d,
                            double* central_weight_d, int* species_d, int* i_beg_d, int* i_end_d, int radial_enhancement,
                            int* species_multiplicity_d, double* W_d, double* S_d, int size_species_1, hipStream_t* stream) {
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

  dim3 nnth = dim3(tpbcnk, 1, 1); // each block does the inner loops over l and m in total k_max per block

  /*hipEventRecord(start);

      for(int lll=1;lll<=1000;lll++){*/
  cuda_get_cnk_one_new_new<<<dim3((n_sites + tpbcnk - 1) / tpbcnk, n_max, k_max), nnth, 0, stream[0]>>>(
      cnk_d, radial_exp_coeff_d, angular_exp_coeff_d, n_neigh_d, k2_start_d, n_atom_pairs, n_sites, k_max, n_max, l_max);
  /*}
      hipEventRecord(stop);
      hipEventSynchronize(stop);
      milliseconds = 0.0;
      hipEventElapsedTime(&milliseconds, start, stop);
      printf("\n Time of the second kernel in s %f\n", milliseconds/1000.0);
      exit(0);*/
  if (bintybint == 1000) {
    dim3 nblocks = dim3((n_sites - 1 + tpb) / tpb, 1, 1);
    dim3 nthreads = dim3(tpb, 1, 1);
    cuda_get_cnk_two<<<nblocks, nthreads, 0, stream[0]>>>(cnk_d, atom_sigma_r_d, atom_sigma_t_d, rcut_hard_d, central_weight_d,
                                                          species_d, i_beg_d, i_end_d, radial_enhancement, species_multiplicity_d,
                                                          W_d, S_d, n_sites, k_max, n_max, size_species_1);
  }
}

__global__ void cuda_get_plm_arrays_one(double* plm_array_global_d, int kmax, int lmax, double* thetas_d, int n_atom_pairs) {
  int k_ij = threadIdx.x + blockIdx.x * blockDim.x;
  if (k_ij < n_atom_pairs) {
    double x = cos(thetas_d[k_ij]);
    // compute the first 6 polynomials to initialize the recursion series
    plm_array_global_d[k_ij] = 1.0;                                   //plm_array_global_d[  k_ij*kmax]=1.0;
    plm_array_global_d[k_ij + 1 * n_atom_pairs] = x;                  //plm_array_global_d[1+k_ij*kmax]=x;
    plm_array_global_d[k_ij + 2 * n_atom_pairs] = -sqrt(1.0 - x * x); // plm_array_global_d[2+k_ij*kmax]=-sqrt(1.0-x*x);
    plm_array_global_d[k_ij + 3 * n_atom_pairs] = 1.5 * x * x - 0.5;  //plm_array_global_d[3+k_ij*kmax]=1.5*x*x-0.5;
    plm_array_global_d[k_ij + 4 * n_atom_pairs] =
        -3.0 * x * sqrt(1.0 - x * x);                                //plm_array_global_d[4+k_ij*kmax]=-3.0*x*sqrt(1.0-x*x);
    plm_array_global_d[k_ij + 5 * n_atom_pairs] = 3.0 - 3.0 * x * x; //plm_array_global_d[5+k_ij*kmax]=3.0 -3.0*x*x;

    for (int l = 3; l <= lmax; l++) {
      int k = 0;
      for (int m = 0; m <= l - 2; m++) {
        k = 1 + l * (l + 1) / 2 + m;
        plm_array_global_d[k_ij + (k - 1) * n_atom_pairs] =
            ((2.0 * l - 1.0) * x * plm_array_global_d[k_ij + (k - l - 1) * n_atom_pairs] - //plm_array_global_d[k-1+k_ij*kmax]
             (l - 1.0 + m) * plm_array_global_d[k_ij + (k - 2 * l + 1 - 1) * n_atom_pairs]) /
            (l - m);
      }
      k = k + 1;
      plm_array_global_d[k_ij + (k - 1) * n_atom_pairs] =
          x * (2.0 * l - 1.0) * plm_array_global_d[k_ij + (k - l - 1) * n_atom_pairs]; //plm_array_global_d[k-1+k_ij*kmax]
      k = k + 1;
      plm_array_global_d[k_ij + (k - 1) * n_atom_pairs] =
          -(2.0 * l - 1.0) * sqrt(1.0 - x * x) *
          plm_array_global_d[k_ij + (k - l - 1 - 1) * n_atom_pairs]; //plm_array_global_d[k-1+k_ij*kmax]
    }
  }
}

extern "C" void gpu_get_plm_array_global(double* plm_array_global_d, int n_atom_pairs, int kmax, int lmax, double* thetas_d,
                                         hipStream_t* stream) {
  dim3 nblocks = dim3((n_atom_pairs - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);
  cuda_get_plm_arrays_one<<<nblocks, nthreads, 0, stream[0]>>>(plm_array_global_d, kmax, lmax, thetas_d, n_atom_pairs);
}


__global__ void cuda_get_exp_coeff_one(hipDoubleComplex* eimphi_global_d, double* rjs_d, double* phis_d, bool* mask_d,
                                       double* atom_sigma_in_d, double* atom_sigma_scaling_d, double rcut, int n_atom_pairs,
                                       int n_species, int lmax, int kmax, double* prefl_global_d, hipDoubleComplex* prefm_global_d,
                                       double* preflm_d, double* plm_array_global_d,
                                       hipDoubleComplex* exp_coeff_d) //, double *fact_array_d)
{
  int k_ij = threadIdx.x + blockIdx.x * blockDim.x;
  double xcut = 1.0e-7;
  if (k_ij < n_atom_pairs) {
    double rj = rjs_d[k_ij];
    double phi = phis_d[k_ij];
    if (rj < rcut) {
      int i_sp = 1;
      for (i_sp = 1; i_sp <= n_species; i_sp++)
        if (mask_d[k_ij + (i_sp - 1) * n_atom_pairs])
          break;
      double atom_sigma = atom_sigma_in_d[i_sp - 1] + atom_sigma_scaling_d[i_sp - 1] * rj;
      double scaling = atom_sigma_scaling_d[i_sp - 1];
      double rjbysigma = rj / atom_sigma;
      double amplitude = (rcut * rcut) / (atom_sigma * atom_sigma);
      double x = rjbysigma;
      double x2 = x * x;
      double x4 = x2 * x2;
      double flm2, flm1, fl;
      double coeff1 = 2.0 * rj / atom_sigma * atom_sigma;
      double coeff2 = 1.0 - scaling * rj / atom_sigma;

      if (x > 0) {
        flm2 = fabs((1.0 - exp(-2.0 * x2)) / 2.0 / x2);
        flm1 = fabs((x2 - 1.0 + exp(-2.0 * x2) * (x2 + 1.0)) / 2.0 / x4);
      } else {
        flm2 = 1.0;
        flm1 = 0.0;
      }
      //  Complex exponential using Euler's formula and Chebyshev recursion
      double cosm2 = cos(phi);
      double cosphi2 = 2.0 * cosm2;
      double sinm2 = -sin(phi);
      double cosm1 = 1.0;
      double sinm1 = 0.0;
      hipDoubleComplex loc_prefm;
      loc_prefm.x = 1.0;
      loc_prefm.y = 0.0; // (1.d0, 0.d0)*
      prefm_global_d[k_ij] = loc_prefm;
      double ilexp = -500000; // no need for this, just makig sure the next lines were working

      double fact = 1.0;
      int k = 0;
      for (int l = 0; l <= lmax; l++) {
        if (l > 0)
          fact = fact * (2.0 * l + 1.0);
        if (l == 0) {
          if (x < xcut)
            ilexp = 1.0 - x2;
          else
            ilexp = flm2;
        } else if (l == 1) {
          if (x2 / 1000.0 < xcut)
            ilexp = (x2 - x4) / fact; //fact_array_d[l-1];
          else
            ilexp = flm1;
        } else {
          if (pow(x2, l) / fact * l < xcut)
            fl = pow(x2, l) / fact;
          else
            fl = fabs(flm2 - (2.0 * l - 1.0) / x2 * flm1);
          flm2 = flm1;
          flm1 = fl;
          ilexp = fl;
        }
        if (l > 0) {
          double cos0 = cosphi2 * cosm1 - cosm2;
          double sin0 = cosphi2 * sinm1 - sinm2;
          cosm2 = cosm1;
          sinm2 = sinm1;
          cosm1 = cos0;
          sinm1 = sin0;
          loc_prefm.x = cos0;
          loc_prefm.y = -sin0;
          prefm_global_d[k_ij + l * n_atom_pairs] = loc_prefm;
        }
        prefl_global_d[k_ij + l * n_atom_pairs] = ilexp;
        for (int m = 0; m <= l; m++) {
          hipDoubleComplex loc_emphi;
          hipDoubleComplex tmp_prefm = prefm_global_d[k_ij + m * n_atom_pairs];
          loc_emphi.x = ilexp * tmp_prefm.x;
          loc_emphi.y = ilexp * tmp_prefm.y;
          eimphi_global_d[k_ij + k * n_atom_pairs] = loc_emphi;
          hipDoubleComplex loc_exp_coeff;
          loc_exp_coeff.x = amplitude * preflm_d[k] * plm_array_global_d[k_ij + k * n_atom_pairs] * loc_emphi.x;
          loc_exp_coeff.y = amplitude * preflm_d[k] * plm_array_global_d[k_ij + k * n_atom_pairs] * loc_emphi.y;
          exp_coeff_d[k_ij + k * n_atom_pairs] = loc_exp_coeff;
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


__global__ void cuda_get_exp_coeff_der_one(hipDoubleComplex* eimphi_global_d, double* rjs_d, double* phis_d, bool* mask_d,
                                           double* atom_sigma_in_d, double* atom_sigma_scaling_d, double rcut, int n_atom_pairs,
                                           int n_species, int lmax, int kmax, double* prefl_global_d,
                                           hipDoubleComplex* prefm_global_d, double* preflm_d, double* prefl_global_der_d,
                                           double* plm_array_global_d, double* plm_array_global_der_d,
                                           hipDoubleComplex* exp_coeff_d, hipDoubleComplex* eimphi_rad_der_global_d,
                                           hipDoubleComplex* eimphi_azi_der_global_d, double* plm_array_div_sin,
                                           double* plm_array_der_mul_sin, hipDoubleComplex* exp_coeff_rad_der_d,
                                           hipDoubleComplex* exp_coeff_azi_der_d, hipDoubleComplex* exp_coeff_pol_der_d) {
  int k_ij = threadIdx.x + blockIdx.x * blockDim.x;

  if (k_ij < n_atom_pairs) {
    double rj = rjs_d[k_ij];
    /*double phi=phis_d[k_ij];*/
    if (rj < rcut) {
      int i_sp = 1;
      for (i_sp = 1; i_sp <= n_species; i_sp++) {
        if (mask_d[k_ij + (i_sp - 1) * n_atom_pairs]) {
          break;
        }
      }
      double atom_sigma = atom_sigma_in_d[i_sp - 1] + atom_sigma_scaling_d[i_sp - 1] * rj;
      double scaling = atom_sigma_scaling_d[i_sp - 1];
      double amplitude = (rcut * rcut) / (atom_sigma * atom_sigma);
      /*double rjbysigma=rj/atom_sigma;
	  double x=rjbysigma;
	  double x2=x*x;
	  double x4=x2*x2;/
	  double flm2, flm1,fl; */

      double coeff1 = 2.0 * rj / (atom_sigma * atom_sigma);
      double coeff2 = 1.0 - scaling * rj / atom_sigma;
      //hipDoubleComplex  loc_prefm;
      //double ilexp=-500000;// no need for this, just makig sure the next lines were working


      int k = 0;

      double ilexp_der;

      for (int l = 0; l <= lmax; l++) {
        if (l == 0) {
          ilexp_der = coeff1 * (prefl_global_d[k_ij + n_atom_pairs] - prefl_global_d[k_ij]);
        } else {
          ilexp_der = (-coeff1 - (2.0 * l + 2.0) / rj) * prefl_global_d[k_ij + l * n_atom_pairs] +
                      coeff1 * prefl_global_d[k_ij + (l - 1) * n_atom_pairs];
        }
        if (rj < 1.0e-5) {
          ilexp_der = 0.0;
        }
        ilexp_der *= coeff2;
        prefl_global_der_d[k_ij + l * n_atom_pairs] = ilexp_der;
        for (int m = 0; m <= l; m++) {
          hipDoubleComplex loc_exp_coeff = exp_coeff_d[k_ij + k * n_atom_pairs];
          double loc_preflm = preflm_d[k];
          hipDoubleComplex loc_emphi = eimphi_global_d[k_ij + k * n_atom_pairs];
          hipDoubleComplex tmp_prefm = prefm_global_d[k_ij + m * n_atom_pairs];
          hipDoubleComplex loc_emphi_rad_der;
          loc_emphi_rad_der.x = ilexp_der * tmp_prefm.x;
          loc_emphi_rad_der.y = ilexp_der * tmp_prefm.y;
          eimphi_rad_der_global_d[k_ij + k * n_atom_pairs] = loc_emphi_rad_der;

          hipDoubleComplex loc_emphi_azi_der;
          loc_emphi_azi_der.x = -loc_emphi.y;
          loc_emphi_azi_der.y = loc_emphi.x;
          eimphi_azi_der_global_d[k_ij + k * n_atom_pairs] = loc_emphi_azi_der;

          hipDoubleComplex loc_e_c_rad_der, loc_e_c_azi_der, loc_e_c_pol_der;
          loc_e_c_rad_der.x = amplitude * loc_preflm * plm_array_global_d[k_ij + k * n_atom_pairs] * loc_emphi_rad_der.x -
                              2.0 * atom_sigma_scaling_d[i_sp - 1] / atom_sigma * loc_exp_coeff.x;
          loc_e_c_rad_der.y = amplitude * loc_preflm * plm_array_global_d[k_ij + k * n_atom_pairs] * loc_emphi_rad_der.y -
                              2.0 * atom_sigma_scaling_d[i_sp - 1] / atom_sigma * loc_exp_coeff.y;

          loc_e_c_azi_der.x = amplitude * loc_preflm * plm_array_div_sin[k_ij + k * n_atom_pairs] * loc_emphi_azi_der.x;
          loc_e_c_azi_der.y = amplitude * loc_preflm * plm_array_div_sin[k_ij + k * n_atom_pairs] * loc_emphi_azi_der.y;

          loc_e_c_pol_der.x = amplitude * loc_preflm * plm_array_der_mul_sin[k_ij + k * n_atom_pairs] * loc_emphi.x;
          loc_e_c_pol_der.y = amplitude * loc_preflm * plm_array_der_mul_sin[k_ij + k * n_atom_pairs] * loc_emphi.y;

          // exp_coeff_rad_der_d[k+k_ij*kmax]=loc_e_c_rad_der;
          // exp_coeff_azi_der_d[k+k_ij*kmax]=loc_e_c_azi_der;
          // exp_coeff_pol_der_d[k+k_ij*kmax]=loc_e_c_pol_der;

          exp_coeff_rad_der_d[k_ij + k * n_atom_pairs] = loc_e_c_rad_der;
          exp_coeff_azi_der_d[k_ij + k * n_atom_pairs] = loc_e_c_azi_der;
          exp_coeff_pol_der_d[k_ij + k * n_atom_pairs] = loc_e_c_pol_der;

          k++;
        }
      }
    }
  }
}


__global__ void cuda_get_plm_arrays_der_one(double* plm_array_global_der_d, int kmax, int lmax, double* thetas_d, int n_atom_pairs,
                                            double* plm_array_div_sin, double* plm_array_der_mul_sin) {
  int k_ij = threadIdx.x + blockIdx.x * blockDim.x;
  if (k_ij < n_atom_pairs) {
    //double x=cos(thetas_d[k_ij]);
    double part1, part2;
    for (int l = 0; l <= lmax; l++) {
      for (int m = 0; m <= l; m++) {
        int k = 1 + l * (l + 1) / 2 + m;
        int k_l_mp1 = k + 1;
        int k_l_mm1 = k - 1;
        int k_temp = -5;
        //       If m = 0 then we are asking for P_l^{-1}, which is not defined. We need
        //      to rewrite in terms of P_l^1:
        if (m == 0) {
          // P_0^1=0
          if (l == 0) {
            part1 = 0.0;
            // P_l^{-1} = - (l-1)!/(l+1)! * P_l^1
          } else {
            k_temp = 1 + l * (l + 1) / 2 + 1;
            part1 = -0.5 * plm_array_global_der_d[k_ij + (k_temp - 1) * n_atom_pairs];
          }
        } else {
          part1 = 0.5 * (l + m) * (l - m + 1) * plm_array_global_der_d[k_ij + (k_l_mm1 - 1) * n_atom_pairs];
        }
        if (m == l) {
          part2 = 0.0;
        } else {
          part2 = -0.5 * plm_array_global_der_d[k_ij + (k_l_mp1 - 1) * n_atom_pairs];
        }
        plm_array_der_mul_sin[k_ij + (k - 1) * n_atom_pairs] = part1 + part2;
      }
    }
    for (int l = 0; l <= lmax; l++) {
      for (int m = 0; m <= l; m++) {
        int k = 1 + l * (l + 1) / 2 + m;
        if (m == 0) {
          plm_array_div_sin[k_ij + (k - 1) * n_atom_pairs] = 0.0;
        } else {
          int k_lp1_mp1 = 1 + (l + 1) * (l + 2) / 2 + m + 1;
          int k_lp1_mm1 = 1 + (l + 1) * (l + 2) / 2 + m - 1;
          part1 = 0.5 * (l - m + 1) * (l - m + 2) * plm_array_global_der_d[k_ij + (k_lp1_mm1 - 1) * n_atom_pairs];
          part2 = 0.5 * plm_array_global_der_d[k_ij + (k_lp1_mp1 - 1) * n_atom_pairs];
          plm_array_div_sin[k_ij + (k - 1) * n_atom_pairs] = part1 + part2;
        }
      }
    }
  }
}

extern "C" void gpu_get_exp_coeff_array(hipDoubleComplex* eimphi_global_d, double* rjs_d, double* phis_d, double* thetas_d,
                                        bool* mask_d, double* atom_sigma_in_d, double* atom_sigma_scaling_d, double rcut,
                                        int n_atom_pairs, int n_species, int lmax, int kmax, double* prefl_global_d,
                                        double* plm_array_global_d, double* plm_array_global_der_d, double* prefl_global_der_d,
                                        double* preflm_d, hipDoubleComplex* exp_coeff_d, bool c_do_derivatives,
                                        hipDoubleComplex* eimphi_rad_der_global_d, hipDoubleComplex* eimphi_azi_der_global_d,
                                        double* plm_array_div_sin, double* plm_array_der_mul_sin,
                                        hipDoubleComplex* exp_coeff_rad_der_d, hipDoubleComplex* exp_coeff_azi_der_d,
                                        hipDoubleComplex* exp_coeff_pol_der_d, hipStream_t* stream) {
  dim3 nblocks = dim3((n_atom_pairs - 1 + tpb) / tpb, 1, 1);
  dim3 nthreads = dim3(tpb, 1, 1);

  hipDoubleComplex* prefm_global_d;

  gpuErrchk(hipMallocAsync(&prefm_global_d, n_atom_pairs * (lmax + 1) * sizeof(hipDoubleComplex), stream[0]));

  cuda_get_exp_coeff_one<<<nblocks, nthreads, 0, stream[0]>>>(
      eimphi_global_d, rjs_d, phis_d, mask_d, atom_sigma_in_d, atom_sigma_scaling_d, rcut, n_atom_pairs, n_species, lmax, kmax,
      prefl_global_d, prefm_global_d, preflm_d, plm_array_global_d, exp_coeff_d); //, fact_array_d);


  if (c_do_derivatives) {
    cuda_get_plm_arrays_der_one<<<nblocks, nthreads, 0, stream[0]>>>(plm_array_global_der_d, kmax, lmax, thetas_d, n_atom_pairs,
                                                                     plm_array_div_sin, plm_array_der_mul_sin);


    cuda_get_exp_coeff_der_one<<<nblocks, nthreads, 0, stream[0]>>>(
        eimphi_global_d, rjs_d, phis_d, mask_d, atom_sigma_in_d, atom_sigma_scaling_d, rcut, n_atom_pairs, n_species, lmax, kmax,
        prefl_global_d, prefm_global_d, preflm_d, prefl_global_der_d, plm_array_global_d, plm_array_global_der_d, exp_coeff_d,
        eimphi_rad_der_global_d, eimphi_azi_der_global_d, plm_array_div_sin, plm_array_der_mul_sin, exp_coeff_rad_der_d,
        exp_coeff_azi_der_d, exp_coeff_pol_der_d);
  }

  /*size_t free, total;
      hipMemGetInfo(& free, & total);
      counter++;
      printf("\nFree memory %zu, from %zu in iteration %d\n", free/1024/1024, total/1024/1024, counter); */
  gpuErrchk(hipFreeAsync(prefm_global_d, stream[0]));
}
