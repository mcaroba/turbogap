// cuda_get_soap_p and cuda_get_soap_der_one against the kernels they replaced,
// kept below, on random input: every output must be bit-identical, with and
// without skipped components and with pairs past the cutoff. The der_one
// reference is the old pair of kernels, cnk_*_der written then read. Prints the
// time of each on the same data. gpu_soap_pair_contract is checked, to a
// tolerance, against the chain it replaces done on the host from der_one, and
// cuda_get_soap_index against the walk it replaced.
#include "../../src/gpu/gap_soap_descriptor.cu"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#define cuda_get_soap_der_one reference_soap_der_one
#define cuda_get_soap_p reference_soap_p
static void cuda_get_soap_der_one(double* soap_rad_der_d, double* soap_azi_der_d, double* soap_pol_der_d,
                                  double* multiplicity_array_d, double* trans_soap_rad_der_d, double* trans_soap_azi_der_d,
                                  double* trans_soap_pol_der_d, hipDoubleComplex* cnk_d, hipDoubleComplex* cnk_rad_der_d,
                                  hipDoubleComplex* cnk_azi_der_d, hipDoubleComplex* cnk_pol_der_d, int* k2_i_site_d,
                                  bool* skip_soap_component_d, int n_sites, int n_atom_pairs, int n_soap, int k_max, int n_max,
                                  int l_max, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_soap_soap_der_one", n_atom_pairs, stream,
      TG_LAMBDA(const int k2) {
        int i_site = k2_i_site_d[k2] - 1;
        int counter = 0;
        int counter2 = 0;
        for (int n = 1; n <= n_max; n++) {
          for (int np = n; np <= n_max; np++) {
            for (int l = 0; l <= l_max; l++) {
              if (!skip_soap_component_d
                      [l +
                       (l_max + 1) *
                           (np - 1 +
                            (n - 1) *
                                n_max)]) { //if( skip_soap_component(l, np, n) )cycle // if it happens lots of time, do it in reverse
                counter++;
                double my_soap_rad_der =
                    0; //trans_soap_rad_der_d[k2+(counter-1)*n_atom_pairs]; //soap_rad_der_d[counter-1+k2*n_soap];
                double my_soap_azi_der =
                    0; //trans_soap_azi_der_d[k2+(counter-1)*n_atom_pairs]; //soap_azi_der_d[counter-1+k2*n_soap];
                double my_soap_pol_der =
                    0; //trans_soap_pol_der_d[k2+(counter-1)*n_atom_pairs]; //soap_pol_der_d[counter-1+k2*n_soap];
                for (int m = 0; m <= l; m++) {
                  int k = 1 + l * (l + 1) / 2 + m;
                  counter2++;
                  /* if(threadIdx.x==121 && blockIdx.x==154){
		   printf("\n Pair  %d \n" , k2, i_site);
		   } */
                  hipDoubleComplex tmp_1_cnk_d = cnk_d
                      [i_site +
                       n_sites *
                           (k - 1 +
                            (n - 1) *
                                k_max)]; //trans_cnk_d[i_site+n_sites*(k-1+(n-1)*k_max)];  //cnk_d[k-1+ k_max*(n-1 +i_site*n_max)];
                  hipDoubleComplex tmp_2_cnk_d = cnk_d
                      [i_site +
                       n_sites *
                           (k - 1 +
                            (np - 1) *
                                k_max)]; //trans_cnk_d[i_site+n_sites*(k-1+(np-1)*k_max)]; //cnk_d[k-1+k_max*(np-1+i_site*n_max)];
                  hipDoubleComplex tmp_1_cnk_rad_d = cnk_rad_der_d
                      [k2 +
                       n_atom_pairs *
                           (k - 1 +
                            (n - 1) *
                                k_max)]; //trans_cnk_rad_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; // cnk_rad_der_d[k-1+k_max*(n-1 +k2*n_max)];
                  hipDoubleComplex tmp_2_cnk_rad_d = cnk_rad_der_d
                      [k2 +
                       n_atom_pairs *
                           (k - 1 +
                            (np - 1) *
                                k_max)]; //trans_cnk_rad_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; // cnk_rad_der_d[k-1+k_max*(np-1+k2*n_max)];
                  hipDoubleComplex tmp_1_cnk_azi_d = cnk_azi_der_d
                      [k2 +
                       n_atom_pairs *
                           (k - 1 +
                            (n - 1) *
                                k_max)]; //trans_cnk_azi_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //cnk_azi_der_d[k-1+k_max*(n-1 +k2*n_max)];
                  hipDoubleComplex tmp_2_cnk_azi_d = cnk_azi_der_d
                      [k2 +
                       n_atom_pairs *
                           (k - 1 +
                            (np - 1) *
                                k_max)]; //trans_cnk_azi_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //cnk_azi_der_d[k-1+k_max*(np-1+k2*n_max)];
                  hipDoubleComplex tmp_1_cnk_pol_d = cnk_pol_der_d
                      [k2 +
                       n_atom_pairs *
                           (k - 1 +
                            (n - 1) *
                                k_max)]; //trans_cnk_pol_der_d[k2+n_atom_pairs*(k-1+(n-1)*k_max) ]; //cnk_pol_der_d[k-1+k_max*(n-1 +k2*n_max)];
                  hipDoubleComplex tmp_2_cnk_pol_d = cnk_pol_der_d
                      [k2 +
                       n_atom_pairs *
                           (k - 1 +
                            (np - 1) *
                                k_max)]; //trans_cnk_pol_der_d[k2+n_atom_pairs*(k-1+(np-1)*k_max)]; //cnk_pol_der_d[k-1+k_max*(np-1+k2*n_max)];
                  my_soap_rad_der +=
                      multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_rad_d.x * tmp_2_cnk_d.x + tmp_1_cnk_rad_d.y * tmp_2_cnk_d.y +
                                                            tmp_1_cnk_d.x * tmp_2_cnk_rad_d.x + tmp_1_cnk_d.y * tmp_2_cnk_rad_d.y);
                  my_soap_azi_der +=
                      multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_azi_d.x * tmp_2_cnk_d.x + tmp_1_cnk_azi_d.y * tmp_2_cnk_d.y +
                                                            tmp_1_cnk_d.x * tmp_2_cnk_azi_d.x + tmp_1_cnk_d.y * tmp_2_cnk_azi_d.y);
                  my_soap_pol_der +=
                      multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_pol_d.x * tmp_2_cnk_d.x + tmp_1_cnk_pol_d.y * tmp_2_cnk_d.y +
                                                            tmp_1_cnk_d.x * tmp_2_cnk_pol_d.x + tmp_1_cnk_d.y * tmp_2_cnk_pol_d.y);
                }
                trans_soap_rad_der_d[k2 + (counter - 1) * n_atom_pairs] =
                    my_soap_rad_der; //soap_rad_der_d[counter-1+k2*n_soap]=my_soap_rad_der;
                trans_soap_azi_der_d[k2 + (counter - 1) * n_atom_pairs] =
                    my_soap_azi_der; //soap_azi_der_d[counter-1+k2*n_soap]=my_soap_azi_der;
                trans_soap_pol_der_d[k2 + (counter - 1) * n_atom_pairs] =
                    my_soap_pol_der; //soap_pol_der_d[counter-1+k2*n_soap]=my_soap_pol_der;
              }
            }
          }
        }
      },
      tpb);
}

static void cuda_get_soap_p(double* soap_d, double* sqrt_dot_p_d, double* multiplicity_array_d, hipDoubleComplex* cnk_d,
                            bool* skip_soap_component_d, int n_sites, int n_soap, int n_max, int l_max, hipStream_t* stream) {
  tg_parallel_for(
      "turbogap_soap_soap_p", n_sites, stream,
      TG_LAMBDA(const int i_site) {
        int k_max = 1 + l_max * (l_max + 1) / 2 + l_max;
        double my_sqrt_dot_p = 0.0;
        int counter = 0;
        int counter2 = 0;
        //int ssc_counter=0;
        for (int n = 1; n <= n_max; n++) {
          for (int np = n; np <= n_max; np++) {
            for (int l = 0; l <= l_max; l++) {
              //if(!skip_soap_component_d[ssc_counter]){ //if( skip_soap_component(l, np, n) )cycle
              bool my_skip = skip_soap_component_d[l + (l_max + 1) * (np - 1 + (n - 1) * n_max)];
              if (!(my_skip)) { //if( skip_soap_component(l, np, n) )cycle

                counter++;
                double my_soap = 0.0; //soap_d[counter-1+i_site*n_soap];
                for (int m = 0; m <= l; m++) {
                  int k = 1 + l * (l + 1) / 2 + m; //k = 1 + l*(l+1)/2 + m
                  counter2++;
                  hipDoubleComplex tmp_1_cnk_d =
                      cnk_d[i_site + n_sites * ((k - 1) + (n - 1) * k_max)]; //cnk_d[k-1+k_max*(n-1 +i_site*n_max)];
                  hipDoubleComplex tmp_2_cnk_d =
                      cnk_d[i_site + n_sites * ((k - 1) + (np - 1) * k_max)]; //cnk_d[k-1+k_max*(np-1+i_site*n_max)];
                  my_soap += multiplicity_array_d[counter2 - 1] * (tmp_1_cnk_d.x * tmp_2_cnk_d.x + tmp_1_cnk_d.y * tmp_2_cnk_d.y);
                  /*               if(isnan(my_soap)){
			       printf("\n my_soap is nan %lf %lf %lf %lf %lf!!\n", multiplicity_array_d[counter2-1], tmp_1_cnk_d.x, tmp_1_cnk_d.y,tmp_2_cnk_d.x,tmp_2_cnk_d.y);
			       } */
                  //soap(counter, i) = soap(counter, i) + multiplicity * real(cnk(k, n, i) * conjg(cnk(k, np, i)))
                }
                soap_d[counter - 1 + i_site * n_soap] = my_soap;
                my_sqrt_dot_p += my_soap * my_soap;
              }
            }
          }
        }
        my_sqrt_dot_p = sqrt(my_sqrt_dot_p);
        if (my_sqrt_dot_p < 1.0e-5) {
          my_sqrt_dot_p = 1.0;
        }
        sqrt_dot_p_d[i_site] = my_sqrt_dot_p;
      },
      tpb);
}
#undef cuda_get_soap_der_one
#undef cuda_get_soap_p

__global__ void reference_get_derivatives(double* radial_exp_coeff_d, hipDoubleComplex* angular_exp_coeff_d,
                                          double* radial_exp_coeff_der_d, hipDoubleComplex* angular_exp_coeff_rad_der_d,
                                          hipDoubleComplex* angular_exp_coeff_azi_der_d,
                                          hipDoubleComplex* angular_exp_coeff_pol_der_d, hipDoubleComplex* cnk_rad_der_d,
                                          hipDoubleComplex* cnk_azi_der_d, hipDoubleComplex* cnk_pol_der_d, double* rjs_d,
                                          double rcut_max, int n_atom_pairs, int n_sites, int n_soap, int k_max, int n_max,
                                          int l_max) {
  int k2 = threadIdx.x + blockDim.x * blockIdx.x;
  int n = blockIdx.y + 1;
  int k = blockIdx.z + 1;
  double Pi = 4.0 * acos(-1.0);
  if (k2 < n_atom_pairs) {
    double my_rjs = rjs_d[k2];
    if (my_rjs < rcut_max) {
      double my_radial_exp_c;
      my_radial_exp_c = radial_exp_coeff_d[n - 1 + k2 * n_max];
      double my_radial_exp_c_der;
      my_radial_exp_c_der = radial_exp_coeff_der_d[n - 1 + k2 * n_max];
      hipDoubleComplex my_cnk_rad_der;
      hipDoubleComplex my_cnk_azi_der;
      hipDoubleComplex my_cnk_pol_der;
      hipDoubleComplex my_ang_exp_c;
      my_ang_exp_c = angular_exp_coeff_d[k2 + n_atom_pairs * (k - 1)];
      hipDoubleComplex my_ang_exp_c_rad_der;
      my_ang_exp_c_rad_der = angular_exp_coeff_rad_der_d[k2 + (k - 1) * n_atom_pairs];
      hipDoubleComplex my_ang_exp_c_azi_der;
      my_ang_exp_c_azi_der = angular_exp_coeff_azi_der_d[k2 + (k - 1) * n_atom_pairs];
      hipDoubleComplex my_ang_exp_c_pol_der;
      my_ang_exp_c_pol_der = angular_exp_coeff_pol_der_d[k2 + (k - 1) * n_atom_pairs];

      my_cnk_rad_der.x = Pi * (my_ang_exp_c.x * my_radial_exp_c_der + my_ang_exp_c_rad_der.x * my_radial_exp_c);
      my_cnk_rad_der.y = Pi * (my_ang_exp_c.y * my_radial_exp_c_der + my_ang_exp_c_rad_der.y * my_radial_exp_c);

      my_cnk_azi_der.x = Pi * (my_ang_exp_c_azi_der.x * my_radial_exp_c);
      my_cnk_azi_der.y = Pi * (my_ang_exp_c_azi_der.y * my_radial_exp_c);

      my_cnk_pol_der.x = Pi * (my_ang_exp_c_pol_der.x * my_radial_exp_c);
      my_cnk_pol_der.y = Pi * (my_ang_exp_c_pol_der.y * my_radial_exp_c);

      cnk_rad_der_d[k2 + n_atom_pairs * (k - 1 + (n - 1) * k_max)] = my_cnk_rad_der;
      cnk_azi_der_d[k2 + n_atom_pairs * (k - 1 + (n - 1) * k_max)] = my_cnk_azi_der;
      cnk_pol_der_d[k2 + n_atom_pairs * (k - 1 + (n - 1) * k_max)] = my_cnk_pol_der;
    }
  }
}

static const double rcut_max = 5.0;
static const double CONTRACT_TOL = 1e-11;

struct problem {
  int n_sites;
  int n_atom_pairs;
  int n_max;
  int l_max;
  int k_max;
  int n_soap;
  std::vector<int> k2_i_site;
  std::vector<char> skip;
  std::vector<double> multiplicity;
  std::vector<hipDoubleComplex> cnk;
  std::vector<double> rjs;
  std::vector<double> rad;
  std::vector<double> rad_der;
  std::vector<hipDoubleComplex> ang[4];
};

static problem make_problem(int n_sites, int neigh, int n_max, int l_max, double skip_fraction, unsigned seed) {
  problem p;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  p.n_sites = n_sites;
  p.n_atom_pairs = n_sites * neigh;
  p.n_max = n_max;
  p.l_max = l_max;
  p.k_max = 1 + l_max * (l_max + 1) / 2 + l_max;
  for (int i = 0; i < n_sites; i++)
    for (int j = 0; j < neigh; j++) p.k2_i_site.push_back(i + 1);
  int n_mult = 0;
  p.n_soap = 0;
  p.skip.assign(n_max * n_max * (l_max + 1), 0);
  for (int n = 1; n <= n_max; n++)
    for (int np = n; np <= n_max; np++)
      for (int l = 0; l <= l_max; l++) {
        int s = l + (l_max + 1) * (np - 1 + (n - 1) * n_max);
        p.skip[s] = u01(rng) < skip_fraction;
        if (!p.skip[s]) {
          p.n_soap++;
          n_mult += l + 1;
        }
      }
  for (int i = 0; i < n_mult; i++) p.multiplicity.push_back(u(rng));
  size_t n_cnk = (size_t) n_sites * p.k_max * n_max;
  size_t n_rad = (size_t) p.n_atom_pairs * n_max;
  size_t n_ang = (size_t) p.n_atom_pairs * p.k_max;
  for (size_t i = 0; i < n_cnk; i++) p.cnk.push_back(make_hipDoubleComplex(u(rng), u(rng)));
  for (int k2 = 0; k2 < p.n_atom_pairs; k2++) p.rjs.push_back(rcut_max * 1.1 * u01(rng));
  for (size_t i = 0; i < n_rad; i++) p.rad.push_back(u(rng));
  for (size_t i = 0; i < n_rad; i++) p.rad_der.push_back(u(rng));
  for (int d = 0; d < 4; d++)
    for (size_t i = 0; i < n_ang; i++) p.ang[d].push_back(make_hipDoubleComplex(u(rng), u(rng)));
  return p;
}

template <typename T> static T* upload(const std::vector<T>& v) {
  T* d;
  gpuErrchk(hipMalloc((void**) &d, sizeof(T) * v.size()));
  gpuErrchk(hipMemcpy(d, v.data(), sizeof(T) * v.size(), hipMemcpyHostToDevice));
  return d;
}

struct device_problem {
  int* k2_i_site;
  bool* skip;
  double* multiplicity;
  hipDoubleComplex* cnk;
  double* rjs;
  double* rad;
  double* rad_der;
  hipDoubleComplex* ang[4];
  hipDoubleComplex* der[3];
  double* out[3];
  double* soap;
  double* sqrt_dot;
  int2* index;
};

static device_problem upload_problem(const problem& p) {
  device_problem d;
  d.k2_i_site = upload(p.k2_i_site);
  d.skip = (bool*) upload(p.skip);
  d.multiplicity = upload(p.multiplicity);
  d.cnk = upload(p.cnk);
  d.rjs = upload(p.rjs);
  d.rad = upload(p.rad);
  d.rad_der = upload(p.rad_der);
  for (int i = 0; i < 4; i++) d.ang[i] = upload(p.ang[i]);
  for (int i = 0; i < 3; i++) {
    gpuErrchk(hipMalloc((void**) &d.der[i], sizeof(hipDoubleComplex) * p.n_atom_pairs * p.k_max * p.n_max));
    gpuErrchk(hipMalloc((void**) &d.out[i], sizeof(double) * p.n_atom_pairs * p.n_soap));
  }
  gpuErrchk(hipMalloc((void**) &d.soap, sizeof(double) * p.n_sites * p.n_soap));
  gpuErrchk(hipMalloc((void**) &d.sqrt_dot, sizeof(double) * p.n_sites));
  gpuErrchk(hipMalloc((void**) &d.index, sizeof(int2) * p.n_max * p.n_max * (p.l_max + 1)));
  return d;
}

static void free_problem(device_problem& d) {
  hipFree(d.k2_i_site);
  hipFree(d.skip);
  hipFree(d.multiplicity);
  hipFree(d.cnk);
  hipFree(d.rjs);
  hipFree(d.rad);
  hipFree(d.rad_der);
  for (int i = 0; i < 4; i++) hipFree(d.ang[i]);
  for (int i = 0; i < 3; i++) {
    hipFree(d.der[i]);
    hipFree(d.out[i]);
  }
  hipFree(d.soap);
  hipFree(d.sqrt_dot);
  hipFree(d.index);
}

// Pairs past the cutoff are left unwritten by the old kernel; zeroed here so
// they read as the zero the new one forms.
static void reference_der(const problem& p, device_problem& d, hipStream_t* s) {
  size_t n_der = (size_t) p.n_atom_pairs * p.k_max * p.n_max;
  for (int i = 0; i < 3; i++) gpuErrchk(hipMemsetAsync(d.der[i], 0, sizeof(hipDoubleComplex) * n_der, *s));
  reference_get_derivatives<<<dim3((p.n_atom_pairs + tpbcnk - 1) / tpbcnk, p.n_max, p.k_max), tpbcnk, 0, *s>>>(
      d.rad, d.ang[0], d.rad_der, d.ang[1], d.ang[2], d.ang[3], d.der[0], d.der[1], d.der[2], d.rjs, rcut_max, p.n_atom_pairs,
      p.n_sites, p.n_soap, p.k_max, p.n_max, p.l_max);
  reference_soap_der_one(nullptr, nullptr, nullptr, d.multiplicity, d.out[0], d.out[1], d.out[2], d.cnk, d.der[0], d.der[1],
                         d.der[2], d.k2_i_site, d.skip, p.n_sites, p.n_atom_pairs, p.n_soap, p.k_max, p.n_max, p.l_max, s);
}

static void new_der(const problem& p, device_problem& d, hipStream_t* s) {
  cuda_get_soap_index(d.index, d.skip, p.n_max, p.l_max, s);
  cuda_get_soap_der_one(d.multiplicity, d.out[0], d.out[1], d.out[2], d.cnk, d.rad, d.rad_der, d.ang[0], d.ang[1], d.ang[2],
                        d.ang[3], d.rjs, rcut_max, d.k2_i_site, d.index, p.n_sites, p.n_atom_pairs, p.k_max, p.n_max, p.l_max,
                        s);
}

static void reference_p(const problem& p, device_problem& d, hipStream_t* s) {
  reference_soap_p(d.soap, d.sqrt_dot, d.multiplicity, d.cnk, d.skip, p.n_sites, p.n_soap, p.n_max, p.l_max, s);
}

static void new_p(const problem& p, device_problem& d, hipStream_t* s) {
  cuda_get_soap_index(d.index, d.skip, p.n_max, p.l_max, s);
  cuda_get_soap_p(d.soap, d.sqrt_dot, d.multiplicity, d.cnk, d.index, p.n_sites, p.n_soap, p.n_max, p.l_max, s);
}

static std::vector<double> fetch(const double* dev, size_t n) {
  std::vector<double> h(n);
  gpuErrchk(hipMemcpy(h.data(), dev, sizeof(double) * n, hipMemcpyDeviceToHost));
  return h;
}

static std::vector<double> der_outputs(const problem& p, const device_problem& d) {
  size_t n = (size_t) p.n_atom_pairs * p.n_soap;
  std::vector<double> h;
  for (int i = 0; i < 3; i++) {
    std::vector<double> part = fetch(d.out[i], n);
    h.insert(h.end(), part.begin(), part.end());
  }
  return h;
}

static std::vector<double> p_outputs(const problem& p, const device_problem& d) {
  std::vector<double> h = fetch(d.soap, (size_t) p.n_sites * p.n_soap);
  std::vector<double> norm = fetch(d.sqrt_dot, p.n_sites);
  h.insert(h.end(), norm.begin(), norm.end());
  return h;
}

static void poison(const problem& p, device_problem& d) {
  for (int i = 0; i < 3; i++) gpuErrchk(hipMemset(d.out[i], 0xff, sizeof(double) * p.n_atom_pairs * p.n_soap));
  gpuErrchk(hipMemset(d.soap, 0xff, sizeof(double) * p.n_sites * p.n_soap));
  gpuErrchk(hipMemset(d.sqrt_dot, 0xff, sizeof(double) * p.n_sites));
}

typedef void (*kernel_fn)(const problem&, device_problem&, hipStream_t*);
typedef std::vector<double> (*outputs_fn)(const problem&, const device_problem&);

static float time_ms(kernel_fn fn, const problem& p, device_problem& d, hipStream_t* s, int reps) {
  hipEvent_t a, b;
  hipEventCreate(&a);
  hipEventCreate(&b);
  fn(p, d, s);
  hipEventRecord(a, *s);
  for (int r = 0; r < reps; r++) fn(p, d, s);
  hipEventRecord(b, *s);
  hipEventSynchronize(b);
  float ms;
  hipEventElapsedTime(&ms, a, b);
  return ms / reps;
}

static long differing_words(const std::vector<double>& ref, const std::vector<double>& got) {
  long n = 0;
  for (size_t i = 0; i < ref.size(); i++)
    if (std::memcmp(&ref[i], &got[i], sizeof(double)) != 0) n++;
  return n;
}

// Returns the number of differing words; prints them and the timings.
static long compare_kernel(const char* name, kernel_fn ref_fn, kernel_fn new_fn, outputs_fn out, const problem& p,
                           device_problem& d, hipStream_t* s, int reps) {
  poison(p, d);
  ref_fn(p, d, s);
  gpuErrchk(hipStreamSynchronize(*s));
  std::vector<double> ref = out(p, d);
  poison(p, d);
  new_fn(p, d, s);
  gpuErrchk(hipStreamSynchronize(*s));
  std::vector<double> got = out(p, d);
  long n_diff = differing_words(ref, got);
  double max_rel = 0.0;
  for (size_t i = 0; i < ref.size(); i++) {
    double scale = std::max(std::fabs(ref[i]), 1e-300);
    max_rel = std::max(max_rel, std::fabs(got[i] - ref[i]) / scale);
  }
  float t_ref = time_ms(ref_fn, p, d, s, reps);
  float t_new = time_ms(new_fn, p, d, s, reps);
  std::printf("      %-8s differing %ld of %zu (max rel %.1e)  ref %8.3f ms  new %8.3f ms  x%.2f\n", name, n_diff, ref.size(),
              max_rel, t_ref, t_new, t_ref / t_new);
  return n_diff;
}

// The chain gpu_soap_pair_contract replaces, on the host from der_one's output:
// normalisation, Cartesian, the central pair as minus the others, the dot with
// G. Returns the largest error over the largest value.
static double check_contract(const problem& p, device_problem& d, hipStream_t* s) {
  int ns = p.n_sites;
  int npairs = p.n_atom_pairs;
  int nsoap = p.n_soap;
  int neigh = npairs / ns;
  std::mt19937 rng(777);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::vector<double> thetas(npairs);
  std::vector<double> phis(npairs);
  std::vector<double> G((size_t) ns * nsoap);
  std::vector<int> k2_start(ns);
  std::vector<int> n_neigh(ns, neigh);
  for (int k2 = 0; k2 < npairs; k2++) {
    thetas[k2] = std::acos(u(rng));
    phis[k2] = M_PI * u(rng);
  }
  for (double& g : G) g = u(rng);
  for (int i = 0; i < ns; i++) k2_start[i] = i * neigh;

  new_p(p, d, s);
  new_der(p, d, s);
  gpuErrchk(hipStreamSynchronize(*s));
  std::vector<double> soap = fetch(d.soap, (size_t) ns * nsoap);
  std::vector<double> norm = fetch(d.sqrt_dot, ns);
  std::vector<double> der[3];
  for (int c = 0; c < 3; c++) der[c] = fetch(d.out[c], (size_t) npairs * nsoap);

  std::vector<double> ref(3 * (size_t) npairs, 0.0);
  for (int k2 = 0; k2 < npairs; k2++) {
    if (k2 % neigh == 0)
      continue;
    int i = k2 / neigh;
    double nrm = norm[i];
    double pd[3] = {0.0, 0.0, 0.0};
    for (int q = 0; q < nsoap; q++)
      for (int c = 0; c < 3; c++) pd[c] += soap[q + (size_t) i * nsoap] * der[c][k2 + (size_t) q * npairs];
    double f[3] = {0.0, 0.0, 0.0};
    for (int q = 0; q < nsoap; q++) {
      double pq = soap[q + (size_t) i * nsoap];
      for (int c = 0; c < 3; c++)
        f[c] += G[i + (size_t) q * ns] * (der[c][k2 + (size_t) q * npairs] / nrm - pq * pd[c] / (nrm * nrm * nrm));
    }
    double st = std::sin(thetas[k2]);
    double ct = std::cos(thetas[k2]);
    double sp = std::sin(phis[k2]);
    double cp = std::cos(phis[k2]);
    double r = p.rjs[k2];
    ref[3 * k2] = st * cp * f[0] - ct * cp / r * f[2] - sp / r * f[1];
    ref[3 * k2 + 1] = st * sp * f[0] - ct * sp / r * f[2] + cp / r * f[1];
    ref[3 * k2 + 2] = ct * f[0] + st / r * f[2];
  }
  for (int i = 0; i < ns; i++)
    for (int c = 0; c < 3; c++)
      for (int j = 1; j < neigh; j++) ref[3 * (i * neigh) + c] -= ref[3 * (i * neigh + j) + c];

  std::vector<double> q_host(soap);
  for (int i = 0; i < ns; i++)
    for (int q = 0; q < nsoap; q++) q_host[q + (size_t) i * nsoap] /= norm[i];
  double* q_d = upload(q_host);
  double* G_d = upload(G);
  double* thetas_d = upload(thetas);
  double* phis_d = upload(phis);
  int* k2_start_d = upload(k2_start);
  int* n_neigh_d = upload(n_neigh);
  double* pair_d;
  gpuErrchk(hipMalloc((void**) &pair_d, sizeof(double) * 3 * npairs));
  gpu_soap_keep_expansion(q_d, d.sqrt_dot, d.cnk, d.rad, d.rad_der, d.ang[0], d.ang[1], d.ang[2], d.ang[3], d.rjs, thetas_d,
                          phis_d, d.k2_i_site, k2_start_d, n_neigh_d, d.multiplicity, d.skip, rcut_max, ns, npairs, nsoap,
                          p.k_max, p.n_max, p.l_max, s);
  gpu_soap_pair_contract(G_d, 1.0, pair_d, s);
  gpuErrchk(hipStreamSynchronize(*s));
  std::vector<double> got = fetch(pair_d, 3 * (size_t) npairs);
  hipFree(kept.index_d);
  kept = {};
  for (void* x : {(void*) q_d, (void*) G_d, (void*) thetas_d, (void*) phis_d, (void*) k2_start_d, (void*) n_neigh_d, (void*) pair_d})
    hipFree(x);

  double err = 0.0;
  double big = 0.0;
  for (size_t k = 0; k < ref.size(); k++) {
    err = std::max(err, std::fabs(got[k] - ref[k]));
    big = std::max(big, std::fabs(ref[k]));
  }
  std::printf("      contract max error / max value %.1e\n", err / big);
  return err / big;
}

// cuda_get_soap_index against the walk it replaced, on the host.
static long check_index(const problem& p, device_problem& d, hipStream_t* s) {
  int n_max = p.n_max;
  int l_max = p.l_max;
  cuda_get_soap_index(d.index, d.skip, n_max, l_max, s);
  gpuErrchk(hipStreamSynchronize(*s));
  std::vector<int2> got(n_max * n_max * (l_max + 1));
  gpuErrchk(hipMemcpy(got.data(), d.index, sizeof(int2) * got.size(), hipMemcpyDeviceToHost));
  long bad = 0;
  int counter = 0;
  int counter2 = 0;
  for (int n = 1; n <= n_max; n++)
    for (int np = n; np <= n_max; np++)
      for (int l = 0; l <= l_max; l++) {
        int q = l + (l_max + 1) * (np - 1 + (n - 1) * n_max);
        int2 want = p.skip[q] ? make_int2(-1, 0) : make_int2(counter, counter2);
        if (!p.skip[q]) {
          counter++;
          counter2 += l + 1;
        }
        if (got[q].x != want.x || got[q].y != want.y)
          bad++;
      }
  std::printf("      index    differing %ld\n", bad);
  return bad;
}

static long check_case(const char* label, int n_sites, int neigh, int n_max, int l_max, double skip_fraction, int reps,
                       hipStream_t* s) {
  problem p = make_problem(n_sites, neigh, n_max, l_max, skip_fraction, 12345);
  device_problem d = upload_problem(p);
  std::printf("    %s: sites %d pairs %d n_soap %d\n", label, p.n_sites, p.n_atom_pairs, p.n_soap);
  long bad = 0;
  bad += check_index(p, d, s);
  bad += compare_kernel("soap_p", reference_p, new_p, p_outputs, p, d, s, reps);
  bad += compare_kernel("der_one", reference_der, new_der, der_outputs, p, d, s, reps);
  if (!(check_contract(p, d, s) < CONTRACT_TOL))
    bad++;
  free_problem(d);
  return bad;
}

int main(int argc, char** argv) {
  int n_sites = argc > 1 ? std::atoi(argv[1]) : 1000;
  int neigh = argc > 2 ? std::atoi(argv[2]) : 40;
  int reps = argc > 3 ? std::atoi(argv[3]) : 5;
  hipStream_t s;
  gpuErrchk(hipStreamCreate(&s));
  long bad = 0;
  bad += check_case("n_max 16 l_max 8", n_sites, neigh, 16, 8, 0.0, reps, &s);
  bad += check_case("n_max 16 l_max 8, 30% skip", n_sites, neigh, 16, 8, 0.3, reps, &s);
  bad += check_case("n_max 5 l_max 3", n_sites, neigh, 5, 3, 0.0, reps, &s);
  bad += check_case("n_max 11 l_max 6, 50% skip", n_sites, neigh, 11, 6, 0.5, reps, &s);
  bad += check_case("n_max 24 l_max 10", n_sites / 4, neigh, 24, 10, 0.0, reps, &s);
  std::printf(bad == 0 ? "    PASS: bit-identical, contraction within tolerance\n" : "    FAIL: outputs differ\n");
  return bad == 0 ? 0 : 1;
}
