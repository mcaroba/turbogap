// cuda_get_soap_p and cuda_get_soap_der_one against the loop orders they
// replaced, kept below, on random input: every output must be bit-identical,
// with and without skipped components. Prints the time of each on the same data.
#include "../../src/gpu/gap_soap_descriptor.cu"

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
  std::vector<hipDoubleComplex> der[3];
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
  size_t n_der = (size_t) p.n_atom_pairs * p.k_max * n_max;
  for (size_t i = 0; i < n_cnk; i++) p.cnk.push_back(make_hipDoubleComplex(u(rng), u(rng)));
  for (int d = 0; d < 3; d++)
    for (size_t i = 0; i < n_der; i++) p.der[d].push_back(make_hipDoubleComplex(u(rng), u(rng)));
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
  for (int i = 0; i < 3; i++) {
    d.der[i] = upload(p.der[i]);
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
  for (int i = 0; i < 3; i++) {
    hipFree(d.der[i]);
    hipFree(d.out[i]);
  }
  hipFree(d.soap);
  hipFree(d.sqrt_dot);
  hipFree(d.index);
}

static void reference_der(const problem& p, device_problem& d, hipStream_t* s) {
  reference_soap_der_one(nullptr, nullptr, nullptr, d.multiplicity, d.out[0], d.out[1], d.out[2], d.cnk, d.der[0], d.der[1],
                         d.der[2], d.k2_i_site, d.skip, p.n_sites, p.n_atom_pairs, p.n_soap, p.k_max, p.n_max, p.l_max, s);
}

static void new_der(const problem& p, device_problem& d, hipStream_t* s) {
  cuda_get_soap_index(d.index, d.skip, p.n_max, p.l_max, s);
  cuda_get_soap_der_one(d.multiplicity, d.out[0], d.out[1], d.out[2], d.cnk, d.der[0], d.der[1], d.der[2], d.k2_i_site, d.index,
                        p.n_sites, p.n_atom_pairs, p.k_max, p.n_max, p.l_max, s);
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
  float t_ref = time_ms(ref_fn, p, d, s, reps);
  float t_new = time_ms(new_fn, p, d, s, reps);
  std::printf("      %-8s differing %ld of %zu  ref %8.3f ms  new %8.3f ms  x%.2f\n", name, n_diff, ref.size(), t_ref, t_new,
              t_ref / t_new);
  return n_diff;
}

static long check_case(const char* label, int n_sites, int neigh, int n_max, int l_max, double skip_fraction, int reps,
                       hipStream_t* s) {
  problem p = make_problem(n_sites, neigh, n_max, l_max, skip_fraction, 12345);
  device_problem d = upload_problem(p);
  std::printf("    %s: sites %d pairs %d n_soap %d\n", label, p.n_sites, p.n_atom_pairs, p.n_soap);
  long bad = 0;
  bad += compare_kernel("soap_p", reference_p, new_p, p_outputs, p, d, s, reps);
  bad += compare_kernel("der_one", reference_der, new_der, der_outputs, p, d, s, reps);
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
  std::printf(bad == 0 ? "    PASS: bit-identical\n" : "    FAIL: outputs differ\n");
  return bad == 0 ? 0 : 1;
}
