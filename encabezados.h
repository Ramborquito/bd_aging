# include <cuda.h>
#include <curand_kernel.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <vector>

# include <gsl/gsl_blas.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# pragma once
# define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)
# define PI 3.14159266f
# define CUTOFF_FAC 1.12246205f
# define PROB 0.2f
# define RAN rannew64
# define EPS_MIN 1.0e-6f
# define EPS 1.0f
# define N_NEIGHS_MAX 15000




typedef struct 
{
  float *force_ls;
} 
return_vectors;


struct Grains
{
    int number_particles;
    float3 *rr;
    float3 *rr0;        // for msd
    float3 *rr_raw;     // for msd
    float3 *vv;
    float3 *ff;
    float *diameter;
    float *mass;
    float *virial;
    float *potential;
};

struct OverlapParameters {
    float side;
    float cell_side;
    float kappa;
    float gamma;
    float v0;
    float overlap_max;
    float time;
    float dt;
    int idum;
    int ngrain_tot;
    int ncell;
    int ntags;
};

struct SystemParameters {

    int ngrain_total;
    int number_species;
    int *N;     //number of particles of each specie
    float side;
    float dt;
    float temperature;
    float bin_size_dplt;
    float bin_size_gder;
    float range_dplt;
    float range_gder;
    int nbins_dplt;
    int nbins_gder;
    int not_finished;
    int NH;
    int WCA_flag;
    int idum;
};

struct PolydispersityParameters {

    float *phi;
    float *diameter;
    float *polydispersity;
    float *max_polydispersity;

};

struct cell{


};

//====================================================================================

float RAN(long *);

void set_vec_int_hst(int *, int, int);
void set_vec_float_hst(float *, int, float);
void set_vec_float3_hst(float3 *, int, float3);

void cell_locate_hst(char, float3 *, int *, int *, SystemParameters);

void update_ermak_hst(Grains &grain_vec, SystemParameters pars, gsl_rng *rand);

void get_forces_same_hst(Grains &, SystemParameters);

void get_forces_diff_hst(Grains &, Grains &, SystemParameters);

void potential_wca49_50_hst(float sigma, float dist_inv, float invT,float AA,
                            float &potential, float &normal_force);

void potential_wca49_50_AO_hst(float sigma, float dist_inv, float invT,float AA, float sigma_pol, float phi_pol,
                               float &potential, float &normal_force);

void potential_wca_hst(float sigma, float dist_inv, float &potential, float &normal_force);

void depletion(Grains, SystemParameters, float *);
int svdcmp(float **, int, int, float *, float **);

void get_gder_hst(int, int, float3 *, float3 *, float *, SystemParameters);

float calculate_temp_linear(float T0, float Tf, float t0, float tf, float t);
float calculate_temp_sine(float T0, float Tf, float t0, float tf, float period, float t);

void fself_isotropic(float3 *rr_raw, float3 *rr0, int ngrain, float *fself, int sample_number, float qmax);

// ===================================================================================
__device__ void potential_wca49_50_dev(float sigma, float dist_inv, float invT,float AA,
                                       float &potential, float &normal_force);

__device__ void potential_wca49_50_AO_dev(float sigma, float dist_inv, float invT,float AA, float sigma_pol, float phi_pol,
                                          float &potential, float &normal_force);

__device__ void potential_wca_dev(float sigma, float dist_inv, float &potential, float &normal_force);

__global__ void set_vec_int_dev(int *, int, int); 
__global__ void set_vec_float_dev(float *, int, float); 
__global__ void set_vec_float3_dev(float3 *, int, float3); 

__global__ void cell_locate_dev(char, float3 *, int *, int *, SystemParameters);

__global__ void setup_rng_kernel(curandState *states, unsigned long seed, int N);
__global__ void update_ermak_dev(Grains &grain_vec, SystemParameters pars, curandState *states);

__global__ void get_forces_same_dev(Grains &, SystemParameters);

__global__ void get_forces_diff_dev(Grains &, Grains &, SystemParameters);

__global__ void get_gder_dev(int, int, const float3 *, const float3 *, float *, SystemParameters);
