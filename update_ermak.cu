#include "encabezados.h"

// ================================================== 2023-05-22 =====================
// HOST			(Ermak-McCammon)
// ===================================================================================


void update_ermak_hst(char type, float3 *rr_vec, float3 *rr_raw_vec,
                      float3 *ff_vec, gsl_rng *rand, parametros pars) {
    int ngrain;
    float3 drr_random, drr_total, ff, rr, rr_raw;
    float side, side_inv, dt, diameter, diameter_inverse;
    float temperature, amplitude;

    side = pars.side;
    side_inv = 1.0f / side;
    temperature = pars.temp_set;
    dt = pars.dt;

    if (type == 'b') {
        ngrain = pars.ngrain_big;
        diameter = pars.sigma_big;
        diameter_inverse = 1.0f / diameter;
    } else {
        ngrain = pars.ngrain_sml;
        diameter = pars.sigma_sml;
        diameter_inverse = 1.0f / diameter;
    }

    // amplitude of the random displacement
    amplitude = sqrt(2.0 * temperature * diameter_inverse * dt);

    for (int mm = 0; mm < ngrain; ++mm) {
        //fetch
        ff = ff_vec[mm];
        rr = rr_vec[mm];
        rr_raw = rr_raw_vec[mm];

        //random displacement

        drr_random.x = amplitude * gsl_ran_gaussian(rand, 1.0);
        drr_random.y = amplitude * gsl_ran_gaussian(rand, 1.0);
        drr_random.z = amplitude * gsl_ran_gaussian(rand, 1.0);

        // added force displacement

        drr_total.x = ff.x * diameter_inverse * dt + drr_random.x;
        drr_total.y = ff.y * diameter_inverse * dt + drr_random.y;
        drr_total.z = ff.z * diameter_inverse * dt + drr_random.z;

        //move position

        rr.x += drr_total.x;
        rr.y += drr_total.y;
        rr.z += drr_total.z;

        rr_raw.x += drr_total.x;
        rr_raw.y += drr_total.y;
        rr_raw.z += drr_total.z;

        //apply periodic boundary conditions

        rr.x -= side * floor(side_inv * rr.x);
        rr.y -= side * floor(side_inv * rr.y);
        rr.z -= side * floor(side_inv * rr.z);

        if (rr.x < 0.0001f && rr.y < 0.0001f && rr.z < 0.0001f)
            printf("Warning: Particle position at zero: (%f, %f, %f)\n", rr.x, rr.y, rr.z);
        //save

        rr_vec[mm] = rr;
        rr_raw_vec[mm] = rr_raw;

    }


}

// ===================================================================================
// DEVICE 		(Ermak-McCammon)
// ===================================================================================

#include <curand_kernel.h>

__global__ void setup_rng_kernel(curandState *states, unsigned long seed, int N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < N)
        curand_init(seed, idx, 0, &states[idx]);
}

__global__ void update_ermak_dev(char type, float3 *rr_vec, float3 *rr_raw_vec,
                      float3 *ff_vec, curandState *states, parametros pars) {
    int ngrain, mm;
    float3 drr_random, drr_total, ff, rr, rr_raw;
    float side, side_inv, dt, diameter, diameter_inverse;
    float temperature, amplitude;

    side = pars.side;
    side_inv = 1.0f / side;
    temperature = pars.temp_set;
    dt = pars.dt;

    if (type == 'b') {
        ngrain = pars.ngrain_big;
        diameter = pars.sigma_big;
        diameter_inverse = 1.0f / diameter;
    } else {
        ngrain = pars.ngrain_sml;
        diameter = pars.sigma_sml;
        diameter_inverse = 1.0f / diameter;
    }

    // amplitude of the random displacement
    amplitude = sqrt(2.0 * temperature * diameter_inverse * dt);

    mm = threadIdx.x + blockIdx.x * blockDim.x;

    if (mm < ngrain) {
        //fetch
        ff = ff_vec[mm];
        rr = rr_vec[mm];
        rr_raw = rr_raw_vec[mm];

        //random displacement

        drr_random.x = amplitude * curand_normal(&states[mm]);
        drr_random.y = amplitude * curand_normal(&states[mm]);
        drr_random.z = amplitude * curand_normal(&states[mm]);

        // added force displacement

        drr_total.x = ff.x * diameter_inverse * dt + drr_random.x;
        drr_total.y = ff.y * diameter_inverse * dt + drr_random.y;
        drr_total.z = ff.z * diameter_inverse * dt + drr_random.z;

        //move position

        rr.x += drr_total.x;
        rr.y += drr_total.y;
        rr.z += drr_total.z;

        rr_raw.x += drr_total.x;
        rr_raw.y += drr_total.y;
        rr_raw.z += drr_total.z;

        //apply periodic boundary conditions

        rr.x -= side * floor(side_inv * rr.x);
        rr.y -= side * floor(side_inv * rr.y);
        rr.z -= side * floor(side_inv * rr.z);

        //save

        rr_vec[mm] = rr;
        rr_raw_vec[mm] = rr_raw;

    }


}
