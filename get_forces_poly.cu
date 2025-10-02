#include "encabezados.h"

// =============================================== 2023-05-26 ========================
// HOST
// ===================================================================================

float Apow_hst(float x) {
    float x25;
    float x8;
    float xx;
    xx = x * x;
    xx = xx * xx;
    x8 = xx * xx;
    xx = x8 * x8;
    x25 = x * xx * x8;
    return x25 * x25;
}

void get_forces_same_hst(Grains &grain_vec, SystemParameters pars) {
    float3 rrm, rrn, drr, ff_pair, ffm;
    float cutoff2, side, side_inv, cutoff, virial, potential, sigma, dist_inv, dist2, normal_force;

    //wca_49-50
//  float invT=0.678575;
//  float BB=50.0/49.0;
//  float AA= 50*Apow_dev(BB);
//  cutoff=BB*sigma;

    side = pars.side;
    side_inv = 1.0f / side;

    // calcula interacciones

    for (int mm = 0; mm < grain_vec.number_particles; mm++) {
        // fetch

//        if (mm == 0) printf("is running get_forces_same_hst\n");


        rrm = grain_vec.rr[mm];
        ffm = grain_vec.ff[mm];
        virial = grain_vec.virial[mm];
        potential = grain_vec.potential[mm];

        for (int nn = 0; nn < grain_vec.number_particles; ++nn) {

            if (nn == mm) continue;

            // fecht another

            rrn = grain_vec.rr[nn];

            // calculate sigma

            sigma = 0.5f * (grain_vec.diameter[mm] + grain_vec.diameter[nn]);
            cutoff = 1.122462048f * sigma;
            cutoff2 = cutoff * cutoff;

            // coordinate differences

            drr.x = rrn.x - rrm.x;
            drr.y = rrn.y - rrm.y;
            drr.z = rrn.z - rrm.z;

            // periodic boundary conditions

            drr.x -= side * floor(side_inv * drr.x + 0.5f);    // per. bound. cond.
            drr.y -= side * floor(side_inv * drr.y + 0.5f);
            drr.z -= side * floor(side_inv * drr.z + 0.5f);

            // distance and normal_force

            dist2 = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
            if (dist2 < cutoff2) {
                dist_inv = sqrt(1.0f / dist2);
                //potential_wca49_50_dev(sigma,dist_inv,invT,AA,potential,normal_force);
                potential_wca_hst(sigma, dist_inv, potential, normal_force);

                // calcula las fuerzas normales para este par

                ff_pair.x = normal_force * dist_inv * drr.x;
                ff_pair.y = normal_force * dist_inv * drr.y;
                ff_pair.z = normal_force * dist_inv * drr.z;

                // suma a la fuerza, al virial y al potencial. se hace una doble suma

                ffm.x += ff_pair.x;
                ffm.y += ff_pair.y;
                ffm.z += ff_pair.z;

                virial -= 0.5 * (drr.x * ff_pair.x + drr.y * ff_pair.y + drr.z * ff_pair.z);
            }
        }

        // save

        grain_vec.ff[mm] = ffm;
        grain_vec.virial[mm] = virial;
        grain_vec.potential[mm] = potential;
    }
    return;
}

// ===================================================================================

void
get_forces_diff_hst(Grains &grain_vec, Grains &grain_other_vec, SystemParameters pars) {
    float3 rrm, rrn, drr, ff_pair, ffm;
    float cutoff2, side, side_inv, cutoff, virial, potential, sigma, dist_inv, dist2, normal_force;

    side = pars.side;
    side_inv = 1.0f / side;

    //wca_49-50
//  float invT=0.678575;
//  float BB=50.0/49.0;
//  float AA= 50*Apow_dev(BB);
//  cutoff=BB*sigma;

    //wca
    cutoff = 1.122462048f * sigma;

    cutoff2 = cutoff * cutoff;
    side = pars.side;
    side_inv = 1.0f / side;

    // calcula interacciones sobre distintos

    for (int mm = 0; mm < grain_vec.number_particles; mm++) {

//        if (mm == 0) printf("is running get_forces_diff_hst\n");

        // fetch

        rrm = grain_vec.rr[mm];
        ffm = grain_vec.ff[mm];
        virial = grain_vec.virial[mm];
        potential = grain_vec.potential[mm];

        // checa con granos del otro taman~o

        for (int nn = 0; nn < grain_other_vec.number_particles; ++nn) {

            // fetch again

            rrn = grain_other_vec.rr[nn];

            // calculate sigma

            sigma = 0.5f * (grain_vec.diameter[mm] + grain_other_vec.diameter[nn]);
            cutoff = 1.122462048f * sigma;
            cutoff2 = cutoff * cutoff;

            // coordinate differences

            drr.x = rrn.x - rrm.x;
            drr.y = rrn.y - rrm.y;
            drr.z = rrn.z - rrm.z;

            // periodic boundary conditions

            drr.x -= side * floor(side_inv * drr.x + 0.5f);    // per. bound. cond.
            drr.y -= side * floor(side_inv * drr.y + 0.5f);
            drr.z -= side * floor(side_inv * drr.z + 0.5f);

            // distance and normal force

            dist2 = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
            if (dist2 < cutoff2) {
                dist_inv = sqrt(1.0f / dist2);
                //potential_wca49_50_dev(sigma,dist_inv,invT,AA,potential,normal_force);
                potential_wca_hst(sigma, dist_inv, potential, normal_force);

                // calcula las fuerzas normales para este par

                ff_pair.x = normal_force * dist_inv * drr.x;
                ff_pair.y = normal_force * dist_inv * drr.y;
                ff_pair.z = normal_force * dist_inv * drr.z;

                // suma a la fuerza, al virial y el potencial. se hace una doble suma

                ffm.x += ff_pair.x;
                ffm.y += ff_pair.y;
                ffm.z += ff_pair.z;

                virial -= 0.5 * (drr.x * ff_pair.x + drr.y * ff_pair.y + drr.z * ff_pair.z);
            }
        }


        // save

        grain_vec.ff[mm] = ffm;
        grain_vec.virial[mm] = virial;
        grain_vec.potential[mm] = potential;
    }
    return;
}

// ===================================================================================
// DEVICE
// ===================================================================================

__device__ float Apow_dev(float x) {
    float x25;
    float x8;
    float xx;
    xx = x * x;
    xx = xx * xx;
    x8 = xx * xx;
    xx = x8 * x8;
    x25 = x * xx * x8;
    return x25 * x25;
}

// ===================================================================================

__global__ void get_forces_same_dev(Grains &grain_vec, SystemParameters pars) {
    float3 rrm, rrn, drr, ff_pair, ffm;
    float cutoff2, side, side_inv, cutoff, virial, potential, sigma, dist_inv, dist2, normal_force;
    float diameter_mm, diameter_nn;
    int mm;

    side = pars.side;
    side_inv = 1.0f / side;

    // calcula interacciones

    mm = threadIdx.x + blockIdx.x * blockDim.x;
    // printf("is running get_forces_same_dev\n");

    if (mm < grain_vec.number_particles) {

        // fetch

        rrm = grain_vec.rr[mm];
        ffm = grain_vec.ff[mm];
        virial = grain_vec.virial[mm];
        potential = grain_vec.potential[mm];
        diameter_mm = grain_vec.diameter[mm];

        for (int nn = 0; nn < grain_vec.number_particles; ++nn) {
            if (nn == mm) continue;

            // another fetch

            rrn = grain_vec.rr[nn];

            // coordinate differences

            drr.x = rrn.x - rrm.x;
            drr.y = rrn.y - rrm.y;
            drr.z = rrn.z - rrm.z;

            // periodic boundary conditions

            drr.x -= side * floor(side_inv * drr.x + 0.5f);    // per. bound. cond.
            drr.y -= side * floor(side_inv * drr.y + 0.5f);
            drr.z -= side * floor(side_inv * drr.z + 0.5f);

            // distance

            dist2 = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;

            if(dist2 > 2.0f) continue;  // skip pairs that are too far away for short range interactions

            // calculate sigma

            diameter_nn = grain_vec.diameter[nn];
            sigma = 0.5f * (diameter_mm + diameter_nn);
            cutoff = 1.122462048f * sigma;
            cutoff2 = cutoff * cutoff;

            if (dist2 < cutoff2) {
                dist_inv = sqrt(1.0f / dist2);
                //potential_wca49_50_dev(sigma,dist_inv,invT,AA,potential,normal_force);
                potential_wca_dev(sigma, dist_inv, potential, normal_force);

                // calcula las fuerzas normales para este par

                ff_pair.x = normal_force * dist_inv * drr.x;
                ff_pair.y = normal_force * dist_inv * drr.y;
                ff_pair.z = normal_force * dist_inv * drr.z;

                // suma al virial y al potencial. el 2. se acepta suma doble en potencial

                ffm.x += ff_pair.x;
                ffm.y += ff_pair.y;
                ffm.z += ff_pair.z;
                /*if (mm == 0)
                    printf("get_forces_same_dev: ff_pair.x: %f, ff_pair.y: %f, ff_pair.z: %f\n", ff_pair.x, ff_pair.y,
                            ff_pair.z);*/

                virial -= 0.5f * (drr.x * ff_pair.x + drr.y * ff_pair.y + drr.z * ff_pair.z);
            }
        }
        // save

        grain_vec.ff[mm] = ffm;
        grain_vec.virial[mm] = virial;
        grain_vec.potential[mm] = potential;
    }

    return;
}

// ===================================================================================

__global__ void get_forces_diff_dev(Grains &grain_vec, Grains &grain_other_vec, SystemParameters pars) {
    float3 rrm, rrn, drr, ff_pair, ffm;
    float cutoff2, side, side_inv, cutoff, virial, potential, sigma, dist_inv, dist2, normal_force;
    float diameter_mm, diameter_nn;
    int mm;


    side = pars.side;
    side_inv = 1.0f / side;

    // calcula interacciones sobre grandes

    mm = threadIdx.x + blockIdx.x * blockDim.x;


    if (mm < grain_vec.number_particles) {
        // fetch

        rrm = grain_vec.rr[mm];
        ffm = grain_vec.ff[mm];
        virial = grain_vec.virial[mm];
        potential = grain_vec.potential[mm];
        diameter_mm = grain_vec.diameter[mm];

        // checa con granos de taman~o diferente

        for (int nn = 0; nn < grain_other_vec.number_particles; ++nn) {

            // another fetch

            rrn = grain_other_vec.rr[nn];

            // coordinate differences

            drr.x = rrn.x - rrm.x;
            drr.y = rrn.y - rrm.y;
            drr.z = rrn.z - rrm.z;

            // periodic boundary conditions

            drr.x -= side * floor(side_inv * drr.x + 0.5f);    // per. bound. cond.
            drr.y -= side * floor(side_inv * drr.y + 0.5f);
            drr.z -= side * floor(side_inv * drr.z + 0.5f);

            // distance

            dist2 = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;

            if(dist2 > 2.0f) continue;  // skip pairs that are too far away for short range interactions

            // calculate sigma

            diameter_nn = grain_other_vec.diameter[nn];
            sigma = 0.5f * (diameter_mm + diameter_nn);
            cutoff = 1.122462048f * sigma;
            cutoff2 = cutoff * cutoff;

            if (dist2 < cutoff2) {
                dist_inv = sqrt(1.0f / dist2);
                //potential_wca49_50_dev(sigma,dist_inv,invT,AA,potential,normal_force);
                potential_wca_dev(sigma, dist_inv, potential, normal_force);

                // calcula las fuerzas normales para este par

                ff_pair.x = normal_force * dist_inv * drr.x;
                ff_pair.y = normal_force * dist_inv * drr.y;
                ff_pair.z = normal_force * dist_inv * drr.z;

                // suma a fuerzas, al potencial y al virial

                ffm.x += ff_pair.x;
                ffm.y += ff_pair.y;
                ffm.z += ff_pair.z;

                virial -= 0.5f * (drr.x * ff_pair.x + drr.y * ff_pair.y + drr.z * ff_pair.z);
            }
        }

        // save

        grain_vec.ff[mm] = ffm;
        grain_vec.virial[mm] = virial;
        grain_vec.potential[mm] = potential;
    }
    return;
}
