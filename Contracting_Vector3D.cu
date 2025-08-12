//
// Created by marcoramirez on 5/8/25.
//

#include "Contracting_Vector3D.cuh"
#include <iostream>

constexpr float3 zero={0, 0, 0};


Contracting_Vector3D::Contracting_Vector3D(int bins, float range)
        : bins(bins), range(range), samples(0), calculated(false) {


    DEBUG_ASSERT(bins >= 0, "Bins number must be greater than zero.");
    DEBUG_ASSERT(range >= 0.0f, "Range must be greater than zero.");

    bin_size = range / (float) bins;
    bin_size_inv = 1.0f / ((float) bin_size);
    range_square = range * range;

    contracted_vec = (float *) malloc(bins * sizeof(float));

    atb_vec = (float3 *) malloc(bins * sizeof(float3));
    hit_vec = (float3 *) malloc(bins * sizeof(float3));

    drr_vec = (float3 *) malloc(NEIGHS_MAX * sizeof(float3));
    dist_vec = (float *) malloc(NEIGHS_MAX * sizeof(float));


    ata_mat_x = (float **) malloc(bins * sizeof(float *));
    for (int i = 0; i < bins; ++i)
        ata_mat_x[i] =
                (float *) malloc(bins * sizeof(float));
    ata_mat_y = (float **) malloc(bins * sizeof(float *));
    for (int i = 0; i < bins; ++i)
        ata_mat_y[i] =
                (float *) malloc(bins * sizeof(float));
    ata_mat_z = (float **) malloc(bins * sizeof(float *));
    for (int i = 0; i < bins; ++i)
        ata_mat_z[i] =
                (float *) malloc(bins * sizeof(float));

    for (int i = 0; i < bins; i++)
        for (int j = 0; j < bins; j++) {
            ata_mat_x[i][j] = 0.0;
            ata_mat_y[i][j] = 0.0;
            ata_mat_z[i][j] = 0.0;
        }
    for (int i = 0; i < bins; i++) {
        contracted_vec[i] = 0.0;
        atb_vec[i] = zero;
        hit_vec[i] = zero;
    }
}

void Contracting_Vector3D::add_statistics(float3 *rr_vec, float3 *xx_vec, int nparticles, std::optional<float> side) {

    DEBUG_ASSERT(nparticles > 0, "number of particles must be greater than zero");

    // essential variables
    float3 rrm, rrn, drr, ff;
    float dist, dist2, aux;
    int neighs, nb;



    // periodic boundary conditions for simulations
    float side_inv;
    if (side.has_value()) side_inv = 1.0f / side.value();

#ifdef DEBUG_CHECKS
    if (side.has_value()) {
        DEBUG_ASSERT(side.value() > 0.0f, "side must be greater than zero");
    }
#endif

    // corre sobre particulas

    for (int mm = 0; mm < nparticles; mm++) {
        rrm = rr_vec[mm];
        ff = xx_vec[mm];

        neighs = 0;
        for (int nn = 0; nn < nparticles; nn++) {

            if (nn == mm) continue;

            // calculo de la distancia

            rrn.x = (float) rr_vec[nn].x;
            rrn.y = (float) rr_vec[nn].y;
            rrn.z = (float) rr_vec[nn].z;

            drr.x = rrm.x - rrn.x;
            drr.y = rrm.y - rrn.y;
            drr.z = rrm.z - rrn.z;

            // considers that particles are in the box [0, side]
            if (side.has_value()) {
                drr.x -= side.value() * (float) floor(side_inv * drr.x + 0.5);
                drr.y -= side.value() * (float) floor(side_inv * drr.y + 0.5);
                drr.z -= side.value() * (float) floor(side_inv * drr.z + 0.5);
            }

            dist2 = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
            if (dist2 < range_square) {
                dist = sqrt(dist2);
                dist_vec[neighs] = dist;
                drr_vec[neighs] = drr;
                neighs++;
                if (neighs == NEIGHS_MAX) {
                    printf("n_neighs %d = N_NEIGHS_MAX %d\n", neighs, NEIGHS_MAX);
                    exit(1);
                }
            }
        }


        // DEPLETION LS

        if (neighs > 0) {
            for (int ng = 0; ng < neighs; ng++)    // obtine vectores de hits
            {
                drr = drr_vec[ng];
                dist = dist_vec[ng];
                nb = (int) (bin_size_inv * dist);
                hit_vec[nb].x += (drr.x / dist);
                hit_vec[nb].y += (drr.y / dist);
                hit_vec[nb].z += (drr.z / dist);
            }

            // los acumula en ata_mat

            for (int i = 0; i < bins; i++)
                for (int j = 0; j <= i; j++) {
                    aux = hit_vec[i].x * hit_vec[j].x;
                    ata_mat_x[i][j] += aux;
                    if (i != j) ata_mat_x[j][i] += aux;

                    aux = hit_vec[i].y * hit_vec[j].y;
                    ata_mat_y[i][j] += aux;
                    if (i != j) ata_mat_y[j][i] += aux;

                    aux = hit_vec[i].z * hit_vec[j].z;
                    ata_mat_z[i][j] += aux;
                    if (i != j) ata_mat_z[j][i] += aux;
                }

            // add to atb vectors

            for (int i = 0; i < bins; i++) {
                atb_vec[i].x += hit_vec[i].x * ff.x;
                atb_vec[i].y += hit_vec[i].y * ff.y;
                atb_vec[i].z += hit_vec[i].z * ff.z;
            }

            // clear hit_vec

            for (int i = 0; i < bins; i++) hit_vec[i] = zero;
        }
    }
    samples++;
}

void Contracting_Vector3D::calculate_vector() {

    if (samples <= 0)
        throw std::runtime_error("Number of samples taken must be greater than zero.");

    gsl_vector *work, *W;
    gsl_matrix *A, *V;
    float ww, w_max, w_min, sum;

    auto *force_ls_loc_vec = (float *) malloc(bins * sizeof(float));


    auto *w_vec = (float *) malloc(bins * sizeof(float));
    auto *w_inv_vec = (float *) malloc(bins * sizeof(float));
    auto *aux_vec = (float *) malloc(bins * sizeof(float));
    auto **v_mat = (float **) malloc(bins * sizeof(float *));
    for (int i = 0; i < bins; i++)
        v_mat[i] = (float *) malloc(bins * sizeof(float));

    // memory for gsl

    work = gsl_vector_alloc(bins);
    W = gsl_vector_alloc(bins);
    A = gsl_matrix_alloc(bins, bins);
    V = gsl_matrix_alloc(bins, bins);

    // now do matrix inversion for the three f using svd for stability

    for (int dim = 0; dim < 3; dim++) {
        if (dim == 0) {
            for (int i = 0; i < bins; i++)
                for (int j = 0; j < bins; j++)
                    gsl_matrix_set(A, i, j, (double) ata_mat_x[i][j]);

            gsl_linalg_SV_decomp(A, V, W, work);

            for (int i = 0; i < bins; i++) {
                w_vec[i] = (float) gsl_vector_get(W, i);
                for (int j = 0; j < bins; j++) {
                    ata_mat_x[i][j] = (float) gsl_matrix_get(A, i, j);
                    v_mat[i][j] = (float) gsl_matrix_get(V, i, j);
                }
            }
        } else if (dim == 1) {
            for (int i = 0; i < bins; i++)
                for (int j = 0; j < bins; j++)
                    gsl_matrix_set(A, i, j, (double) ata_mat_y[i][j]);

            gsl_linalg_SV_decomp(A, V, W, work);

            for (int i = 0; i < bins; i++) {
                w_vec[i] = (float) gsl_vector_get(W, i);
                for (int j = 0; j < bins; j++) {
                    ata_mat_y[i][j] = (float) gsl_matrix_get(A, i, j);
                    v_mat[i][j] = (float) gsl_matrix_get(V, i, j);
                }
            }
        } else if (dim == 2) {
            for (int i = 0; i < bins; i++)
                for (int j = 0; j < bins; j++)
                    gsl_matrix_set(A, i, j, (double) ata_mat_z[i][j]);

            gsl_linalg_SV_decomp(A, V, W, work);

            for (int i = 0; i < bins; i++) {
                w_vec[i] = (float) gsl_vector_get(W, i);
                for (int j = 0; j < bins; j++) {
                    ata_mat_z[i][j] = (float) gsl_matrix_get(A, i, j);
                    v_mat[i][j] = (float) gsl_matrix_get(V, i, j);
                }
            }
        }

        // check singular values

        w_min = w_max = fabs(w_vec[0]);
        for (int i = 0; i < bins; i++) {
            ww = fabs(w_vec[i]);
            if (w_max < ww) w_max = ww;
            if (w_min > ww) w_min = ww;
        }

        // calcula w_inv_vec

        for (int i = 0; i < bins; i++) {
            if (fabs(w_vec[i] / w_max) > MIN_EPS) w_inv_vec[i] = 1.0f / w_vec[i];
            else w_inv_vec[i] = 0.0;
        }

        // calcula la solucion

        for (int j = 0; j < bins; j++) {
            sum = 0.0;
            if (dim == 0)
                for (int i = 0; i < bins; i++)
                    sum += ata_mat_x[i][j] * atb_vec[i].x;
            else if (dim == 1)
                for (int i = 0; i < bins; i++)
                    sum += ata_mat_y[i][j] * atb_vec[i].y;
            else if (dim == 2)
                for (int i = 0; i < bins; i++)
                    sum += ata_mat_z[i][j] * atb_vec[i].z;
            aux_vec[j] = sum;
        }
        for (int kk = 0; kk < bins; kk++) aux_vec[kk] *= w_inv_vec[kk];
        for (int i = 0; i < bins; i++) {
            sum = 0.0;
            for (int kk = 0; kk < bins; kk++) sum += v_mat[i][kk] * aux_vec[kk];
            force_ls_loc_vec[i] = sum;
        }

        // acumula force_table

        for (int nb = 0; nb < bins; nb++) {
            contracted_vec[nb] += force_ls_loc_vec[nb];
        }
    }

    //average over different dimensions
    for (int nb = 0; nb < bins; ++nb) contracted_vec[nb] /= 3;

    calculated = true;

    free(w_vec);
    free(w_inv_vec);
    free(aux_vec);
    free(force_ls_loc_vec);

    for (int i = 0; i < bins; ++i) free(v_mat[i]);
    free(v_mat);

}

void Contracting_Vector3D::get_vector(float *xx_contracted_vec) const {
    if (!calculated)
        throw std::runtime_error("Retrieved contracted vector must be firstly calculated.\n");

    memcpy(xx_contracted_vec, contracted_vec, bins * sizeof(float));
}
