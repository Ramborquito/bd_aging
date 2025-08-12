//
// Created by marcoramirez on 5/8/25.
//

#ifndef CONTRACTION_VECTOR_CONTRACTING_VECTOR3D_H
#define CONTRACTION_VECTOR_CONTRACTING_VECTOR3D_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <optional>
#include <cstring>

# include <gsl/gsl_blas.h>
# include <gsl/gsl_vector.h>
# include <gsl/gsl_matrix.h>
# include <gsl/gsl_linalg.h>

// Habilita o deshabilita validaciones
//#define DEBUG_CHECKS

#ifdef DEBUG_CHECKS
#define DEBUG_ASSERT(cond, msg)                                         \
        do {                                                                \
            if (!(cond)) {                                                 \
                std::cerr << "DEBUG ASSERTION FAILED: " << msg             \
                          << "\nFile: " << __FILE__                         \
                          << "\nLine: " << __LINE__ << std::endl;          \
                std::abort();                                               \
            }                                                               \
        } while (0)
#else
#define DEBUG_ASSERT(cond, msg) do {} while (0)
#endif

constexpr int NEIGHS_MAX = 10000;
constexpr float MIN_EPS = 1.0e-6f;

class Contracting_Vector3D {

public:

    Contracting_Vector3D(int bins, float range);

    void add_statistics(float3 *rr_vec, float3 *xx_vec, int nparticles, std::optional<float> side = std::nullopt);

    void calculate_vector();

    void get_vector(float *xx_contracted_vec) const;


private:

    int samples;
    bool calculated;
    int bins;
    float range;
    float bin_size;
    float3 *atb_vec;
    float3 *hit_vec;
    float3 *drr_vec;
    float *dist_vec;
    float **ata_mat_x;
    float **ata_mat_y;
    float **ata_mat_z;

    // auxiliar
    float bin_size_inv;
    float range_square;

    //result
    float *contracted_vec;


};


#endif //CONTRACTION_VECTOR_CONTRACTING_VECTOR2D_H
