#include "encabezados.h"

void fself_isotropic(float3 *rr_raw, float3 *rr0, int ngrain, float *fself, int sample_number, float qmax) {

    float sum = 0.0f;
    // Loop over each grain
    for (int i = 0; i < ngrain; i++) {
        float3 drr;
        drr.x = rr_raw[i].x - rr0[i].x;
        drr.y = rr_raw[i].y - rr0[i].y;
        drr.z = rr_raw[i].z - rr0[i].z;

        float drr_magnitude = sqrt(drr.x * drr.x + drr.y * drr.y + drr.z * drr.z);
        sum += sin(qmax * drr_magnitude) / (qmax * drr_magnitude); // Sinc function


    }

    sum /= ngrain;
    fself[sample_number] += sum;

}