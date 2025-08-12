//
// Created by marcoramirez on 8/7/25.
//

#include "encabezados.h"
# define FACTOR 2.0


bool between(const float3 rr1, const float3 rr2, const float3 *rr_sml_vec, parametros pars) {

    //unitary vector
    float3 u, rr, dr1, dr2, uxrr;
    float u_abs, dr1Norm, dr2Norm, dr_line;
    float max_distance_big = 0.5 * pars.sigma_big + FACTOR * pars.sigma_sml;
    float max_distance_line = 0.35 * pars.sigma_big;
    float side, side_inv;
    side = pars.side;
    side_inv = 1 / side;

    u.x = rr1.x - rr2.x;
    u.y = rr1.y - rr2.y;
    u.z = rr1.z - rr2.z;
    // periodic conditions
    // periodic boundary conditions
    u.x -= side * floor(side_inv * u.x + 0.5f);
    u.y -= side * floor(side_inv * u.y + 0.5f);
    u.z -= side * floor(side_inv * u.z + 0.5f);

    u_abs = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    u.x /= u_abs;
    u.y /= u_abs;
    u.z /= u_abs;

    for (int i = 0; i < pars.ngrain_sml; ++i) {
        rr = rr_sml_vec[i];
        dr1.x = rr.x - rr1.x;
        dr1.y = rr.y - rr1.y;
        dr1.z = rr.z - rr1.z;
        // periodic boundary conditions
        dr1.x -= side * floor(side_inv * dr1.x + 0.5f);
        dr1.y -= side * floor(side_inv * dr1.y + 0.5f);
        dr1.z -= side * floor(side_inv * dr1.z + 0.5f);
        dr1Norm = sqrt(dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z);
        if (dr1Norm > max_distance_big) continue;
        dr2.x = rr.x - rr2.x;
        dr2.y = rr.y - rr2.y;
        dr2.z = rr.z - rr2.z;
        // periodic boundary conditions
        dr2.x -= side * floor(side_inv * dr2.x + 0.5f);
        dr2.y -= side * floor(side_inv * dr2.y + 0.5f);
        dr2.z -= side * floor(side_inv * dr2.z + 0.5f);
        dr2Norm = sqrt(dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z);
        if (dr2Norm > max_distance_big) continue;
        // calculate distance to the line
        uxrr.x = dr1.y * u.z - dr1.z * u.y;
        uxrr.y = -dr1.x * u.z + dr1.z * u.x;
        uxrr.z = dr1.x * u.y - dr1.y * u.x;
        dr_line = sqrt(uxrr.x * uxrr.x + uxrr.y * uxrr.y + uxrr.z * uxrr.z);
        if (dr_line > max_distance_line) continue;
        return true;
    }
    return false;
}

void
calculate_color(int mm, float3 *rr_big_vec, float3 *rr_sml_vec, char *color, int &blue_count, int &green_count, int &red_count,
                parametros pars) {

    float3 drr, rrm, rrn;
    float side, side_inv, dist, sigma_big, sigma_sml;
    int ngrain_big = pars.ngrain_big;
    side = pars.side;
    side_inv = 1.0 / side;
    sigma_big = pars.sigma_big;
    sigma_sml = pars.sigma_sml;

        rrm = rr_big_vec[mm];
        int color_id = 0;
        bool some_near = false;
        for (int nn = 0; nn < ngrain_big; nn++) {
            if (nn == mm) continue;
            rrn = rr_big_vec[nn];
            drr.x = rrn.x - rrm.x;
            drr.y = rrn.y - rrm.y;
            drr.z = rrn.z - rrm.z;

            // periodic boundary conditions
            drr.x -= side * floor(side_inv * drr.x + 0.5f);
            drr.y -= side * floor(side_inv * drr.y + 0.5f);
            drr.z -= side * floor(side_inv * drr.z + 0.5f);

            dist = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
            if (dist > sigma_big + 2.0 * sigma_sml) continue;
            between(rrm, rrn, rr_sml_vec, pars) ? color_id++ : color_id--;
            some_near = true;
        }

        if (some_near) {
            if (color_id > 0) {
                sprintf(color, "Red");
                red_count++;
            } else {
                sprintf(color, "Green");
                green_count++;
            }
        } else {
            sprintf(color, "Blue");
            blue_count++;
        }

}