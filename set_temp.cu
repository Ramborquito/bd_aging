#include "encabezados.h"

// ===================================================================================
//   ONLY HOST
// ===================================================================================

float calculate_temp_linear(float T0, float Tf, float t0, float tf, float t){
    if(t < t0) return T0;
    if(t > tf) return Tf;
    float m = (Tf-T0)/(tf-t0);
    return T0 + m * (t - t0);
}

float calculate_temp_sine(float T0, float Tf, float t0, float tf, float period, float t){
    // a + b sin [2pi ( t / tau - 1 / 4 )]
    float a, b, argument;
    if(t < t0) return T0;
    if(t > tf) return Tf;
    a = 0.5 * (T0 + Tf);
    b = 0.5 * (Tf - T0);
    argument = 2 * PI * (t / period - 0.25);

    return a + b * sin(argument);
}