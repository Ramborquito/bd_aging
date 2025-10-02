#include "encabezados.h"

// =============================================== 2023-05-26 ========================
// 				HOST
// ===================================================================================

void get_gder_hst(int specie_1, int specie_2, float3 *rr_vec, float3 *rr_other_vec,
                  float *gder_vec, SystemParameters pars)
{
  float3 rrm, rrn, drr;
  float dist, range_gder, bin_size, side, side_inv;
  int NH, blockidxx, ngrain, ngrain_other, bin, bin_index, nbins_gder;

  side = pars.side;
  range_gder = pars.range_gder;
  nbins_gder = pars.nbins_gder;
  NH = pars.NH;   // NH = numero de hilos en un bloque
  side_inv = 1.0f/side;

  bin_size = range_gder/((float) nbins_gder);

  // colects distances

  for (int mm = 0; mm < pars.N[specie_1]; mm++)
  {
    // fetch

    rrm = rr_vec[mm];

    // initialize

    blockidxx = mm/NH;

    for (int nn = 0; nn < pars.N[specie_2]; nn++)
    {
      // cartesian distances, modulo side

      if (specie_1 == specie_2 && mm == nn) continue;
      rrn = rr_other_vec[nn];
 
      drr.x = rrm.x - rrn.x;
      drr.y = rrm.y - rrn.y;
      drr.z = rrm.z - rrn.z;
      drr.x -= side*floorf(side_inv*drr.x + 0.5);
      drr.y -= side*floorf(side_inv*drr.y + 0.5);
      drr.z -= side*floorf(side_inv*drr.z + 0.5);
      dist = sqrt(drr.x*drr.x + drr.y*drr.y + drr.z*drr.z);
      if (dist < range_gder)
      {
        bin = (int) (dist/bin_size);
        if (bin == nbins_gder) bin--;
        bin_index = blockidxx*nbins_gder + bin;  // "blocks" with own histogram
        gder_vec[bin_index] += 1.0f;
      }
    }
  } 
  return;
}

// ===================================================================================
// 				DEVICE
// ===================================================================================

__global__ void get_gder_dev(int specie_1, int specie_2, const float3 *rr_vec,
                             const float3 *rr_other_vec, float *gder_vec, SystemParameters pars)
{
  float3 rrm, rrn, drr;
  float dist, range_gder, bin_size, side, side_inv, aux; 
  int mm, bin, nbins_gder, bin_index;



  side = pars.side;
  range_gder = pars.range_gder;
  nbins_gder = pars.nbins_gder;
  side_inv = 1.0f/side;

  bin_size = range_gder/((float) nbins_gder);

  // get thread index

  mm = threadIdx.x + blockIdx.x*blockDim.x; 

  // colects distances

  if (mm < pars.N[specie_1])
  {
    // fetch

    rrm = rr_vec[mm];

    // run over other particles

    for (int nn = 0; nn < pars.N[specie_2]; nn++)
    {
      // cartesian distances, modulo side

      if (specie_1 == specie_2 && mm == nn) continue;

      rrn = rr_other_vec[nn];
 
      drr.x = rrm.x - rrn.x;
      drr.y = rrm.y - rrn.y;
      drr.z = rrm.z - rrn.z;
      drr.x -= side*floorf(side_inv*drr.x + 0.5);
      drr.y -= side*floorf(side_inv*drr.y + 0.5);
      drr.z -= side*floorf(side_inv*drr.z + 0.5);
      dist = sqrt(drr.x*drr.x + drr.y*drr.y + drr.z*drr.z);
      if (dist < range_gder)
      {
        bin = (int) (dist/bin_size);
        if (bin == nbins_gder) bin--;
        bin_index = blockIdx.x*nbins_gder + bin;  // "blocks" with own histogram
        aux = atomicAdd(&(gder_vec[bin_index]), 1.0f);
      }
    }
  } 
  return;
}
