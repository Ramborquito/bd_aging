#include "encabezados.h"

// =============================================== 2023-05-26 ========================
// 				HOST
// ===================================================================================

void get_gder_hst(char type, char type_other, float3 *rr_vec, float3 *rr_other_vec, 
         float *gder_vec, parametros pars)
{
  static float3 rrm, rrn, drr;
  static float dist, range_gder, bin_size, side, side_inv; 
  static int mm, nn, NH, blockidxx, ngrain, ngrain_other, bin, bin_index, nbins_gder;

  if (type == 'b') ngrain = pars.ngrain_big;
  else ngrain = pars.ngrain_sml;
  if (type_other == 'b') ngrain_other = pars.ngrain_big;
  else ngrain_other = pars.ngrain_sml;

  side = pars.side;
  range_gder = pars.range_gder;
  nbins_gder = pars.nbins_gder;
  NH = pars.NH;   // NH = numero de hilos en un bloque
  side_inv = 1.0/side;

  bin_size = range_gder/((float) nbins_gder);

  // colects distances

  for (mm = 0; mm < ngrain; mm++)
  {
    // fetch

    rrm = rr_vec[mm];

    // initialize

    blockidxx = mm/NH;

    for (nn = 0; nn < ngrain_other; nn++)
    {
      // cartesian distances, modulo side

      if (type == type_other && mm == nn) continue;
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

__global__ void get_gder_dev(char type, char type_other, float3 *rr_vec, 
                    float3 *rr_other_vec, float *gder_vec, parametros pars)
{
  float3 rrm, rrn, drr;
  float dist, range_gder, bin_size, side, side_inv, aux; 
  int mm, nn, ngrain, ngrain_other, bin, nbins_gder, bin_index;

  if (type == 'b') ngrain = pars.ngrain_big;
  else ngrain = pars.ngrain_sml;
  if (type_other == 'b') ngrain_other = pars.ngrain_big;
  else ngrain_other = pars.ngrain_sml;

  side = pars.side;
  range_gder = pars.range_gder;
  nbins_gder = pars.nbins_gder;
  side_inv = 1.0/side;

  bin_size = range_gder/((float) nbins_gder);

  // get thread index

  mm = threadIdx.x + blockIdx.x*blockDim.x; 

  // colects distances

  if (mm < ngrain)
  {
    // fetch

    rrm = rr_vec[mm];

    // run over other particles

    for (nn = 0; nn < ngrain_other; nn++)
    {
      // cartesian distances, modulo side

      if (type == type_other && mm == nn) continue;

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
