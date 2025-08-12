#include "encabezados.h"

// ======================================= 2023-05-22 ================================
//   HOST
// ===================================================================================

void set_vec_int_hst(int *vec, int size, int value)
{
  static int mm;

  // reset vec

  for (mm = 0; mm < size; mm++) vec[mm] = value;

  return;
}

// ===================================================================================

void set_vec_float_hst(float *vec, int size, float value)
{
  static int mm;

  // reset vec

  for (mm = 0; mm < size; mm++) vec[mm] = value;

  return;
}

// ===================================================================================

void set_vec_float3_hst(float3 *vec, int size, float3 value)
{
  static int mm;

  // reset vec

  for (mm = 0; mm < size; mm++) vec[mm] = value;

  return;
}

// ===================================================================================

void cell_locate_hst(char type, float3 *rr_vec, int *nocup_vec, int *cell_vec,
         parametros pars)
{
  static float3 rr;
  static float cell_side_inv;
  static int mm, ii, jj, kk, cell_index, ngrain, ntags, ncell, shift;

  // escoje celdas

  if (type == 'b')
  {
    cell_side_inv = 1.0/pars.cell_side_big;
    ngrain = pars.ngrain_big;
    ntags = pars.ntags_big;
    ncell = pars.ncell_big;
  }
  else
  {
    cell_side_inv = 1.0/pars.cell_side_sml;
    ngrain = pars.ngrain_sml;
    ntags = pars.ntags_sml;
    ncell = pars.ncell_sml;
  }

  // distribuye particulas en celdas

  for (mm = 0; mm < ngrain; mm++)
  {
    // fetch 

    rr = rr_vec[mm];

    // get cell location

    ii = (int) (cell_side_inv*rr.x);
    jj = (int) (cell_side_inv*rr.y);
    kk = (int) (cell_side_inv*rr.z);
    if (ii == ncell) ii--;  // 100. - 1.0x10^-6 = 100. Ouch!
    if (jj == ncell) jj--;  // ditto
    if (kk == ncell) kk--;  // ditto
    cell_index = ii + ncell*(jj + ncell*kk);
    
    // add one to nocup 

    shift = nocup_vec[cell_index];
    if (shift == ntags-1) 
    {
      printf("error: cell_index %d  shift %d  ntags %d\n", cell_index, shift, ntags);
      exit (1);
    }
    nocup_vec[cell_index]++;

    // save grain index in cell vector

    cell_vec[cell_index*ntags + shift] = mm;
  }
  
  return;
}

// ===================================================================================
//   DEVICE
// ===================================================================================

__global__ void set_vec_int_dev(int *vec, int size, int value)
{
  int mm;

  mm = threadIdx.x + blockIdx.x*blockDim.x;

  // reset vec

  if (mm < size) vec[mm] = value;

  return;
}

// ===================================================================================

__global__ void set_vec_float_dev(float *vec, int size, float value)
{
  int mm;

  mm = threadIdx.x + blockIdx.x*blockDim.x;

  // reset vec

  if (mm < size) vec[mm] = value;

  return;
}

// ===================================================================================

__global__ void set_vec_float3_dev(float3 *vec, int size, float3 value)
{
  int mm;

  mm = threadIdx.x + blockIdx.x*blockDim.x;

  // reset vec

  if (mm < size) vec[mm] = value;

  return;
}

// ===================================================================================

__global__ void cell_locate_dev(char type, float3 *rr_vec_dev, int *nocup_vec_dev,
	            int *cell_vec_dev, parametros pars)
{
  float3 rr;
  float cell_side_inv;
  int mm, ii, jj, kk, cell_index, ngrain, ncell, ntags, shift;

  // escoje entre celdas big o sml

  if (type == 'b')
  {
    cell_side_inv = 1.0/pars.cell_side_big;
    ngrain = pars.ngrain_big;
    ncell = pars.ncell_big;
    ntags = pars.ntags_big;
  }
  else
  {
    cell_side_inv = 1.0/pars.cell_side_sml;
    ngrain = pars.ngrain_sml;
    ncell = pars.ncell_sml;
    ntags = pars.ntags_sml;
  }

  mm = threadIdx.x + blockIdx.x*blockDim.x;

  // fill up cells

  if (mm < ngrain)
  {
    // fetch

    rr = rr_vec_dev[mm];

    // get cell location. first reduce to unit cell

    ii = (int) (cell_side_inv*rr.x);
    jj = (int) (cell_side_inv*rr.y);
    kk = (int) (cell_side_inv*rr.z);
    if (ii == ncell) ii--;
    if (jj == ncell) jj--;
    if (kk == ncell) kk--;
    cell_index = ii + ncell*(jj + ncell*kk);
    
    // add one to nocup 

    shift = atomicAdd(&(nocup_vec_dev[cell_index]), 1);

    // save grain index in cell vector

    cell_vec_dev[cell_index*ntags + shift] = mm;
  }
  
  return;
}
