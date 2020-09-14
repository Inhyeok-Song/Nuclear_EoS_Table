#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "nuc_eos.h"


void cubinterp_find(int ix, int iz, double xrho, double xye, double f_in, double *f, double *fs
                    int keytemp, int *keyerr, double rfeps) {

  int iv;
  
  // determine spacing parameters of equidistant tabldx = (
  
  
  // set mode (energy, temp, entropy, pressure)
  if (keytemp == 0) {
    iv == 1;
  }
  else if (keytemp == 2) {
    iv == 2;
  }
  else if (keytemp == 3) {
    iv == 3; 
  }
  
  double fa[ntemp-1];
  double xval;
  
  double lrho = log10(xrho);
  
  if (ix>1 || ix<nrho-1 || iz>1 || iz<nye-1) {
    for (int j = 1; j < ntemp-1; j++) {
      // cubic interpolation in Ye
      double fv[4]  = {0.0}; // array for bicubic interpolation     
      double fv2[4] = {0.0}; // array for bicubic interpolation
      for (int ip = 0; ip < 4; ip++) {
  
        fv[0]   = f[iv + NTABLES*((ix-1+ip) + nrho*(j + ntemp*(iz-1)))];
        fv[1]   = f[iv + NTABLES*((ix-1+ip) + nrho*(j + ntemp*iz))];
        fv[2]   = f[iv + NTABLES*((ix-1+ip) + nrho*(j + ntemp*(iz+1)))];
        fv[3]   = f[iv + NTABLES*((ix-1+ip) + nrho*(j + ntemp*(iz+2)))];
  
        xval    = (xye - yes[iz])*dyei;
  
        fv2[ip] =   (-0.5*fv[0] + 1.5*fv[1] - 1.5*fv[2] + 0.5*fv[3])*xval**3
                  + ( 1.0*fv[0] - 2.5*fv[1] + 2.0*fv[2] - 0.5*fv[3])*xval**2
                  + (-0.5*fv[0]             + 0.5*fv[2]            )*xval
                  + fv[1];
      }
      // cubic interpolation in ye
      x_val  = (lrho - logrho[ix])*drhoi;
      fa[j]  =   (-0.5*fv2[0] + 1.5*fv2[1] - 1.5*fv2[2] + 0.5*fv2[3])*xval**3
               + ( 1.0*fv2[0] - 2.5*fv2[1] + 2.0*fv2[2] - 0.5*fv2[3])*xval**2
               + (-0.5*fv2[0]              + 0.5*fv2[2]             )*xval
               + fv2[1];
    }
  // linear interpolation at boundaries
    double fv[8]  = {0.0}; // array for bilinear interpolation     
    double fv2[4] = {0.0}; // array for bilinear interpolation
  
    fv[0] = f[iv + NTABLES*(ix + nrho*(ntemp*iz))];
    fv[1] = f[iv + NTABLES*(ix + nrho*(ntemp*(iz+1)))];
    fv[2] = f[iv + NTABLES*(ix+1 + nrho*(ntemp*iz))];
    fv[3] = f[iv + NTABLES*(ix+1 + nrho*(ntemp*(iz+1)))];
  
    fv[4] = f[iv + NTABLES*(ix + nrho*(ntemp-1 + ntemp*iz))];
    fv[5] = f[iv + NTABLES*(ix + nrho*(ntemp-1 + ntemp*(iz+1)))];
    fv[6] = f[iv + NTABLES*(ix+1 + nrho*(ntemp-1 + ntemp*iz))];
    fv[7] = f[iv + NTABLES*(ix+1 + nrho*(ntemp-1 + ntemp*(iz+1)))];
  
    fv2[0] = (fv[1] - fv[0])*dyei*(xye - yes[iz]);
    fv2[1] = (fv[3] - fv[2])*dyei*(xye - yes[iz]);
  
    fv2[2] = (fv[5] - fv[4])*dyei*(xye - yes[iz]);
    fv2[3] = (fv[7] - fv[6])*dyei*(xye - yes[iz]);
  
    fa[0]       = (fv2[1] - fv2[0])*drhoi*(lrho - logrho[ix]);
    fa[ntemp-2] = (fv2[1] - fv2[0])*drhoi*(lrho - logrho[ix]);
  
  }
  else {
    for (int j = 0, j < ntemp-1, j++) {
      // linear interpolation
      double fv[4]  = {0.0}; // array for bilinear interpolation     
      double fv2[2] = {0.0}; // array for bilinear interpolation
  
      fv[0] = f[iv + NTABLES*(ix + nrho*(j + ntemp*iz))];
      fv[1] = f[iv + NTABLES*(ix + nrho*(j + ntemp*(iz+1)))];
      fv[2] = f[iv + NTABLES*(ix+1 + nrho*(j + ntemp*iz))];
      fv[3] = f[iv + NTABLES*(ix+1 + nrho*(j + ntemp*(iz+1)))];
  
      fv2[0] = (fv[1] - fv[0])*dyei*(xye - yes[iz]);
      fv2[1] = (fv[3] - fv[2])*dyei*(xye - yes[iz]);
  
      fa[j]  = (fv2[1] - fv2[0])*drhoi*(lrho - logrho[ix]);
    
    }
  }


  // root finding
  int iy;
  for (int j = 0; j < ntemp-1; j++) {
    if (f_in > fa[j] && f_in < fa[j+1]) {
      iy = j + 1;
    }
    else {
      keyerrt == 668;
    }
  }
  if (iy > 1 || iy < ntemp-1) {
    for (int j = 1; j < ntemp-1; j++) {
      // cubic interpolation in Ye
      double fv[4]  = {0.0}; // array for bicubic interpolation     
      for (int ip = 0; ip < 4; ip++) {
  
        fv[0]   = fa[iy-1]
        fv[1]   = fa[iy];
        fv[2]   = fa[iy+1];
        fv[3]   = fa[iy+2];
  
        xval    = ( - fv[1])*dtempi;
  
        fv2[ip] =   (-0.5*fv[0] + 1.5*fv[1] - 1.5*fv[2] + 0.5*fv[3])*xval**3
                  + ( 1.0*fv[0] - 2.5*fv[1] + 2.0*fv[2] - 0.5*fv[3])*xval**2
                  + (-0.5*fv[0]             + 0.5*fv[2]            )*xval
                  + fv[1];
 

 
  
  return;
}
