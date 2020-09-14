#include "NuclearEos.h"


//-------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_linterp_some
// Description :  Finding thermodynamic varibales
//                using linear interpolation
//                by searching the tabulated nuclear EoS
//
// Parameter   :  x           : input vector of first  variable (rho)
//                y           : input vector of second variable (eps)
//                z           : input vector of third  variable (Y_e)
//                output_vars : output variables of interpolated function values
//                alltables   : 3d array of tabulated variables
//                nx          : x-dimension of table
//                ny          : y-dimension of table
//                nz          : z-dimension of table
//                nvars       : number of variables we will find
//                xt          : vector of x-coordinates of table
//                yt          : vector of y-coordinates of table
//                zt          : vector of z-coordinates of table
//-------------------------------------------------------------------------------------
void nuc_eos_C_linterp_some( double x, double y, double z,
                             double *output_vars, double *alltables,
                             int nx, int ny, int nz, int nvars,
                             double *xt, double *yt, double *zt ) 
{

   // helper variables
   double d1, d2, d3;
   double fh[8][nvars], delx, dely, delz, a[8];
   double dx, dy, dz, dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi;
   int n, ix, iy, iz;
   

// determine spacing parameters of equidistant (!!!) table
#if 1
   dx = ( xt[nx-1] - xt[0] ) / ( 1.0*(nx-1) );
   dy = ( yt[ny-1] - yt[0] ) / ( 1.0*(ny-1) );
   dz = ( zt[nz-1] - zt[0] ) / ( 1.0*(nz-1) );
   
   dxi = 1.0 / dx;
   dyi = 1.0 / dy;
   dzi = 1.0 / dz;
#endif

#if 0
   dx = drho;
   dy = deps;
   dz = dye;

   dxi = drhoi;
   dyi = depsi;   
   dzi = dyei;    
#endif            
                  
   dxyi = dxi * dyi;
   dxzi = dxi * dzi;
   dyzi = dyi * dzi;
               
   dxyzi = dxi * dyi * dzi;
               
   // determine location in table
               
   ix = 1 + (int)( (x - xt[0] - 1.0e-10) * dxi );
   iy = 1 + (int)( (y - yt[0] - 1.0e-10) * dyi );
   iz = 1 + (int)( (z - zt[0] - 1.0e-10) * dzi );
               
   ix = MAX( 1, MIN( ix, nx-1 ) );
   iy = MAX( 1, MIN( iy, ny-1 ) );
   iz = MAX( 1, MIN( iz, nz-1 ) );
   
   // set up aux vars for interpolation
   
   delx = xt[ix] - x;
   dely = yt[iy] - y;
   delz = zt[iz] - z;
   
   int idx[8];
   
   idx[0] = NTABLES*(ix + nx*(iy + ny*iz));
   idx[1] = NTABLES*((ix-1) + nx*(iy + ny*iz));
   idx[2] = NTABLES*(ix + nx*((iy-1) + ny*iz));
   idx[3] = NTABLES*(ix + nx*(iy + ny*(iz-1)));
   idx[4] = NTABLES*((ix-1) + nx*((iy-1) + ny*iz));
   idx[5] = NTABLES*((ix-1) + nx*(iy + ny*(iz-1)));
   idx[6] = NTABLES*(ix + nx*((iy-1) + ny*(iz-1)));
   idx[7] = NTABLES*((ix-1) + nx*((iy-1) + ny*(iz-1)));
   
   
   for(int iv=0;iv<nvars;iv++) 
   {
     // set up aux vars for interpolation
     // assuming array ordering (iv, ix, iy, iz)
   
      fh[0][iv] = alltables[iv+idx[0]];
      fh[1][iv] = alltables[iv+idx[1]];
      fh[2][iv] = alltables[iv+idx[2]];
      fh[3][iv] = alltables[iv+idx[3]];
      fh[4][iv] = alltables[iv+idx[4]];
      fh[5][iv] = alltables[iv+idx[5]];
      fh[6][iv] = alltables[iv+idx[6]];
      fh[7][iv] = alltables[iv+idx[7]];
   
      // set up coeffs of interpolation polynomical and
      // evaluate function values
   
      a[0] = fh[0][iv];
      a[1] = dxi *   ( fh[1][iv] - fh[0][iv] );
      a[2] = dyi *   ( fh[2][iv] - fh[0][iv] );
      a[3] = dzi *   ( fh[3][iv] - fh[0][iv] );
      a[4] = dxyi *  ( fh[4][iv] - fh[1][iv] - fh[2][iv] + fh[0][iv] );
      a[5] = dxzi *  ( fh[5][iv] - fh[1][iv] - fh[3][iv] + fh[0][iv] );
      a[6] = dyzi *  ( fh[6][iv] - fh[2][iv] - fh[3][iv] + fh[0][iv] );
      a[7] = dxyzi * ( fh[7][iv] - fh[0][iv] + fh[1][iv] + fh[2][iv] + 
                       fh[3][iv] - fh[4][iv] - fh[5][iv] - fh[6][iv] );
   
      output_vars[iv] = a[0] + a[1] * delx
                      + a[2] * dely
                      + a[3] * delz
                      + a[4] * delx * dely
                      + a[5] * delx * delz
                      + a[6] * dely * delz
                      + a[7] * delx * dely * delz;
   }


   return;

} // FUNCTION : nuc_eos_linterp_some
