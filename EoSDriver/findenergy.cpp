#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "NuclearEos.h"


void find_energy_bdry( double x, double y, double z, double *found_leps, double *alltables_mode,
                       int nx, int ny, int nz, int neps, double *xt, double *yt, double *zt, double *logeps,
                       int keymode, int *keyerr );


//-------------------------------------------------------------------------------------
// Function    :  find_energy
// Description :  Finding energy from temperature (1)
//                                    entropy     (2)
//                                    pressure    (3) mode
//                Using 3D Catmull-Rom cubic interpolation formula
//                to search the corresponding energy given (rho, (T, e, P), Y_e)
//
// Parameter   :  x              : input vector of first  variable (rho)
//                y              : input vector of second variable (e, T, P)
//                z              : input vector of third  variable (Y_e)
//                found_leps     : output log(eps) of interpolated function values
//                alltables_mode : 3d array of tabulated logenergy
//                nx             : x-dimension of table
//                ny             : y-dimension of table
//                nz             : z-dimension of table
//                neps           : size of energy array in the Nuclear Eos table
//                xt             : vector of x-coordinates of table
//                yt             : vector of y-coordinates of table
//                zt             : vector of z-coordinates of table
//                logeps         : log(eps) array in the table
//                keymode        : which mode we will use
//                                 1: temperature mode (coming in with T)
//                                 2: entropy mode     (coming in with entropy)
//                                 3: pressure mode    (coming in with P)
//                keyerr         : output error
//-------------------------------------------------------------------------------------
void find_energy( double x, double y, double z, double *found_leps, double *alltables_mode,
                  int nx, int ny, int nz, int neps, double *xt, double *yt, double *zt, double *logeps,
                  int keymode, int *keyerr )
{


# define CUBE(x) ((x)*(x)*(x))
# define SQR(x) ((x)*(x))

   double   dx, dy, dz, dxi, dyi, dzi;
   double   delx, dely, delz;
   double   *pv;
   double   u[4], v[4], w[4];
   double   r[4], q[4];
   double   vox;
   int      ix, iy, iz;
   int      nxy;

   nxy = nx*ny;

// determine spacing parameters of equidistant (!!!) table
#if 1
   dx = ( xt[nx-1] - xt[0] ) / ( 1.0*(nx-1) );
   dy = ( yt[ny-1] - yt[0] ) / ( 1.0*(ny-1) );
   dz = ( zt[nz-1] - zt[0] ) / ( 1.0*(nz-1) );
   
   dxi = 1.0/dx;
   dyi = 1.0/dy;
   dzi = 1.0/dz;
#endif

#if 0
   dx = drho;
   dy = deps;
   dz = dye;
   dxi = drhoi;
   dyi = depsi;
   dzi = dyei;
#endif

   // determine location in table
   ix = (int)( (x - xt[0] + 1.0e-10) * dxi );
   iy = (int)( (y - yt[0] + 1.0e-10) * dyi );
   iz = (int)( (z - zt[0] + 1.0e-10) * dzi );
   
   if ( ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz )
   {
      *keyerr = 667; // out of range
      return;
   }
   
   // linear interpolation at boundaries
   if ( ix == 0 || iy == 0 || iz == 0 ||
        ix == nx-2 || iy == ny-2 || iz == nz-2 )
   {
      find_energy_bdry( x, y, z, found_leps, alltables_mode,
                        nx, ny, nz, neps, xt, yt, zt, logeps, keymode, keyerr );
      return;
   }

   // differences
   delx = ( x - xt[ix] ) * dxi;
   dely = ( y - yt[iy] ) * dyi;
   delz = ( z - zt[iz] ) * dzi;
   
   // factors for Catmull-Rom interpolation
   u[0] = -0.5*CUBE(delx) + SQR(delx) - 0.5*delx;
   u[1] = 1.5*CUBE(delx) - 2.5*SQR(delx) + 1.0;
   u[2] = -1.5*CUBE(delx) + 2.0*SQR(delx) + 0.5*delx;
   u[3] = 0.5*CUBE(delx) - 0.5*SQR(delx);
   
   v[0] = -0.5*CUBE(dely) + SQR(dely) - 0.5*dely;
   v[1] = 1.5*CUBE(dely) - 2.5*SQR(dely) + 1.0;
   v[2] = -1.5*CUBE(dely) + 2.0*SQR(dely) + 0.5*dely;
   v[3] = 0.5*CUBE(dely) - 0.5*SQR(dely);
   
   w[0] = -0.5*CUBE(delz) + SQR(delz) - 0.5*delz;
   w[1] = 1.5*CUBE(delz) - 2.5*SQR(delz) + 1.0;
   w[2] = -1.5*CUBE(delz) + 2.0*SQR(delz) + 0.5*delz;
   w[3] = 0.5*CUBE(delz) - 0.5*SQR(delz);
   
   int iv;
   if      ( keymode == 1 ) iv = 0; // energy table for the temperature mode
   else if ( keymode == 2 ) iv = 1; // energy table for the entropy mode
   else if ( keymode == 3 ) iv = 2; // energy table for the pressure mode
   
   vox = 0.0;
    
   pv = alltables_mode + iv + 3 * ((ix-1) + (iy-1)*nx + (iz-1)*nxy);
      for (int k=0; k<4; k++) 
      {
         q[k] = 0.0;
         for (int j=0; j<4; j++)
         {
            r[j] = 0.0;
            for (int i=0; i<4; i++)
            {
               r[j] += u[i]* *pv;
               pv += 3;
            }
            q[k] += v[j]*r[j];
            pv += 3*nx - 4*3;
         }
         vox += w[k]*q[k];
         pv += nxy*3 - 4*3*nx;
      }
   *found_leps = vox;


   // linear interpolation when cubic interpolations failed
   if ( isnan(vox) )
   {
      find_energy_bdry( x, y, z, found_leps, alltables_mode,
                        nx, ny, nz, neps, xt, yt, zt, logeps, keymode, keyerr );
      return;
   }


  return;
} // FUMCTIOM : find_energy


//-------------------------------------------------------------------------------------
// Function    :  find_energy_bdry
// Description :  Finding energy from temperature (1)
//                                    entropy     (2)
//                                    pressure    (3) mode
//                Using linear interpolation at boundaries of table
//                to search the corresponding energy given
//                (rho, (T, e, P), Y_e)
//
// Parameter   :  x              : input vector of first  variable (rho)
//                y              : input vector of second variable (e, T, P)
//                z              : input vector of third  variable (Y_e)
//                found_leps     : output log(eps) of interpolated function values
//                alltables_mode : 3d array of tabulated logenergy
//                nx             : x-dimension of table
//                ny             : y-dimension of table
//                nz             : z-dimension of table
//                neps           : size of energy array in the Nuclear Eos table
//                xt             : vector of x-coordinates of table
//                yt             : vector of y-coordinates of table
//                zt             : vector of z-coordinates of table
//                logeps         : log(eps) array in the table
//                keymode        : which mode we will use
//                                 1: temperature mode (coming in with T)
//                                 2: entropy mode     (coming in with entropy)
//                                 3: pressure mode    (coming in with P)
//                keyerr         : output error
//-------------------------------------------------------------------------------------
void find_energy_bdry( double x, double y, double z, double *found_leps, double *alltables_mode,
                       int nx, int ny, int nz, int neps, double *xt, double *yt, double *zt, double *logeps,
                       int keymode, int *keyerr )
{


  // helper variables
  double d1, d2, d3;
  double fh[8], delx, dely, delz, a[8];
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
   
   idx[0] = 3 * (ix + nx * (iy+ny*iz));
   idx[1] = 3 * ((ix-1) + nx * (iy + ny*iz));
   idx[2] = 3 * (ix + nx * ((iy-1) + ny*iz));
   idx[3] = 3 * (ix + nx * (iy + ny * (iz-1))) ;
   idx[4] = 3 * ((ix-1) + nx * ((iy-1) + ny*iz));
   idx[5] = 3 * ((ix-1)  + nx * (iy + ny * (iz-1)));
   idx[6] = 3 * (ix + nx * ((iy-1) + ny * (iz-1))) ;
   idx[7] = 3 * ((ix-1) + nx * ((iy-1) + ny * (iz-1)));
   
   int iv;
   if      ( keymode == 1 ) iv = 0;
   else if ( keymode == 2 ) iv = 1;
   else if ( keymode == 3 ) iv = 2;
   // set up aux vars for interpolation
   // assuming array ordering (iv, ix, iy, iz)
   
   fh[0] = alltables_mode[iv+idx[0]];
   fh[1] = alltables_mode[iv+idx[1]];
   fh[2] = alltables_mode[iv+idx[2]];
   fh[3] = alltables_mode[iv+idx[3]];
   fh[4] = alltables_mode[iv+idx[4]];
   fh[5] = alltables_mode[iv+idx[5]];
   fh[6] = alltables_mode[iv+idx[6]];
   fh[7] = alltables_mode[iv+idx[7]];

   // set up coeffs of interpolation polynomical and
   // evaluate function values
   
   a[0] = fh[0];
   a[1] = dxi   * ( fh[1] - fh[0] );
   a[2] = dyi   * ( fh[2] - fh[0] );
   a[3] = dzi   * ( fh[3] - fh[0] );
   a[4] = dxyi  * ( fh[4] - fh[1] - fh[2] + fh[0] );
   a[5] = dxzi  * ( fh[5] - fh[1] - fh[3] + fh[0] );
   a[6] = dyzi  * ( fh[6] - fh[2] - fh[3] + fh[0] );
   a[7] = dxyzi * ( fh[7] - fh[0] + fh[1] + fh[2] + 
                    fh[3] - fh[4] - fh[5] - fh[6] );
   
   *found_leps = a[0] + a[1] * delx
               + a[2] * dely
               + a[3] * delz
               + a[4] * delx * dely
               + a[5] * delx * delz
               + a[6] * dely * delz
               + a[7] * delx * dely * delz;
   
   if ( isnan(*found_leps) || !( *found_leps>logeps[0] && *found_leps<logeps[neps-1] ) )
   {
      *keyerr = 667;
   }


  return;

} // FUNCTION : find_energy_bdry
