#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nuc_eos.h"


void nuc_eos_C_cubinterp_some(double x, double y, double z, 
                              double *f, double *ft, 
                              int *ivs, 
                              int nx, int ny, int nz, int nvars, 
                              double *xt, double *yt, double *zt)
{


//!
//!     purpose: interpolation of a function of three variables in an
//!              equidistant(!!!) table.
//!
//!     method:  Catmull-Rom cubic interpolation formula
//!
//!     x        input vector of first  variable
//!     y        input vector of second variable
//!     z        input vector of third  variable
//!
//!     f        output vector of interpolated function values
//!
//!     ft       3d array of tabulated function values
//!     nx       x-dimension of table
//!     ny       y-dimension of table
//!     nz       z-dimension of table
//!     xt       vector of x-coordinates of table
//!     yt       vector of y-coordinates of table
//!     zt       vector of z-coordinates of table
//!
//!---------------------------------------------------------------------


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
   dx = (xt[nx-1] - xt[0]) / (1.0*(nx-1));
   dy = (yt[ny-1] - yt[0]) / (1.0*(ny-1));
   dz = (zt[nz-1] - zt[0]) / (1.0*(nz-1));

   dxi = 1.0/dx;
   dyi = 1.0/dy;
   dzi = 1.0/dz;
#endif

#if 0
   dx = drho;
   dy = dtemp;
   dz = dye;
   dxi = drhoi;
   dyi = dtempi;
   dzi = dyei;
#endif


   // determine location in table
   ix = (int)( (x - xt[0] - 1.0e-10) * dxi );
   iy = (int)( (y - yt[0] - 1.0e-10) * dyi );
   iz = (int)( (z - zt[0] - 1.0e-10) * dzi );

   ix = MAX( 1, MIN( ix, nx-2 ) );
   iy = MAX( 1, MIN( iy, ny-2 ) );
   iz = MAX( 1, MIN( iz, nz-2 ) );

   
   // linear interpolation at boundaries
   if (ix == 1 || iy == 1 || iz == 1 || 
       ix == nx-2 || iy == ny-2 || iz == nz-2)
   {
      nuc_eos_C_linterp_some(x, y, z, f, ft, ivs,
                             nx, ny, nz, nvars, xt, yt, zt);
      return;
   }


   // difference 
   delx = (x - xt[ix])*dxi;
   dely = (y - yt[iy])*dyi;
   delz = (z - zt[iz])*dzi;

   // factors for Catmull-Rom interpolation
   u[0] = -0.5*CUBE(delx) + SQR(delx) - 0.5*delx;
   u[1] = 1.5*CUBE(delx) - 2.5*SQR(delx) + 1.0;
   u[2] = -1.5*CUBE(delx) + 2.0*SQR(delx) + 0.5*delx;
   u[3] = 0.5*CUBE(delx) - 0.5*SQR(delx);

   v[0] = -0.5*CUBE(dely) + SQR(dely) - 0.5*dely;
   v[1] = 1.5*CUBE(dely) - 2.5*SQR(dely) + 1.0;
   v[2] = -1.5*CUBE(dely) + 2*SQR(dely) + 0.5*dely;
   v[3] = 0.5*CUBE(dely) - 0.5*SQR(dely);

   w[0] = -0.5*CUBE(delz) + SQR(delz) - 0.5*delz;
   w[1] = 1.5*CUBE(delz) - 2.5*SQR(delz) + 1.0;
   w[2] = -1.5*CUBE(delz) + 2.0*SQR(delz) + 0.5*delz;
   w[3] = 0.5*CUBE(delz) - 0.5*SQR(delz);
   
   for (int iv = 0; iv < nvars; iv++)
   {
      vox = 0.0;
      pv = ft + iv + NTABLES*((ix-1) + (iy-1)*nx + (iz-1)*nxy);
      for (int k = 0; k < 4; k++)
      {
         q[k] = 0.0;
         for (int j = 0; j < 4; j++)
         {
            r[j] = 0.0;
            for (int i = 0; i < 4; i++)
            {
               r[j] += u[i]* *pv;
               pv += NTABLES;
            }
            q[k] += v[j]*r[j];
            pv += NTABLES*nx - 4*NTABLES;
         }
         vox += w[k]*q[k];
         pv += nxy*NTABLES - 4*NTABLES*nx;
      }
      f[iv] = vox;      
   }

   if (isnan(*f)) 
   {
      nuc_eos_C_linterp_some(x, y, z, f, ft, ivs,
                             nx, ny, nz, nvars, xt, yt, zt);
      return;
   }


   return;
}

