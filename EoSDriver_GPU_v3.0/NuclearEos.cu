#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NuclearEos.h"
#include "findenergy.cu"
#include "cubinterp_some.cu"


//__device__ static
//void nuc_eos_C_short( double xrho, double *xenr, double xye,
//                      double *xtemp, double *xent, double *xprs,
//                      double *xcs2, double *xmunu, double energy_shift,
//                      int nrho, int neps, int nye, int nmode,
//                      double *alltables, double *alltables_mode,
//                      double *logrho, double *logeps, double *yes,
//                      double *logtemp_mode, double *entr_mode, double *logprss_mode,
//                      int keymode, int *keyerr, double rfeps );


//-----------------------------------------------------------------------------------------------
// Function    :  nuc_eos_C_short
// Description :  Function to find thermodynamic varibles by searching
//                a pre-calculated nuclear equation of state table
// 
// Note        :  It will strictly return values in cgs or MeV
//                Four modes are supported
//                The defalut mode is energy (0) mode
//                In case three other modes are available for finding energy
//                temperature (1) mode
//                entropy     (2) mode
//                pressure    (3) mode
//
// Parameter   :  xrho            : input density (rho (g/cm^3))
//                xenr            : specific internal energy (eps)
//                xye             : electron fraction (Y_e)
//                xtemp           : input (temperature mode)
//                                  or ouput temperature in MeV
//                xent            : input (entropy mode)
//                                  or output specific entropy (e)
//                xprs            : input (pressure mode)
//                                  or output pressure
//                xcs2            : output sound speed
//                xmunu           : output chemcial potential
//                energy_shift    : energy_shift
//                nrho            : size of density array in the Nuclear EoS table
//                neps            : size of energy  array in the Nuclear EoS table
//                nye             : size of Y_e     array in the Nuclear EoS table
//                nmode           : size of log(T)     (1)
//                                          entropy    (2)
//                                          log(P)     (3) array in the Nuclear EoS table
//                                                         for each mode
//                alltables       : Nuclear EoS table
//                alltables_mode  : Auxiliary log(eps) arrays for temperature mode
//                                                                entropy mode
//                                                                pressure mode
//                logrho          : log(rho) array in the table
//                logeps          : log(eps) array in the table
//                yes             : Y_e      array in the table
//                logtemp_mode    : log(T)   array for temperature mode
//                entr_mode       : entropy  array for entropy mode
//                logprss_mode    : log(P)   array for pressure mode
//                keymode         : which mode we will use
//                                  0 : energy mode      (coming in with eps)
//                                  1 : temperature mode (coming in with T)
//                                  2 : entropy mode     (coming in with entropy)
//                                  3 : pressure mode    (coming in with P)
//                keyerr          : output error
//                                  667 : fail in finding energy (T, e, P modes)
//                                  101 : Y_e too high
//                                  102 : Y_e too low
//                                  103 : eps too high (if keymode = 0)
//                                  104 : eps too low  (if keymode = 0)
//                                  105 : rho too high
//                                  106 : rho too low
//                rfeps           : tolerence for interpolations
//-----------------------------------------------------------------------------------------------
__device__ inline
void nuc_eos_C_short( double xrho, double *xenr, double xye,
                      double *xtemp, double *xent, double *xprs,
                      double *xcs2, double *xmunu, double energy_shift,
                      int nrho, int neps, int nye, int nmode,
                      double *alltables, double *alltables_mode,
                      double *logrho, double *logeps, double *yes,
                      double *logtemp_mode, double *entr_mode, double *logprss_mode,
                      int keymode, int *keyerr, double rfeps )
{
   printf("sdf");
   *keyerr = 0;
   // check whether (rho, eps, Y_e) is within the table
   if ( log10(xrho) > logrho[nrho-1] )
   {
      *keyerr = 105;
      return;
   }
   if ( log10(xrho) < logrho[0] )
   {
      *keyerr = 106;
      return;
   }
   if ( xye > yes[nye-1] )
   {
      *keyerr = 101;
      return;
   }
   if ( xye < yes[0] )
   {
      *keyerr = 102;
      return;
   }
   
   if ( keymode == 0 )
   {
      if ( log10(*xenr) > logeps[neps-1] )
      {
         *keyerr = 103;
         return;
      }
      if ( log10(*xenr) < logeps[0] )
      {
         *keyerr = 104;
         return;
      }
   }

   
   // set up local vars
   double lr = log10(xrho);
   double xeps = *xenr + energy_shift;
   double leps = log10(MAX(xeps, 1.0));
    
   // find energy if need be
   // ( keymode = 0: energy mode      )
   // ( keymode = 1: temperature mode )
   // ( keymode = 2: entropy mode     )
   // ( keymode = 3  pressure mode    )
   if ( keymode == 1 )
   {
      double lt = log10(*xtemp);
      find_energy( lr, lt, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                   logrho, logtemp_mode, yes, logeps, keymode, keyerr );
      if ( *keyerr != 0 ) return;
      *xenr = pow(10.0, leps) - energy_shift;
   }
   else if ( keymode == 2 )
   {
      double entr = *xent;
      find_energy( lr, entr, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                   logrho, entr_mode, yes, logeps, keymode, keyerr );
      if ( *keyerr != 0 ) return;
      *xenr = pow(10.0, leps) - energy_shift;
   }
   else if ( keymode == 3 )
   {
      double lprs = log10(*xprs);
      find_energy( lr, lprs, xye, &leps, alltables_mode, nrho, nmode, nye, neps,
                   logrho, logprss_mode, yes, logeps, keymode, keyerr );
      if ( *keyerr != 0 ) return;
      *xenr = pow(10.0, leps) - energy_shift;
   }

   double res[5]; // result array

   // linear interolation for other variables
   //nuc_eos_C_linterp_some(lr, leps, xye, res, alltables,
   //                       nrho, neps, nye, 5, logrho, logeps, yes);
   
   // cubic interpolation for other variables
   nuc_eos_C_cubinterp_some( lr, leps, xye, res, alltables,
                             nrho, neps, nye, 5, logrho, logeps, yes );
   

   // assign results
   *xprs  = pow(10.0, res[0]);
   *xtemp = pow(10.0, res[1]);
   *xent  = res[2];
   *xmunu = res[3];
   *xcs2  = res[4];


   return;

} // FUNCTION : nuc_eos_C_short
