#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NuclearEos.h"


//-------------------------------------------------------
// Function    :  nuc_eos_C_testing
// Description :  testing Nuclear equation of state
//                including each mode
//-------------------------------------------------------
void nuc_eos_C_testing()
{
   int npoints          = 50000;
   int npoints_valid[4];
   for (int i=0; i<4; i++) npoints_valid[i] = npoints;

   // Randomly draw npoints from
   // range of logrho, logeps, Y_e 
   double lrhomin    = log10(pow(10.0, logrho[0])*1.02);
   double lrhomax    = log10(pow(10.0, logrho[nrho-1])*0.98);
   double lepsmin    = log10(1.02*1E17); //log10(pow(10.0,logeps[0])*1.02);
   double lepsmax    = log10(0.98*1E21); //log10(pow(10.0,logeps[neps-1])*0.98);
   double yemin      = yes[0]*1.02;
   double yemax      = yes[nye-1]*0.98;
   double lr_range   = lrhomax-lrhomin;
   double leps_range = lepsmax-lepsmin;
   double ye_range   = yemax-yemin;
   
   // output files for writing errors
   FILE *output_error[6];

   FILE *ener_err1 = fopen("ener_err1.txt", "w"); 
   fprintf(ener_err1, "#    Num        logrho          logenergy             Ye\n", "w");

   FILE *temp_err1 = fopen("temp_err1.txt", "w"); 
   FILE *temp_err2 = fopen("temp_err2.txt", "w");

   FILE *entr_err1 = fopen("entr_err1.txt", "w"); 
   FILE *entr_err2 = fopen("entr_err2.txt", "w");

   FILE *prss_err1 = fopen("prss_err1.txt", "w"); 
   FILE *prss_err2 = fopen("prss_err2.txt", "w");
   
   output_error[0] = temp_err1;
   output_error[1] = entr_err1;
   output_error[2] = prss_err1;
   output_error[3] = temp_err2;
   output_error[4] = entr_err2;
   output_error[5] = prss_err2;
   
   for ( int i=0; i<6; i++ )
   {
      if ( i < 3 )
      {
         fprintf(output_error[i], "#    Num        logrho          logenergy             Ye\n", "w");
      }
      else if ( i >= 3 )
      {
         fprintf(output_error[i], "#    Num        logrho          logenergy             Ye            rel_error\n", "w");
      }
   }

   // variables for writing errors
   int failcount1[4] = {0};   // count error1: fail in finding energy
   int failcount2[4] = {0};   // count error2: relative error > 10^-3
   double max_err[4] = {0.0};
      
   // tolerence
   const double rfeps = 1.0E-10;
   
   for (int i=0; i<npoints; i++)
   {
   
      // test energy mode
      double xrho, xenr, xye, xtemp, xent, xprs, xcs2, xmunu;
      int    keymode = 0;
      int    keyerr;
      
      // set random points (rho, eps, Ye) 
      double rand1 = rand() / (double)RAND_MAX; 
      double rand2 = rand() / (double)RAND_MAX;
      double rand3 = rand() / (double)RAND_MAX;
   
      double rand_lr   = lrhomin + lr_range   * rand1;
      double rand_leps = lepsmin + leps_range * rand2;
      double rand_ye   = yemin   + ye_range   * rand3;
     
      xrho = pow(10.0, rand_lr);
      xenr = pow(10.0, rand_leps);
      xye  = rand_ye;
      
      double xlr       = rand_lr;
      double savedleps = rand_leps; // save the original log(eps)
      
      // NuclearEos function
      nuc_eos_C_short( xrho, &xenr, xye, &xtemp, &xent, &xprs, 
                       &xcs2, &xmunu, energy_shift,
                       nrho, neps, nye, nmode,
                       alltables, alltables_mode,
                       logrho, logeps, yes,
                       logtemp_mode, entr_mode, logprss_mode,
                       keymode, &keyerr, rfeps );
      double res[4];
      res[0] = xenr;
      res[1] = xtemp;
      res[2] = xent;
      res[3] = xprs;

      if ( keyerr != 0 )
      {
         // write the erorr when energy mode failes
         failcount1[0]++;
         npoints_valid[0]--;
         fprintf(ener_err1, "%8d    %10.8e    %10.8e    %10.8e    %3d\n", 
                 failcount1[0], rand_lr, rand_leps, rand_ye, keyerr);
      }
      else if ( keyerr == 0 )
      {
         for ( int i_mode=1; i_mode<4; i_mode++ ) // try each mode
         {
            if ( isnan(res[i_mode]) ) npoints_valid[i_mode]--; // invalid region
                                                               // (no temp, entr or prss)
            else
            {
               keymode = i_mode; // [1=temp, 2=entr, 3=prss] 
               keyerr  = 0;
               xenr    = pow(10.0, savedleps);
               double leps;
               double rel_err = 0.0;  
               // NuclearEos function
               nuc_eos_C_short( xrho, &xenr, xye, &xtemp, &xent, &xprs, 
                                &xcs2, &xmunu, energy_shift,
                                nrho, neps, nye, nmode,
                                alltables, alltables_mode,
                                logrho, logeps, yes,
                                logtemp_mode, entr_mode, logprss_mode,
                                keymode, &keyerr, rfeps );
               leps    = log10(xenr);
               rel_err = fabs((leps-savedleps)/savedleps);

               if ( keyerr != 0 )
               {
                  // error 1: fail in finding energy from the [i_mode]
                  failcount1[i_mode]++;
                  fprintf(output_error[i_mode-1], "%8d    %10.8e    %10.8e    %10.8e\n", 
                          failcount1[i_mode], xlr, savedleps, xye);
               } 
               else if ( rel_err > 1.0e-3 )
               {
                  // error 2: relative errors between the energy mode and [i_mode] > 0.1%
                  failcount2[i_mode]++;
                  fprintf(output_error[i_mode+2], "%8d    %10.8e    %10.8e    %10.8e    %10.8e\n", 
                          failcount2[i_mode], xlr, savedleps, xye, rel_err);
                  max_err[i_mode] = (rel_err > max_err[i_mode]) ? rel_err: max_err[i_mode];
               }
            }
         }  
      }
   }


   // close the output files
   fclose(ener_err1);
   for (int i=0; i<6; i++) fclose(output_error[i]);


   // print results
   fprintf(stdout, "Energy mode\n");
   fprintf(stdout, "failcount1: %d fraction: %15.6E\n", failcount1[0], failcount1[0]/(1.0*npoints));
   fprintf(stdout, "-----------------------------------------------------------------------------\n");

   fprintf(stdout, "Temperature mode\n");
   fprintf(stdout, "npoints_valid: %d\n", npoints_valid[1]);
   fprintf(stdout, "failcount1: %d fraction: %15.6E\n", failcount1[1], failcount1[1]/(1.0*npoints_valid[1]));
   fprintf(stdout, "failcount2: %d fraction: %15.6E\n", failcount2[1], failcount2[1]/(1.0*npoints_valid[1]));
   fprintf(stdout, "Maximum relative error is %f\n", max_err[1]);
   fprintf(stdout, "-----------------------------------------------------------------------------\n");

   fprintf(stdout, "Entropy mode\n");
   fprintf(stdout, "npoints_valid: %d\n", npoints_valid[2]);
   fprintf(stdout, "failcount1: %d fraction: %15.6E\n", failcount1[2], failcount1[2]/(1.0*npoints_valid[2]));
   fprintf(stdout, "failcount2: %d fraction: %15.6E\n", failcount2[2], failcount2[2]/(1.0*npoints_valid[2]));
   fprintf(stdout, "Maximum relative error is %f\n", max_err[2]);
   fprintf(stdout, "-----------------------------------------------------------------------------\n");

   fprintf(stdout, "Pressure mode\n");
   fprintf(stdout, "npoints_valid: %d\n", npoints_valid[3]);
   fprintf(stdout, "failcount1: %d fraction: %15.6E\n", failcount1[3], failcount1[3]/(1.0*npoints_valid[3]));
   fprintf(stdout, "failcount2: %d fraction: %15.6E\n", failcount2[3], failcount2[3]/(1.0*npoints_valid[3]));
   fprintf(stdout, "Maximum relative error is %f\n", max_err[3]);


   return;
} // FUNCTION : nuc_eos_C_testing
