#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NuclearEos.h"
#include "NuclearEos.cuh"

double *Rand_Vars;


void Pick_Rand_Points();


void PassRandPoints2GPU() {

   Pick_Rand_Points();


   unsigned long Rand_Size = sizeof(double)*NUC_TABLE_NPTR*npoints;
      

   CUDA_CHECK_ERROR( cudaMalloc( (void**)&d_Rand_Vars, Rand_Size ) );
   CUDA_CHECK_ERROR( cudaMemcpy( d_Rand_Vars, Rand_Vars, Rand_Size, cudaMemcpyHostToDevice ) );

   free(Rand_Vars);

} //FUNCTION : PassRandPoints2GPU

void Pick_Rand_Points()
{
   
   // set min and max values
   double lrhomin    = log10(pow(10.0, g_logrho[0])*1.02);
   double lrhomax    = log10(pow(10.0, g_logrho[g_nrho-1])*0.98);
   double lepsmin    = log10(1.02*1E17); //log10(pow(10.0,logeps[0])*1.02);
   double lepsmax    = log10(0.98*1E21); //log10(pow(10.0,logeps[neps-1])*0.98);
   double yemin      = g_yes[0]*1.02;
   double yemax      = g_yes[g_nye-1]*0.98;
   double lr_range   = lrhomax-lrhomin;
   double leps_range = lepsmax-lepsmin;
   double ye_range   = yemax-yemin;

   double lprssmax    = log10(pow(10.0, 0.98 * 18.4581859998271) );
   double lprssmin    = log10(pow(10.0, 1.02 * 37.84367789977789));
   double lprss_range = lprssmax - lprssmin;

   unsigned long Rand_Size = sizeof(double)*NUC_TABLE_NPTR*npoints;
   Rand_Vars = (double*)malloc( Rand_Size );

   for (int i=0; i<npoints; i++) {
      // set random points (rho, eps, Ye) 
      double rand1 = rand() / (double)RAND_MAX; 
      double rand2 = rand() / (double)RAND_MAX;
      double rand3 = rand() / (double)RAND_MAX;
      double rand4 = rand() / (double)RAND_MAX;

      double rand_lr    = lrhomin  + lr_range    * rand1;
      double rand_leps  = lepsmin  + leps_range  * rand2;
      double rand_ye    = yemin    + ye_range    * rand3;
   
      double rand_lprss = lprssmin + lprss_range * rand4;

      double array[NUC_TABLE_NPTR];
      array[0] = pow(10.0, rand_lr);
      array[1] = pow(10.0, rand_leps);                      
      array[2] = rand_ye;
      array[3] = 0.1; 
      array[4] = 0.1;
      array[5] = pow(10.0, rand_lprss);
      array[6] = 0.1;
      array[7] = 0.1;
      
      for (int j=0; j<NUC_TABLE_NPTR; j++) {
         Rand_Vars[j + NUC_TABLE_NPTR*i] = array[j];
      }

   }
 

} // FUNCTION : Pick_Rand_Points