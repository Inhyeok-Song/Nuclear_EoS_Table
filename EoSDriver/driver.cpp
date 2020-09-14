#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "NuclearEos.h"

int nrho;
int neps;
int nye;
int nmode;

double *alltables;
double *alltables_mode;
double *logrho;
double *logeps;
double *yes;
double energy_shift;

double *logtemp_mode;
double *entr_mode;
double *logprss_mode;


int main(void) 
{
   // Read the tabulate EoS
   nuc_eos_C_ReadTable( "../LS220_eps.h5" );
  
   
   printf("nrho: %d\n", nrho);
   printf("neps: %d\n", neps);
   printf("nye: %d\n",  nye);
 
   // testing
   nuc_eos_C_testing();
 
 
   return 0;
}
