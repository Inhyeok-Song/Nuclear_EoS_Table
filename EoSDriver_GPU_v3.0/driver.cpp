#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "NuclearEos.h"

//CPU variables
int    g_nrho;
int    g_neps;
int    g_nye;
int    g_nmode;
double g_energy_shift;

double *g_alltables;
double *g_alltables_mode;
double *g_logrho;
double *g_logeps;
double *g_yes;
double *g_logtemp_mode;
double *g_logprss_mode;
double *g_entr_mode;

double EoS_AuxArray_Flt[EOS_NAUX_MAX];
double EoS_AuxArray_Int[EOS_NAUX_MAX];

double *h_EoS_Table[EOS_NTABLE_MAX];

//GPU variables
int    *d_nrho         = NULL;
int    *d_neps         = NULL;
int    *d_nye          = NULL;
int    *d_nmode        = NULL;
double *d_energy_shift = NULL;

double *d_EoS_Table[EOS_NTABLE_MAX] = {NULL};

double *d_Rand_Vars = NULL;


void PassNuclearEoSTable2GPU(void);
void PassRandPoints2GPU(void);
void NuclearEoS_Init(void);
void MemFree_NuclearEoS(void);


int main( int argc, char* argv[] ) 
{
   // Read the tabulate EoS
   nuc_eos_C_ReadTable( "../LS220_eps.h5" );
  
   printf("nrho: %d\n", g_nrho);
   printf("neps: %d\n", g_neps);
   printf("nye: %d\n",  g_nye);

   clock_t time_init, time_end;
   double  excute_time;
   // send CPU table to GPU
   PassNuclearEoSTable2GPU();
 
   // send Rand_Vars to GPU
   time_init = clock();
   PassRandPoints2GPU();
   time_end = clock();
   excute_time = (double)(time_end - time_init) / (CLOCKS_PER_SEC);

   // test NuclearEos in GPU
   NuclearEoS_Init();

   // excuete time
   //printf("%5.2e\n", excute_time);

   // free GPU memory for EoS table
   MemFree_NuclearEoS();
 
   return 0;
}
