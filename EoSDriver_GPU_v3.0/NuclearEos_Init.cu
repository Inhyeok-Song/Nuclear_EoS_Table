#include "NuclearEos.h"
#include "NuclearEos.cu"
#include "stdlib.h"
#include "stdio.h"


__global__ static
void NuclearEoS_kernel(int *d_nrho, int *d_neps, int *d_nye,
                       int *d_nmode, double *d_energy_shift,
                       double *d_EoS_Table[EOS_NTABLE_MAX],
                       int *d_keymode, int *d_keyerr, double *d_rfeps,
                       double *d_Rand_Vars, int *d_err_check);
 

void NuclearEoS_Init() {


   int h_keymode, h_keyerr;
   h_keymode = NUC_MODE_PRES;
   h_keyerr  = 0;
   double h_rfeps = 1.E-10;

   int    *d_keymode = NULL;
   int    *d_keyerr  = NULL;
   double *d_rfeps   = NULL;
   unsigned long Err_Check_Size = sizeof(int) * npoints;

   cudaMalloc( (void**)&d_keymode, sizeof(int) );
   cudaMalloc( (void**)&d_keyerr,  sizeof(int) );
   cudaMalloc( (void**)&d_rfeps,   sizeof(double) );
   cudaMemcpy( d_keymode, &h_keymode, sizeof(int),    cudaMemcpyHostToDevice);
   cudaMemcpy( d_keyerr,  &h_keyerr,  sizeof(int),    cudaMemcpyHostToDevice);
   cudaMemcpy( d_rfeps,   &h_rfeps,   sizeof(double), cudaMemcpyHostToDevice);


   int *h_err_check = NULL;
   int *d_err_check = NULL;
   h_err_check = (int*)malloc( Err_Check_Size );
   cudaMalloc( (void**)&d_err_check, Err_Check_Size );


   NuclearEoS_kernel<<<16, 16>>>(d_nrho, d_neps, d_nye, d_nmode, d_energy_shift, 
                                 &d_EoS_Table[EOS_NTABLE_MAX], d_keymode, d_keyerr, 
                                 d_rfeps, d_Rand_Vars, d_err_check);                              
   

   cudaMemcpy( h_err_check, d_err_check, Err_Check_Size, cudaMemcpyDeviceToHost );
   int err_n;
   double err_rate;
   
   for (int i=0; i<npoints; i++) {
      //if (h_err_check[i] != 0) {
         //printf("%d\n", h_err_check[i]);
      //   err_n++;
      //}
   }
    err_rate = double(1.0*err_n)/double(1.0*npoints);
    printf("Error rate for prss mode: %6.4f\n", err_rate);

   //cudaFree(d_keymode);
   //cudaFree(d_keyerr);
   //cudaFree(d_rfeps);
//
   //free(h_err_check);
   //cudaFree(d_err_check);

}

__global__
void NuclearEoS_kernel(int *d_nrho, int *d_neps, int *d_nye,
                       int *d_nmode, double *d_energy_shift,
                       double *d_EoS_Table[EOS_NTABLE_MAX],
                       int *d_keymode, int *d_keyerr, double *d_rfeps,
                       double *d_Rand_Vars, int *d_err_check) {
   
   int t = blockDim.x*blockIdx.x + threadIdx.x;


   //printf("%5.2E\n", d_Rand_Vars[0 + NUC_TABLE_NPTR * t]);

   nuc_eos_C_short( d_Rand_Vars[0 + NUC_TABLE_NPTR * t], &d_Rand_Vars[1 + NUC_TABLE_NPTR * t], 
                    d_Rand_Vars[2 + NUC_TABLE_NPTR * t], &d_Rand_Vars[3 + NUC_TABLE_NPTR * t], 
                   &d_Rand_Vars[4 + NUC_TABLE_NPTR * t], &d_Rand_Vars[5 + NUC_TABLE_NPTR * t], 
                   &d_Rand_Vars[6 + NUC_TABLE_NPTR * t], &d_Rand_Vars[7 + NUC_TABLE_NPTR * t],
                   *d_energy_shift, *d_nrho, *d_neps, *d_nye, 
                   *d_nmode, d_EoS_Table[NUC_TAB_ALL], 
                   d_EoS_Table[NUC_TAB_ALL_MODE],  d_EoS_Table[NUC_TAB_RHO], 
                   d_EoS_Table[NUC_TAB_EPS],       d_EoS_Table[NUC_TAB_YE], 
                   d_EoS_Table[NUC_TAB_TEMP_MODE], d_EoS_Table[NUC_TAB_ENTR_MODE], 
                   d_EoS_Table[NUC_TAB_PRES_MODE], *d_keymode, d_keyerr, *d_rfeps );

   //nuc_eos_C_short( xrho, &xenr, xye, &xtemp, &xent, &xprs, &xcs2, &xmunu,
   //                 *d_energy_shift, *d_nrho, *d_neps, *d_nye, 
   //                 *d_nmode, d_EoS_Table[NUC_TAB_ALL], 
   //                 d_EoS_Table[NUC_TAB_ALL_MODE],  d_EoS_Table[NUC_TAB_RHO], 
   //                 d_EoS_Table[NUC_TAB_EPS],       d_EoS_Table[NUC_TAB_YE], 
   //                 d_EoS_Table[NUC_TAB_TEMP_MODE], d_EoS_Table[NUC_TAB_ENTR_MODE], 
   //                 d_EoS_Table[NUC_TAB_PRES_MODE], keymode, &keyerr, rfeps );
   
   //printf("\n");
   d_err_check[t] = *d_keyerr;

}