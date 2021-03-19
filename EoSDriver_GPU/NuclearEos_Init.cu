#include "NuclearEos.h"
#include "NuclearEos.cu"
#include "stdlib.h"
#include "stdio.h"


__global__ static
void test_kernel(double *d_EoS_Table[EOS_NTABLE_MAX], int *d_nye);

__global__ static
void NuclearEoS_kernel(int *d_nrho, int *d_neps, int *d_nye,
                       int *d_nmode, double *d_energy_shift,
                       double *d_EoS_Table[EOS_NTABLE_MAX],
                       double *d_Rand_Vars[NUC_TABLE_NPTR], int *d_err_check);
 

void NuclearEoS_Init() {


   long Err_Check_Size = sizeof(int) * npoints;

   int *h_err_check = NULL;
   int *d_err_check = NULL;
   h_err_check = (int*)malloc( Err_Check_Size );
   cudaMalloc( (void**)&d_err_check, Err_Check_Size );

   // call nuclear EoS kernel
   //NuclearEoS_kernel<<<10, 10>>>( d_nrho, d_neps, d_nye, d_nmode, d_energy_shift, 
   //                               &d_EoS_Table[EOS_NTABLE_MAX], 
   //                               &d_Rand_Vars[NUC_TABLE_NPTR], 
   //                               d_err_check);                              

   test_kernel<<<1, 10>>>(&d_EoS_Table[EOS_NTABLE_MAX], d_nye);


   cudaMemcpy( h_err_check, d_err_check, Err_Check_Size, cudaMemcpyDeviceToHost );
   
   int err_n;
   double err_rate;
   
   for ( int i=0; i<npoints; i++ ) {
      if ( h_err_check[i] != 0 ) 
      {
      //   printf("%d\n", h_err_check[i]);
         err_n++;
      }
   }
   
   err_rate = double(1.0*err_n)/double(1.0*npoints); 
   //printf("Error rate for prss mode: %6.4f\n", err_rate);
   cudaFree( d_err_check );
   d_err_check = NULL;
   free( h_err_check );
   h_err_check = NULL;


}

__global__
void test_kernel(double *d_EoS_Table[EOS_NTABLE_MAX], int *d_nye) {

   //printf("Hello\n");
   const int ye_size = *d_nye;
   //printf("%d\n", ye_size);
   double *ye;
   ye = (double*)malloc(sizeof(double)* *d_nye);
   ye = d_EoS_Table[NUC_TAB_YE];
   printf("%5.2f\n", ye[0]);

}

__global__
void NuclearEoS_kernel(int *d_nrho, int *d_neps, int *d_nye,
                       int *d_nmode, double *d_energy_shift,
                       double *d_EoS_Table[EOS_NTABLE_MAX],
                       double *d_Rand_Vars[NUC_TABLE_NPTR], 
                       int *d_err_check) {
   
   long t = blockDim.x*blockIdx.x + threadIdx.x;

   if (t >= npoints) return;

   int keymode = NUC_MODE_PRES;
   int keyerr  = 0;
   
   //printf("%d\n", keyerr);
   //printf("%d\n", *d_nrho);
//   printf("%5.2E\n", d_Rand_Vars[0][t]);
   //printf("%5.2E\n", *d_EoS_Table[NUC_TAB_YE][0]);

   //nuc_eos_C_short( d_Rand_Vars[0 + NUC_TABLE_NPTR * t], &d_Rand_Vars[1 + NUC_TABLE_NPTR * t], 
   //                 d_Rand_Vars[2 + NUC_TABLE_NPTR * t], &d_Rand_Vars[3 + NUC_TABLE_NPTR * t], 
   //                &d_Rand_Vars[4 + NUC_TABLE_NPTR * t], &d_Rand_Vars[5 + NUC_TABLE_NPTR * t], 
   //                &d_Rand_Vars[6 + NUC_TABLE_NPTR * t], &d_Rand_Vars[7 + NUC_TABLE_NPTR * t],
   //                *d_energy_shift, *d_nrho, *d_neps, *d_nye, 
   //                *d_nmode, d_EoS_Table[NUC_TAB_ALL], 
   //                d_EoS_Table[NUC_TAB_ALL_MODE],  d_EoS_Table[NUC_TAB_RHO], 
   //                d_EoS_Table[NUC_TAB_EPS],       d_EoS_Table[NUC_TAB_YE], 
   //                d_EoS_Table[NUC_TAB_TEMP_MODE], d_EoS_Table[NUC_TAB_ENTR_MODE], 
   //                d_EoS_Table[NUC_TAB_PRES_MODE], keymode, &keyerr );

   //nuc_eos_C_short( xrho, &xenr, xye, &xtemp, &xent, &xprs, &xcs2, &xmunu,
   //                 *d_energy_shift, *d_nrho, *d_neps, *d_nye, 
   //                 *d_nmode, d_EoS_Table[NUC_TAB_ALL], 
   //                 d_EoS_Table[NUC_TAB_ALL_MODE],  d_EoS_Table[NUC_TAB_RHO], 
   //                 d_EoS_Table[NUC_TAB_EPS],       d_EoS_Table[NUC_TAB_YE], 
   //                 d_EoS_Table[NUC_TAB_TEMP_MODE], d_EoS_Table[NUC_TAB_ENTR_MODE], 
   //                 d_EoS_Table[NUC_TAB_PRES_MODE], keymode, &keyerr, rfeps );

   //nuc_eos_C_short( d_Rand_Vars[0][t], &d_Rand_Vars[1][t], 
   //                 d_Rand_Vars[2][t], &d_Rand_Vars[3][t], 
   //                &d_Rand_Vars[4][t], &d_Rand_Vars[5][t], 
   //                &d_Rand_Vars[6][t], &d_Rand_Vars[7][t], 
   //                *d_energy_shift, *d_nrho, *d_neps, *d_nye, 
   //                *d_nmode, d_EoS_Table[NUC_TAB_ALL], 
   //                d_EoS_Table[NUC_TAB_ALL_MODE],  d_EoS_Table[NUC_TAB_RHO], 
   //                d_EoS_Table[NUC_TAB_EPS],       d_EoS_Table[NUC_TAB_YE], 
   //                d_EoS_Table[NUC_TAB_TEMP_MODE], d_EoS_Table[NUC_TAB_ENTR_MODE], 
   //                d_EoS_Table[NUC_TAB_PRES_MODE], keymode, &keyerr );


   //printf("%d\n", keyerr);
   //d_err_check[t] = keyerr;

}