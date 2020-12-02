#include "NuclearEos.h"
#include "NuclearEos.cuh"
#include <stdio.h>

void PassNuclearEoSTable2GPU() {

   long EoS_TableSize[NUC_TABLE_NPTR];

   EoS_TableSize[NUC_TAB_ALL      ] = sizeof(double)*g_nrho*g_neps*g_nye*NUC_TABLE_NVAR;
   EoS_TableSize[NUC_TAB_ALL_MODE ] = sizeof(double)*g_nrho*g_nmode*g_nye*3;
   EoS_TableSize[NUC_TAB_RHO      ] = sizeof(double)*g_nrho;
   EoS_TableSize[NUC_TAB_EPS      ] = sizeof(double)*g_neps;
   EoS_TableSize[NUC_TAB_YE       ] = sizeof(double)*g_nye;
   EoS_TableSize[NUC_TAB_TEMP_MODE] = sizeof(double)*g_nmode;
   EoS_TableSize[NUC_TAB_ENTR_MODE] = sizeof(double)*g_nmode;
   EoS_TableSize[NUC_TAB_PRES_MODE] = sizeof(double)*g_nmode;


   for( int t=0; t<NUC_TABLE_NPTR; t++ )
   {
      CUDA_CHECK_ERROR( cudaMalloc( (void**) &d_EoS_Table[t], EoS_TableSize[t] ) );
      CUDA_CHECK_ERROR( cudaMemcpy( d_EoS_Table[t], h_EoS_Table[t], EoS_TableSize[t], cudaMemcpyHostToDevice ) );
   }
 
   CUDA_CHECK_ERROR( cudaMalloc( (void**)&d_energy_shift, sizeof(double) ) );
   CUDA_CHECK_ERROR( cudaMemcpy( d_energy_shift, &g_energy_shift, sizeof(double), cudaMemcpyHostToDevice ) );

   CUDA_CHECK_ERROR( cudaMalloc( (void**)&d_nrho, sizeof(int) ) );
   CUDA_CHECK_ERROR( cudaMemcpy( d_nrho, &g_nrho, sizeof(int), cudaMemcpyHostToDevice ) );
   
   CUDA_CHECK_ERROR( cudaMalloc( (void**)&d_neps, sizeof(int) ) );
   CUDA_CHECK_ERROR( cudaMemcpy( d_neps, &g_neps, sizeof(int), cudaMemcpyHostToDevice ) );
   
   CUDA_CHECK_ERROR( cudaMalloc( (void**)&d_nye, sizeof(int) ) );
   CUDA_CHECK_ERROR( cudaMemcpy( d_nye, &g_nye, sizeof(int), cudaMemcpyHostToDevice ) );
   
   CUDA_CHECK_ERROR( cudaMalloc( (void**)d_nmode, sizeof(int) ) );
   CUDA_CHECK_ERROR( cudaMemcpy( d_nmode, &g_nmode, sizeof(int), cudaMemcpyHostToDevice ) );
}
