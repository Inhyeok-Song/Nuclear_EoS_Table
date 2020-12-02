#include <stdio.h>
#include "NuclearEos.h"


void MemFree_NuclearEoS() {
   

   for (int t=0; t<NUC_TABLE_NPTR; t++)
   {
      if ( d_EoS_Table[t] != NULL ) 
      {
         cudaFree( d_EoS_Table[t] );
         d_EoS_Table[t] = NULL;
      }
   }
   
   cudaFree( d_energy_shift );
   d_energy_shift = NULL;
   cudaFree( d_nrho );
   d_nrho = NULL;
   cudaFree( d_neps );
   d_neps = NULL;
   cudaFree( d_nye );
   d_nye = NULL;
   cudaFree( d_nmode );
   d_nmode = NULL;

   cudaFree( d_Rand_Vars );
   d_Rand_Vars = NULL;

}

