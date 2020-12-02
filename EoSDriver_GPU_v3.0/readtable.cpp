#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define H5_USE_16_API 1
#include "hdf5.h"
//#include "GAMER.h"
#include "NuclearEos.h"

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                \
   do {                                                    \
      int _error_code = fn_call;                           \
      if (_error_code < 0)                                 \
      {                                                    \
         fprintf(stderr,                                   \
                 "HDF5 call '%s' returned error code %d",  \
                  #fn_call, _error_code);                  \
         abort();                                          \
      }                                                    \
   } while (0)

static int file_is_readable( const char *filename );
static int file_is_readable( const char *filename )
{
   FILE* fp = NULL;
      fp = fopen(filename, "r");
      if(fp != NULL)
      {
         fclose(fp);
         return 1;
      }
   return 0;
}

// Cactus calls this function. It reads in the table and calls a fortran
// function to setup values for the fortran eos module
void nuc_eos_C_ReadTable( char *nuceos_table_name )
{

   fprintf(stdout,"*******************************\n");
   fprintf(stdout,"Reading nuc_eos_eps table file:\n");
   fprintf(stdout,"%s\n", nuceos_table_name);
   fprintf(stdout,"*******************************\n");

   hid_t file;
   if (!file_is_readable(nuceos_table_name)) 
   {
      fprintf(stderr,
              "Could not read nuceos_table_name %s \n",
              nuceos_table_name);
      abort();
   }
   HDF5_ERROR( file = H5Fopen( nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT ) );


// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_EOS_HDF5(NAME,VAR,TYPE,MEM)                                             \
   do {                                                                              \
      hid_t dataset;                                                                 \
      HDF5_ERROR( dataset = H5Dopen( file, NAME ) );                                 \
      HDF5_ERROR( H5Dread( dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR ) );        \
      HDF5_ERROR( H5Dclose( dataset ) );                                             \
   } while (0)
#define READ_EOSTABLE_HDF5( NAME, OFF )                                              \
   do {                                                                              \
      hsize_t offset[2]     = { OFF, 0 };                                            \
      H5Sselect_hyperslab( mem3, H5S_SELECT_SET, offset, NULL, var3, NULL );         \
      READ_EOS_HDF5( NAME, alltables_temp, H5T_NATIVE_DOUBLE, mem3 );                \
   } while (0)
#define READ_EOSTABLE_MODE_HDF5( NAME,OFF )                                          \
   do {                                                                              \
      hsize_t offset[2]     = { OFF,0 };                                             \
      H5Sselect_hyperslab(mem3_mode, H5S_SELECT_SET, offset, NULL, var3_mode, NULL );\
      READ_EOS_HDF5( NAME, alltables_mode_temp, H5T_NATIVE_DOUBLE, mem3_mode );      \
   } while (0)

  
   // Read size of tables
   READ_EOS_HDF5("pointsrho",    &g_nrho,  H5T_NATIVE_INT, H5S_ALL);
   READ_EOS_HDF5("pointsenergy", &g_neps,  H5T_NATIVE_INT, H5S_ALL);
   READ_EOS_HDF5("pointsye",     &g_nye,   H5T_NATIVE_INT, H5S_ALL);
   READ_EOS_HDF5("points_mode",  &g_nmode, H5T_NATIVE_INT, H5S_ALL);


   // Allocate memory for tables
   double *alltables_temp;
   if ( !(alltables_temp = (double*)malloc(g_nrho * g_neps * g_nye * NUC_TABLE_NVAR * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   double *alltables_mode_temp;
   if ( !(alltables_mode_temp = (double*)malloc(g_nrho * g_nmode * g_nye * 3 * sizeof(double))) )
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   if ( !(g_logrho = (double*)malloc(g_nrho * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   if ( !(g_logeps = (double*)malloc(g_neps * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   if ( !(g_yes = (double*)malloc(g_nye * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   if ( !(g_logtemp_mode = (double*)malloc(g_nmode * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   if ( !(g_entr_mode = (double*)malloc(g_nmode * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   if ( !(g_logprss_mode = (double*)malloc(g_nmode * sizeof(double))) )
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   
   
   // Prepare HDF5 to read hyperslabs into alltables_temp
   hsize_t table_dims[2] = {NUC_TABLE_NVAR, g_nrho * g_neps * g_nye};
   hsize_t var3[2]       = {1, g_nrho * g_neps * g_nye};
   hid_t   mem3          = H5Screate_simple(2, table_dims, NULL);
   
   hsize_t table_dims_mode[2] = {3, g_nrho * g_nmode * g_nye};
   hsize_t var3_mode[2]       = {1, g_nrho * g_nmode * g_nye};
   hid_t   mem3_mode          = H5Screate_simple(2, table_dims_mode, NULL);
   
   // Read alltables_temp
   READ_EOSTABLE_HDF5("logpress",  0);
   READ_EOSTABLE_HDF5("logtemp",   1);
   READ_EOSTABLE_HDF5("entropy",   2);
   READ_EOSTABLE_HDF5("munu",      3);
   READ_EOSTABLE_HDF5("cs2",       4);
   // chemical potentials
   READ_EOSTABLE_HDF5("muhat",     5);
   READ_EOSTABLE_HDF5("mu_e",      6);
   READ_EOSTABLE_HDF5("mu_p",      7);
   READ_EOSTABLE_HDF5("mu_n",      8);
   // compositions
   READ_EOSTABLE_HDF5("Xa",        9);
   READ_EOSTABLE_HDF5("Xh",       10);
   READ_EOSTABLE_HDF5("Xn",       11);
   READ_EOSTABLE_HDF5("Xp",       12);
   // average nucleus
   READ_EOSTABLE_HDF5("Abar",     13);
   READ_EOSTABLE_HDF5("Zbar",     14);
   // Gamma
   READ_EOSTABLE_HDF5("gamma",    15);
   
   // energy for temp, entr modes
   READ_EOSTABLE_MODE_HDF5("logenergy_temp", 0);
   READ_EOSTABLE_MODE_HDF5("logenergy_entr", 1);
   READ_EOSTABLE_MODE_HDF5("logenergy_prss", 2);
   
   // Read additional tables and variables
   READ_EOS_HDF5("logrho",       g_logrho,        H5T_NATIVE_DOUBLE, H5S_ALL);
   READ_EOS_HDF5("logenergy",    g_logeps,        H5T_NATIVE_DOUBLE, H5S_ALL);
   READ_EOS_HDF5("ye",           g_yes,           H5T_NATIVE_DOUBLE, H5S_ALL);
   READ_EOS_HDF5("logtemp_mode", g_logtemp_mode,  H5T_NATIVE_DOUBLE, H5S_ALL);
   READ_EOS_HDF5("entr_mode",    g_entr_mode,     H5T_NATIVE_DOUBLE, H5S_ALL);
   READ_EOS_HDF5("logprss_mode", g_logprss_mode,  H5T_NATIVE_DOUBLE, H5S_ALL);
   READ_EOS_HDF5("energy_shift", &g_energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL);
   
   HDF5_ERROR(H5Sclose(mem3));
   HDF5_ERROR(H5Sclose(mem3_mode));
   HDF5_ERROR(H5Fclose(file));
   
   
   // Change ordering of alltables array so that
   // the table kind is the fastest changing index
   if ( !(g_alltables = (double*)malloc(g_nrho * g_neps * g_nye * NUC_TABLE_NVAR * sizeof(double))) ) 
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   
   }
   for(int iv=0; iv<NUC_TABLE_NVAR; iv++) 
      for(int k=0; k<g_nye; k++) 
         for(int j=0; j<g_neps; j++) 
   for(int i=0; i<g_nrho ; i++) 
   {
      int indold = i + g_nrho*(j + g_neps*(k + g_nye*iv));
      int indnew = iv + NUC_TABLE_NVAR*(i + g_nrho*(j + g_neps*k));
      g_alltables[indnew] = alltables_temp[indold];
   }
   
   if ( !(g_alltables_mode = (double*)malloc(g_nrho * g_nmode * g_nye * 3 * sizeof(double))) )
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   for(int iv=0; iv<3; iv++) 
      for(int k=0; k<g_nye; k++) 
         for(int j=0; j<g_nmode; j++)
   for(int i = 0; i < g_nrho ; i++)
   {
      int indold = i + g_nrho*(j + g_nmode*(k + g_nye*iv));
      int indnew = iv + 3*(i + g_nrho*(j + g_nmode*k));
      g_alltables_mode[indnew] = alltables_mode_temp[indold];
   }
   
   
   // Free memory of temporary array
   free(alltables_temp);
   free(alltables_mode_temp);

   // set the EoS table pointers
   h_EoS_Table[NUC_TAB_ALL      ] = g_alltables;
   h_EoS_Table[NUC_TAB_ALL_MODE ] = g_alltables_mode;
   h_EoS_Table[NUC_TAB_RHO      ] = g_logrho;
   h_EoS_Table[NUC_TAB_EPS      ] = g_logeps;
   h_EoS_Table[NUC_TAB_YE       ] = g_yes;
   h_EoS_Table[NUC_TAB_TEMP_MODE] = g_logtemp_mode;
   h_EoS_Table[NUC_TAB_ENTR_MODE] = g_entr_mode;
   h_EoS_Table[NUC_TAB_PRES_MODE] = g_logprss_mode;


}
