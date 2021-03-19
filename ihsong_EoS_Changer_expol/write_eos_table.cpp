#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "nuc_eos.h"

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                          \
   do                                                                \
   {                                                                 \
      int _error_code = fn_call;                                     \
      if (_error_code < 0)                                           \
      {                                                              \
         fprintf(stderr,                                             \
                 "HDF5 call '%s' returned error code %d",            \
                  #fn_call, _error_code);                            \
         abort();                                                    \
      }                                                              \
   } while (0)


void write_eos_table( char *nuceos_eps_table_name ) 
{


   fprintf(stdout,"*******************************\n");
   fprintf(stdout,"Writing nuc_eos table file:\n");
   fprintf(stdout,"%s\n", nuceos_eps_table_name);
   fprintf(stdout,"*******************************\n");


   hid_t     file;
   hid_t     dataset;
   hid_t     space;
   hsize_t   dims[3]      = {nye, neps, nrho};
   hsize_t   dims_mode[3] = {nye_mode, nmode, nrho_mode};

   
   // Allocating eos table
   double ****eos_table = NULL;  // Write buffer
   if (!(eos_table = (double****)malloc(nye * sizeof(double***))))
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   for (int i=0; i<nye; i++)
   {
      if (!(eos_table[i] = (double***)malloc(neps * sizeof(double**))))
      {
         fprintf(stderr, "Cannot allocate memory for EOS table\n");
         abort();
      }
      for (int j=0; j<neps; j++)
      {
         if (!(eos_table[i][j] = (double**)malloc(nrho * sizeof(double*))))
         {
            fprintf(stderr, "Cannot allocate memory for EOS table\n");
            abort();
         }
         for (int k=0; k<nrho; k++)
         {
            if (!(eos_table[i][j][k] = (double*)malloc((NTABLES) * sizeof(double))))
            {
               fprintf(stderr, "Cannot allocate memory for EOS table\n");
               abort();
            }
         }
      }
   }

   // Allocating eos table for modes (temp, entr, prss)
   double ****eos_table_mode = NULL;   // Write buffer
   if (!(eos_table_mode = (double****)malloc(nye_mode * sizeof(double***))))
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }
   for (int i=0; i<nye_mode; i++)
   {
      if (!(eos_table_mode[i] = (double***)malloc(nmode * sizeof(double**))))
      {
         fprintf(stderr, "Cannot allocate memory for EOS table\n");
         abort();
      }
      for (int j=0; j<nmode; j++)
      {
         if (!(eos_table_mode[i][j] = (double**)malloc(nrho_mode * sizeof(double*))))
         {
            fprintf(stderr, "Cannot allocate memory for EOS table\n");
            abort();
         }
         for (int k=0; k<nrho_mode; k++)
         {
            if (!(eos_table_mode[i][j][k] = (double*)malloc(3 * sizeof(double))))
            {
               fprintf(stderr, "Cannot allocate memory for EOS table\n");
               abort();
            }
         }
      }
   }


   double *wdata;
   if (!(wdata = (double*)malloc(nye * neps * nrho * sizeof(double))))
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }

   double *wdata_mode;
   if (!(wdata_mode = (double*)malloc(nye_mode * nmode * nrho_mode * sizeof(double))))
   {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
   }


   // Find variables. The last row is the elements in the dataspace, i, j 
   // and k are the elements within the array datatype
   for (int i = 0; i<nye; i++)
   for (int j = 0; j<neps; j++)
   for (int k = 0; k<nrho; k++)
   for (int iv= 0; iv<NTABLES; iv++)
   eos_table[i][j][k][iv] = alltables[iv + NTABLES*(k + nrho*(j + neps*i))];

   for (int i = 0; i<nye_mode; i++)
   for (int j = 0; j<nmode; j++)
   for (int k = 0; k<nrho_mode; k++)
   for (int iv= 0; iv<3; iv++)
   eos_table_mode[i][j][k][iv] = alltables_mode[iv + 3*(k + nrho_mode*(j + nmode*i))];

   int N_exp_max = 3;

   for (int N_exp=0; N_exp<N_exp_max; N_exp++) {
      for (int iv=0; iv<NTABLES; iv++)
      {
         for (int k=0; k<nye; k++)
         {
            for (int j=2; j<neps-2; j++)
            {
               for (int i=0; i<nrho; i++)
               {
                  int idx   = iv + NTABLES*(i + nrho*(j + neps*k));
   
                  int idxm3 = iv + NTABLES*(i + nrho*(j-3 + neps*k));
                  int idxm2 = iv + NTABLES*(i + nrho*(j-2 + neps*k));
                  int idxm1 = iv + NTABLES*(i + nrho*(j-1 + neps*k));
                  int idxp1 = iv + NTABLES*(i + nrho*(j+1 + neps*k));
                  int idxp2 = iv + NTABLES*(i + nrho*(j+2 + neps*k));
                  int idxp3 = iv + NTABLES*(i + nrho*(j+3 + neps*k));
                  
                  double eps_expol = 0.0;

                  if ( isnan(alltables[idx]) && !isnan(alltables[idxp1]) )
                  {
                     double rho  = logeps[j];
                     double rho0 = logeps[j+1];
                     double rho1 = logeps[j+2];
                     double rho2 = logeps[j+3];
                     double eps0 = alltables[idxp1];
                     double eps1 = alltables[idxp2];
                     double eps2 = alltables[idxp3];
                     eps_expol = eps0;
                     //eps_expol = (eps1-eps0)/(rho1-rho0)*(rho-rho0) + eps0;
                     //eps_expol   = eps0 * (rho - rho1)*(rho - rho2) / ( (rho0 - rho1)*(rho0 - rho2) )
                     //            + eps1 * (rho - rho0)*(rho - rho2) / ( (rho1 - rho0)*(rho1 - rho2) )
                     //            + eps2 * (rho - rho0)*(rho - rho1) / ( (rho2 - rho0)*(rho2 - rho1) );
                     eos_table[k][j][i][iv] = eps_expol;
                     alltables[idx] = eps_expol;           
   
                  }
                  else if ( isnan(alltables[idx]) && !isnan(alltables[idxm1]) )
                  {
                     double rho  = logeps[j];
                     double rho0 = logeps[j-1];
                     double rho1 = logeps[j-2];
                     double rho2 = logeps[j-3];
                     double eps0 = alltables[idxm1];
                     double eps1 = alltables[idxm2];
                     double eps2 = alltables[idxm3];
                     //eps_expol = eps0;
                     eps_expol = (eps1-eps0)/(rho1-rho0)*(rho-rho0) + eps0;
                     //eps_expol   = eps0 * (rho - rho1)*(rho - rho2) / ( (rho0 - rho1)*(rho0 - rho2) )
                     //            + eps1 * (rho - rho0)*(rho - rho2) / ( (rho1 - rho0)*(rho1 - rho2) )
                     //            + eps2 * (rho - rho0)*(rho - rho1) / ( (rho2 - rho0)*(rho2 - rho1) );
                     
                     eos_table[k][j][i][iv] = eps_expol;
                     alltables[idx] = eps_expol;
                  }
               }
            }
         }
      }
   }
   
   for (int N_exp=0; N_exp<N_exp_max; N_exp++) {
      for (int iv=0; iv<3; iv++) 
      {
         for (int k=0; k<nye_mode; k++)
         {
            for (int j=2; j<nmode-2; j++)
            {
               for (int i=0; i<nrho_mode; i++)
               {
                  int idx   = iv + 3*(i + nrho_mode*(j + nmode*k));
                  
                  int idxm3 = iv + 3*(i + nrho_mode*(j-3 + nmode*k));
                  int idxm2 = iv + 3*(i + nrho_mode*(j-2 + nmode*k));
                  int idxm1 = iv + 3*(i + nrho_mode*(j-1 + nmode*k));
                  int idxp1 = iv + 3*(i + nrho_mode*(j+1 + nmode*k));
                  int idxp2 = iv + 3*(i + nrho_mode*(j+2 + nmode*k));
                  int idxp3 = iv + 3*(i + nrho_mode*(j+3 + nmode*k));
                  
                  double eps_expol = 0.0;
   
                  if ( isnan(alltables_mode[idx]) && !isnan(alltables_mode[idxp1]) )
                  {
                     double rho  = logrho_mode[j];
                     double rho0 = logrho_mode[j+1];
                     double rho1 = logrho_mode[j+2];
                     double rho2 = logrho_mode[j+3];
                     double eps0 = alltables_mode[idxp1];
                     double eps1 = alltables_mode[idxp2];
                     double eps2 = alltables_mode[idxp3];
                     //eps_expol = eps0;
                     //eps_expol   = eps0 * (rho - rho1)*(rho - rho2) / ( (rho0 - rho1)*(rho0 - rho2) )
                     //            + eps1 * (rho - rho0)*(rho - rho2) / ( (rho1 - rho0)*(rho1 - rho2) )
                     //            + eps2 * (rho - rho0)*(rho - rho1) / ( (rho2 - rho0)*(rho2 - rho1) );
                     eps_expol = (eps1-eps0)/(rho1-rho0)*(rho-rho0) + eps0;
                     eos_table_mode[k][j][i][iv] = eps_expol;
                     alltables_mode[idx] = eps_expol;             
   
                  }
                  else if ( isnan(alltables_mode[idx]) && !isnan(alltables_mode[idxm1]) )
                  {
                     double rho  = logrho_mode[i];
                     double rho0 = logrho_mode[i-1];
                     double rho1 = logrho_mode[i-2];
                     double rho2 = logrho_mode[i-3];
                     double eps0 = alltables_mode[idxm1];
                     double eps1 = alltables_mode[idxm2];
                     double eps2 = alltables_mode[idxm3];
                     //eps_expol = eps0;
                     //eps_expol   = eps0 * (rho - rho1)*(rho - rho2) / ( (rho0 - rho1)*(rho0 - rho2) )
                     //             + eps1 * (rho - rho0)*(rho - rho2) / ( (rho1 - rho0)*(rho1 - rho2) )
                     //             + eps2 * (rho - rho0)*(rho - rho1) / ( (rho2 - rho0)*(rho2 - rho1) );
                     eps_expol = (eps1-eps0)/(rho1-rho0)*(rho-rho0) + eps0;
                     eos_table_mode[k][j][i][iv] = eps_expol;
                     alltables_mode[idx] = eps_expol;
                  }
               }
            }
         }
      }
   }
   

   // Create a new file using the default properties
  HDF5_ERROR(file = H5Fcreate(nuceos_eps_table_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  

// Use these defines to easily write a lot of variables in the same way
#define WRITE_EOSTABLE_HDF5(NAME, OFF)                                                                \
   do                                                                                                 \
   {                                                                                                  \
      for(int i = 0; i < nye; i++)                                                                    \
         for(int j = 0; j < neps; j++)                                                                \
            for(int k = 0; k < nrho; k++)                                                             \
               wdata[i*neps*nrho + j*nrho + k] = eos_table[i][j][k][OFF];                             \
               HDF5_ERROR(H5LTmake_dataset(file, NAME, 3, dims, H5T_NATIVE_DOUBLE, wdata));           \
   } while (0)

#define WRITE_EOSTABLE_MODE_HDF5(NAME, OFF)                                                           \
   do                                                                                                 \
   {                                                                                                  \
      for(int i = 0; i < nye_mode; i++)                                                               \
         for(int j = 0; j < nmode; j++)                                                               \
            for(int k = 0; k < nrho_mode; k++)                                                        \
               wdata_mode[i*nmode*nrho_mode + j*nrho_mode + k] = eos_table_mode[i][j][k][OFF];        \
               HDF5_ERROR(H5LTmake_dataset(file, NAME, 3, dims_mode, H5T_NATIVE_DOUBLE, wdata_mode)); \
   } while (0)

#define WRITE_EOS_HDF5(NAME, var)                                                                     \
   do                                                                                                 \
   {                                                                                                  \
      hsize_t dims[1] = {1};                                                                          \
      HDF5_ERROR(space = H5Screate_simple(1, dims, NULL));                                            \
      HDF5_ERROR(dataset = H5Dcreate(file, NAME, H5T_INTEL_I32, space, H5P_DEFAULT,                   \
                 H5P_DEFAULT, H5P_DEFAULT));                                                          \
      HDF5_ERROR(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, var));              \
      HDF5_ERROR(H5Dclose(dataset));                                                                  \
      HDF5_ERROR(H5Sclose(space));                                                                    \
   } while (0)

#define WRITE_EOS_INIT_HDF5(NAME, nvar, var)                                                          \
   do                                                                                                 \
   {                                                                                                  \
      hsize_t dims[1] = {nvar};                                                                       \
      HDF5_ERROR(space = H5Screate_simple(1, dims, NULL));                                            \
      HDF5_ERROR(dataset = H5Dcreate(file, NAME, H5T_IEEE_F64LE, space, H5P_DEFAULT,                  \
                 H5P_DEFAULT, H5P_DEFAULT));                                                          \
      HDF5_ERROR(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var));           \
      HDF5_ERROR(H5Dclose(dataset));                                                                  \
      HDF5_ERROR(H5Sclose(space));                                                                    \
   } while (0)

   
   // Write size of tables
   WRITE_EOS_HDF5("pointsrho",       &nrho);
   WRITE_EOS_HDF5("pointsenergy",    &neps);
   WRITE_EOS_HDF5("pointsye",        &nye);

   WRITE_EOS_HDF5("pointsrho_mode",  &nrho_mode);
   WRITE_EOS_HDF5("points_mode",     &nmode);
   WRITE_EOS_HDF5("pointsye_mode",   &nye_mode);
   
   // writes alltables
   WRITE_EOSTABLE_HDF5("logpress",  0);
   WRITE_EOSTABLE_HDF5("logtemp",   1);
   WRITE_EOSTABLE_HDF5("entropy",   2);
   WRITE_EOSTABLE_HDF5("munu",      3);
   WRITE_EOSTABLE_HDF5("cs2",       4);
   // chemical potentials
   WRITE_EOSTABLE_HDF5("muhat",     5);
   WRITE_EOSTABLE_HDF5("mu_e",      6);
   WRITE_EOSTABLE_HDF5("mu_p",      7);
   WRITE_EOSTABLE_HDF5("mu_n",      8);
   // compositions
   WRITE_EOSTABLE_HDF5("Xa",        9);
   WRITE_EOSTABLE_HDF5("Xh",       10);
   WRITE_EOSTABLE_HDF5("Xn",       11);
   WRITE_EOSTABLE_HDF5("Xp",       12);
   // average nucleus
   WRITE_EOSTABLE_HDF5("Abar",     13);
   WRITE_EOSTABLE_HDF5("Zbar",     14);
   // gamma
   WRITE_EOSTABLE_HDF5("gamma",    15);
    
   // energy for the temperature mode
   WRITE_EOSTABLE_MODE_HDF5("logenergy_temp", 0);

   // energy for the entropy mode
   WRITE_EOSTABLE_MODE_HDF5("logenergy_entr", 1);

   // energy for the pressure mode
   WRITE_EOSTABLE_MODE_HDF5("logenergy_prss", 2);


   // Write additional tables and variable
   WRITE_EOS_INIT_HDF5("logrho",        nrho,       logrho);
   WRITE_EOS_INIT_HDF5("logenergy",     neps,       logeps);
   WRITE_EOS_INIT_HDF5("ye",            nye,        yes);
   WRITE_EOS_INIT_HDF5("energy_shift" , 1,          &energy_shift);
   WRITE_EOS_INIT_HDF5("logrho_mode",   nrho_mode,  logrho_mode);
   WRITE_EOS_INIT_HDF5("logtemp_mode",  nmode,      logtemp_mode);
   WRITE_EOS_INIT_HDF5("entr_mode",     nmode,      entr_mode);
   WRITE_EOS_INIT_HDF5("logprss_mode",  nmode,      logprss_mode);
   WRITE_EOS_INIT_HDF5("ye_mode",       nye_mode,   yes_mode);

   // Close and release resources
   HDF5_ERROR(H5Fclose(file));


   // Free memory
   free(wdata);
   wdata = NULL;

   free(wdata_mode);
   wdata_mode = NULL;
   
   for   (int i=0; i<nye; i++)
   {
      for   (int j=0; j<neps; j++)
      {
         for   (int k=0; k<nrho; k++)
         {
            free(eos_table[i][j][k]);
            eos_table[i][j][k] = NULL;
         }
         free(eos_table[i][j]);
         eos_table[i][j] = NULL;
      }
      free(eos_table[i]);
      eos_table[i] = NULL;
   }
   free(eos_table);
   eos_table = NULL;

   for (int i=0; i<nye_mode; i++)
   {
      for (int j=0; j<nmode; j++)
      {
         for (int k=0; k<nrho_mode; k++)
         {
         free(eos_table_mode[i][j][k]);
            eos_table_mode[i][j][k] = NULL;
         }
         free(eos_table_mode[i][j]);
         eos_table_mode[i][j] = NULL;
      }
      free(eos_table_mode[i]);
      eos_table_mode[i] = NULL;
   }
   free(eos_table_mode);
   eos_table_mode = NULL;

   return;
}
