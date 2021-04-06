#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "nuc_eos.h"

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                \
  do                                                       \
  {                                                        \
    int _error_code = fn_call;                             \
      if (_error_code < 0)                                 \
      {                                                    \
        fprintf(stderr,                                    \
                "HDF5 call '%s' returned error code %d",   \
                #fn_call, _error_code);                    \
        abort();                                           \
      }                                                    \
  } while (0)


void write_eos_table( char *nuceos_eps_table_name ) 
{

  fprintf(stdout, "*******************************\n");
  fprintf(stdout, "Writing nuc_eos table file:\n");
  fprintf(stdout, "%s\n", nuceos_eps_table_name);
  fprintf(stdout, "*******************************\n");


  hid_t    file;
  hid_t    dataset;
  hid_t    space;
  hsize_t  dims[3]      = {nye2, neps2, nrho2};
  hsize_t  dims_mode[3] = {nye_mode, nmode, nrho_mode};

  
  // Allocating eos table
  double ****eos_table = NULL;  // Write buffer
  if (!(eos_table = (double****)malloc(nye2 * sizeof(double***))))
  {
    fprintf(stderr, "Cannot allocate memory for EOS table\n");
    abort();
  }
  for (int i=0; i<nye2; i++)
  {
    if (!(eos_table[i] = (double***)malloc(neps2 * sizeof(double**))))
    {
      fprintf(stderr, "Cannot allocate memory for EOS table\n");
      abort();
    }
    for (int j=0; j<neps2; j++)
    {
      if (!(eos_table[i][j] = (double**)malloc(nrho2 * sizeof(double*))))
      {
        fprintf(stderr, "Cannot allocate memory for EOS table\n");
        abort();
      }
      for (int k=0; k<nrho2; k++)
      {
        if (!(eos_table[i][j][k] = (double*)malloc((NTABLES-3) * sizeof(double))))
        {
          fprintf(stderr, "Cannot allocate memory for EOS table\n");
          abort();
        }
      }
    }
  }

  // Allocating eos table for modes (temp, entr, prss)
  double ****eos_table_mode = NULL;  // Write buffer
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
  if (!(wdata = (double*)malloc(nye2 * neps2 * nrho2 * sizeof(double))))
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
  int keytemp = 0; // Using energy mode
  int keyerr;
  const double prec = 1.e-10; // tolerance

  xtemp  = pow(10.0, logtemp[5]);
  for (int i=0; i<nye2; i++)
  {
    for (int j=0; j<neps2; j++)
    {
      for (int k=0; k<nrho2; k++)
      {
        xrho  = pow(10.0, logrho2[k]);
        xeps  = pow(10.0, logeps2[j]) - energy_shift;
        xye   = yes2[i];
        nuc_eos_C_short( xrho, &xtemp, xye, &xeps, &xprs, &xent, &xcs2, &xmunu,
                         &xmuhat, &xmu_e, &xmu_p, &xmu_n, &xXa, &xXh, &xXn, &xXp,
                         &xAbar, &xZbar, &xgamma, &dpde, keytemp, &keyerr, prec );
        if (keyerr == 0)
        {
          eos_table[i][j][k][0]  = log10(xprs);
          eos_table[i][j][k][1]  = log10(xtemp);
          eos_table[i][j][k][2]  = xent;
          eos_table[i][j][k][3]  = xmunu;
          eos_table[i][j][k][4]  = xcs2;
          eos_table[i][j][k][5]  = xmuhat;
          eos_table[i][j][k][6]  = xmu_e;
          eos_table[i][j][k][7]  = xmu_p;
          eos_table[i][j][k][8]  = xmu_n;
          eos_table[i][j][k][9]  = xXa;
          eos_table[i][j][k][10] = xXh;
          eos_table[i][j][k][11] = xXn;
          eos_table[i][j][k][12] = xXp;
          eos_table[i][j][k][13] = xAbar;
          eos_table[i][j][k][14] = xZbar;
          eos_table[i][j][k][15] = xgamma;
        }  
        else //if (ii == ntemp-2 && keyerr != 0) 
        {
          //  wrinting NAN for invailid regions
          eos_table[i][j][k][0]  = NAN;  // 0.0;
          eos_table[i][j][k][1]  = NAN;  // 0.0;
          eos_table[i][j][k][2]  = NAN;  // 0.0;
          eos_table[i][j][k][3]  = NAN;  // 0.0;
          eos_table[i][j][k][4]  = NAN;  // 0.0;
          eos_table[i][j][k][5]  = NAN;  // 0.0;
          eos_table[i][j][k][6]  = NAN;  // 0.0;
          eos_table[i][j][k][7]  = NAN;  // 0.0;
          eos_table[i][j][k][8]  = NAN;  // 0.0;
          eos_table[i][j][k][9]  = NAN;  // 0.0;
          eos_table[i][j][k][10] = NAN;  // 0.0;
          eos_table[i][j][k][11] = NAN;  // 0.0;
          eos_table[i][j][k][12] = NAN;  // 0.0;
          eos_table[i][j][k][13] = NAN;  // 0.0;
          eos_table[i][j][k][14] = NAN;  // 0.0;
          eos_table[i][j][k][15] = NAN;  // 0.0;
          xtemp  = pow(10.0, logtemp[5]);
        }
      }
    }
  }
  fprintf(stdout, "energy mode is done\n");

  //  energy array for temperature mode 
  keytemp = 1;
  xtemp  = pow(10.0, logtemp[5]);
  for (int i=0; i<nye_mode; i++)
  {
    for (int j=0; j<nmode; j++)
    {
      for (int k=0; k<nrho_mode; k++)
      {
        xrho  = pow(10.0, logrho_mode[k]);
        xtemp = pow(10.0, logtemp_mode[j]);
        xye   = yes_mode[i];
        nuc_eos_C_short( xrho, &xtemp, xye, &xeps, &xprs, &xent, &xcs2, &xmunu,
                         &xmuhat, &xmu_e, &xmu_p, &xmu_n, &xXa, &xXh, &xXn, &xXp,
                         &xAbar, &xZbar, &xgamma, &dpde, keytemp, &keyerr, prec );
        if (keyerr == 0)
        {
          eos_table_mode[i][j][k][0] = log10(MAX(1.0, xeps + energy_shift));
        }
        else 
        {
          eos_table_mode[i][j][k][0] = NAN;
          xtemp  = pow(10.0, logtemp[5]);
        }
      }
    }
  }
  fprintf(stdout, "temperature mode is done\n");

  //  energy array for entropy mode
  keytemp = 2;
  xtemp  = pow(10.0, logtemp[5]);
  for (int i=0; i<nye_mode; i++)
  {
    for (int j=0; j<nmode; j++)
    {
      for (int k=0; k<nrho_mode; k++)
      {
        xrho = pow( 10.0, logrho_mode[k] );
        xent = entr_mode[j];
        xye  = yes_mode[i];
        nuc_eos_C_short( xrho, &xtemp, xye, &xeps, &xprs, &xent, &xcs2, &xmunu,
                    &xmuhat, &xmu_e, &xmu_p, &xmu_n, &xXa, &xXh, &xXn, &xXp,
                    &xAbar, &xZbar, &xgamma, &dpde, keytemp, &keyerr, prec );
        if (keyerr == 0)
        {
          eos_table_mode[i][j][k][1] = log10(xeps + energy_shift);
        }
        else //if (ii == ntemp-2 && keyerr !=0)
        {
          eos_table_mode[i][j][k][1] = NAN;
          xtemp  = pow(10.0, logtemp[5]);
        }
      }
    }
  }
  fprintf(stdout, "entropy mode is done\n");

  // energy array for pressure mode
  keytemp = 3;
  xtemp = pow(10.0, logtemp[5]);
  for (int i=0; i<nye_mode; i++)
  {
    for (int j=0; j<nmode; j++)
    {
      for (int k=0; k<nrho_mode; k++)
      {
        xrho  = pow( 10.0, logrho_mode[k] );
        xprs  = pow( 10.0, logprss_mode[j] );
        xye   = yes_mode[i];
        nuc_eos_C_short( xrho, &xtemp, xye, &xeps, &xprs, &xent, &xcs2, &xmunu,
                         &xmuhat, &xmu_e, &xmu_p, &xmu_n, &xXa, &xXh, &xXn, &xXp,
                         &xAbar, &xZbar, &xgamma, &dpde, keytemp, &keyerr, prec );
        if (keyerr == 0)
        {
          eos_table_mode[i][j][k][2] = log10(xeps + energy_shift);
        }
        else
        {
          eos_table_mode[i][j][k][2] = NAN;
          xtemp = pow(10.0, logtemp[1]);
        }
      }
    }
  }
  fprintf(stdout, "pressure mode is done\n");


  // extend pressure mode boundaries
  //int count_it = 0;
  int iii = 0;
  while (iii<3) {
    iii++;
    xtemp = pow(10.0, logtemp[5]);
    for (int i=0; i<nye_mode; i++)
    {
      for (int j=3; j<nmode-4; j++)
      {
        for (int k=0; k<nrho_mode; k++)
        {
          if ( isnan(eos_table_mode[i][j][k][2]) && !isnan(eos_table_mode[i][j+1][k][2]) )
          {
            keyerr = 0;
            xrho  = pow( 10.0, logrho_mode[k] );
            xprs  = pow( 10.0, logprss_mode[j+iii] );
            xtemp = pow( 10.0, logtemp[5] );
            xye   = yes_mode[i];
            nuc_eos_C_short( xrho, &xtemp, xye, &xeps, &xprs, &xent, &xcs2, &xmunu,
                             &xmuhat, &xmu_e, &xmu_p, &xmu_n, &xXa, &xXh, &xXn, &xXp,
                             &xAbar, &xZbar, &xgamma, &dpde, 3, &keyerr, prec );
            if (keyerr == 0)
            {
              double dp = pow(10.0, logprss_mode[j+iii]) - pow(10.0, logprss_mode[j]);
              eos_table_mode[i][j][k][2] = log10(pow(10.0, eos_table_mode[i][j+iii][k][2]) - dp/dpde);
            }
            else
            {
              eos_table_mode[i][j][k][2] = NAN;
            }

          }
          else if ( isnan(eos_table_mode[i][j][k][2]) && !isnan(eos_table_mode[i][j-1][k][2]) )
          {
            keyerr = 0;
            xrho  = pow( 10.0, logrho_mode[k] );
            xprs  = pow( 10.0, logprss_mode[j-iii] );
            xtemp = pow( 10.0, logtemp[ntemp-5] );
            xye   = yes_mode[i];
            nuc_eos_C_short( xrho, &xtemp, xye, &xeps, &xprs, &xent, &xcs2, &xmunu,
                             &xmuhat, &xmu_e, &xmu_p, &xmu_n, &xXa, &xXh, &xXn, &xXp,
                             &xAbar, &xZbar, &xgamma, &dpde, 3, &keyerr, prec );
            if (keyerr == 0) 
            {
              double dp = pow(10.0, logprss_mode[j]) - pow(10.0, logprss_mode[j-iii]);
              eos_table_mode[i][j][k][2] = log10(pow(10.0, eos_table_mode[i][j-iii][k][2]) + dp/dpde);
            }
            else
            {
              eos_table_mode[i][j][k][2] = NAN;
            }
          }
        }
      }
    }
  }

  // Create a new file using the default properties
  HDF5_ERROR(file = H5Fcreate(nuceos_eps_table_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  

// Use these defines to easily write a lot of variables in the same way
#define WRITE_EOSTABLE_HDF5(NAME, OFF)                                                           \
  do                                                                                             \
  {                                                                                              \
    for(int i = 0; i < nye2; i++)                                                                \
      for(int j = 0; j < neps2; j++)                                                             \
        for(int k = 0; k < nrho2; k++)                                                           \
          wdata[i*neps2*nrho2 + j*nrho2 + k] = eos_table[i][j][k][OFF];                          \
          HDF5_ERROR(H5LTmake_dataset(file, NAME, 3, dims, H5T_NATIVE_DOUBLE, wdata));           \
  } while (0)

#define WRITE_EOSTABLE_MODE_HDF5(NAME, OFF)                                                      \
  do                                                                                             \
  {                                                                                              \
    for(int i = 0; i < nye_mode; i++)                                                            \
      for(int j = 0; j < nmode; j++)                                                             \
        for(int k = 0; k < nrho_mode; k++)                                                       \
          wdata_mode[i*nmode*nrho_mode + j*nrho_mode + k] = eos_table_mode[i][j][k][OFF];        \
          HDF5_ERROR(H5LTmake_dataset(file, NAME, 3, dims_mode, H5T_NATIVE_DOUBLE, wdata_mode)); \
  } while (0)

#define WRITE_EOS_HDF5(NAME, var)                                                                \
  do                                                                                             \
  {                                                                                              \
    hsize_t dims[1] = {1};                                                                       \
    HDF5_ERROR(space = H5Screate_simple(1, dims, NULL));                                         \
    HDF5_ERROR(dataset = H5Dcreate(file, NAME, H5T_INTEL_I32, space, H5P_DEFAULT,                \
               H5P_DEFAULT, H5P_DEFAULT));                                                       \
    HDF5_ERROR(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, var));           \
    HDF5_ERROR(H5Dclose(dataset));                                                               \
    HDF5_ERROR(H5Sclose(space));                                                                 \
  } while (0)

#define WRITE_EOS_INIT_HDF5(NAME, nvar, var)                                                     \
  do                                                                                             \
  {                                                                                              \
    hsize_t dims[1] = {nvar};                                                                    \
    HDF5_ERROR(space = H5Screate_simple(1, dims, NULL));                                         \
    HDF5_ERROR(dataset = H5Dcreate(file, NAME, H5T_IEEE_F64LE, space, H5P_DEFAULT,               \
               H5P_DEFAULT, H5P_DEFAULT));                                                       \
    HDF5_ERROR(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var));        \
    HDF5_ERROR(H5Dclose(dataset));                                                               \
    HDF5_ERROR(H5Sclose(space));                                                                 \
  } while (0)

   
  // Write size of tables
  WRITE_EOS_HDF5("pointsrho",      &nrho2);
  WRITE_EOS_HDF5("pointsenergy",   &neps2);
  WRITE_EOS_HDF5("pointsye",       &nye2);
  WRITE_EOS_HDF5("pointsrho_mode", &nrho_mode);
  WRITE_EOS_HDF5("points_mode",    &nmode);
  WRITE_EOS_HDF5("pointsye_mode",  &nye_mode);

  // writes alltables
  WRITE_EOSTABLE_HDF5("logpress", 0);
  WRITE_EOSTABLE_HDF5("logtemp",  1);
  WRITE_EOSTABLE_HDF5("entropy",  2);
  WRITE_EOSTABLE_HDF5("munu",     3);
  WRITE_EOSTABLE_HDF5("cs2",      4);
  // chemical potentials
  WRITE_EOSTABLE_HDF5("muhat",    5);
  WRITE_EOSTABLE_HDF5("mu_e",     6);
  WRITE_EOSTABLE_HDF5("mu_p",     7);
  WRITE_EOSTABLE_HDF5("mu_n",     8);
  // compositions
  WRITE_EOSTABLE_HDF5("Xa",       9);
  WRITE_EOSTABLE_HDF5("Xh",      10);
  WRITE_EOSTABLE_HDF5("Xn",      11);
  WRITE_EOSTABLE_HDF5("Xp",      12);
  // average nucleus
  WRITE_EOSTABLE_HDF5("Abar",    13);
  WRITE_EOSTABLE_HDF5("Zbar",    14);
  // gamma
  WRITE_EOSTABLE_HDF5("gamma",   15);
    
  // energy for the temperature mode
  WRITE_EOSTABLE_MODE_HDF5("logenergy_temp", 0);
  
  // energy for the entropy mode
  WRITE_EOSTABLE_MODE_HDF5("logenergy_entr", 1);
  
  // energy for the pressure mode
  WRITE_EOSTABLE_MODE_HDF5("logenergy_prss", 2);


  // Write additional tables and variable
  WRITE_EOS_INIT_HDF5("logrho",        nrho2,      logrho2);
  WRITE_EOS_INIT_HDF5("logenergy",     neps2,      logeps2);
  WRITE_EOS_INIT_HDF5("ye",            nye2,       yes2);
  WRITE_EOS_INIT_HDF5("energy_shift" , 1,          &energy_shift2);
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
   
  for   (int i=0; i<nye2; i++)
  {
    for   (int j=0; j<neps2; j++)
    {
      for   (int k=0; k<nrho2; k++)
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
