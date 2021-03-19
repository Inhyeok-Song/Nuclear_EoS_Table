#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nuc_eos.h"

int nrho;
int nrho2;
int ntemp;
int neps2;
int nye;
int nye2;

int nrho_mode;
int ntemp_mode;
int nentr_mode;
int nprss_mode;
int nye_mode;
int nmode;

double *alltables;
double *logrho;
double *logrho2;
double *logtemp;
double *logeps2;
double *yes;
double *yes2;
double energy_shift;
double energy_shift2;
double drho, drhoi;
double dtemp, dtempi;
double deps, depsi;
double dye, dyei;

double *logrho_mode;
double *logtemp_mode;
double *entr_mode;
double *logprss_mode;
double *yes_mode;

// min and max values
double eos_rhomax, eos_rhomin;
double eos_tempmin, eos_tempmax;
double eos_epsmin, eos_epsmax;
double eos_yemin, eos_yemax;

// outvars
double xrho, xtemp, xye, xeps;
double xent, xprs, xcs2;
double xmunu, xmuhat, xmu_e, xmu_p, xmu_n;
double xXa, xXh, xXn, xXp;
double xAbar, xZbar;
double xgamma;
double dpde;

int ivs_short[19];

int main(void) {

  // Read the orginal table
  nuc_eos_C_ReadTable("../LS220.h5");
  printf("nrho: %d\n", nrho);
  printf("ntemp: %d\n", ntemp);
  printf("nye: %d\n", nye);

  // set the table size
  nrho2  = nrho;
  neps2  = 2*ntemp;
  nye2   = nye;

  // double rho_min = 1.E3;
  // double rho_max = 1.E16;
  // double eps_min = 1.E17; // pow(10.0, 15.889969377212918);
  // double eps_max = 1.E23; // pow(10.0, 33.17670042883674);
  // double ye_min  = 0.035;
  // double ye_max  = 0.55;

  double rho_min = pow( 10.0, logrho[0] );
  double rho_max = pow( 10.0, logrho[nrho-1] );
  double eps_min = 1.E17; //pow(10.0, 15.889969377212918); //
  double eps_max = 1.E21; //pow(10.0, 33.17670042883674); //
  double ye_min  = yes[0];
  double ye_max  = yes[nye-1];

  // set density, energy, Ye arrays for
  // the new table  
	set_bins( nrho2, neps2, nye2, rho_min, rho_max, 
					  eps_min, eps_max, ye_min, ye_max );

  // set tables for temp, entr modes
  double rho_mode_min  = pow( 10.0, logrho[0] );
  double rho_mode_max  = pow( 10.0, logrho[nrho-1] );
  double temp_mode_min = pow( 10.0, logtemp[0] );
  double temp_mode_max = pow( 10.0, logtemp[ntemp-1] );
  double entr_mode_min = 1.0E-3;
  double entr_mode_max = 3.0E1;
  double prss_mode_min = pow( 10.0, 20.0 );
  double prss_mode_max = pow( 10.0, 35.5 );
  double ye_mode_min   = yes[0];
  double ye_mode_max   = yes[nye-1];

  nmode        = 2*ntemp;
  nrho_mode    = nrho;
  ntemp_mode   = nmode;
  nentr_mode   = nmode;
  nprss_mode   = nmode;
  nye_mode     = nye;
  //logrho_mode  = logrho;
  //yes_mode     = yes;
  
  // set density, temperature, entropy, pressure and Ye array for
  // temperature, entropy and pressure mode
  set_bins2( nrho_mode, ntemp_mode, nentr_mode, nprss_mode, nye_mode,
             rho_mode_min, rho_mode_max, temp_mode_min, temp_mode_max,
             entr_mode_min, entr_mode_max, prss_mode_min, prss_mode_max,
             ye_mode_min, ye_mode_max );


  // Set new energy_shift
  energy_shift2 = energy_shift;
  
  
  // Write the energy-based EOS
  write_eos_table("../LS220_eps2.h5");


  return 0;
}
