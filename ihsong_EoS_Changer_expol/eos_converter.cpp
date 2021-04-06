#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nuc_eos.h"

int nrho;
int neps;
int nye;

int nrho_mode;
int nmode;
int nye_mode;

double *alltables;
double *alltables_mode;
double *logrho;
double *logeps;
double *yes;
double energy_shift;
double drho, drhoi;
double deps, depsi;
double dye, dyei;

double *logrho_mode;
double *logtemp_mode;
double *entr_mode;
double *logprss_mode;
double *yes_mode;


// outvars
double xrho, xtemp, xye, xeps;
double xent, xprs, xcs2;
double xmunu, xmuhat, xmu_e, xmu_p, xmu_n;
double xXa, xXh, xXn, xXp;
double xAbar, xZbar;
double xgamma;

int ivs_short[16];

int main(void) {
	
   // Read the orginal table
   nuc_eos_C_ReadTable("../LS220_eps3.h5");
   printf("nrho:  %d\n", nrho);
   printf("neps:  %d\n", neps);
   printf("nye:   %d\n", nye);
   
   // Write the energy-based EOS with one point extrapolation
   write_eos_table("../LS220_eps3_expol3.h5");
 
	return 0;
}

