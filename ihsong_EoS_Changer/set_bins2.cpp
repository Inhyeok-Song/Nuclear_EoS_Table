#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nuc_eos.h"

void set_bins2(int rho_bs, int temp_bs, int entr_bs, int prss_bs, int ye_bs, 
               double eos_rhomin,  double eos_rhomax,
					double eos_tempmin, double eos_tempmax, 
               double eos_entrmin, double eos_entrmax,
               double eos_prssmin, double eos_prssmax,
               double eos_yemin,   double eos_yemax) {
  

  // set the table for temp, entr modes
//	if (!(logrho_mode = (double*)malloc(rho_bs * sizeof(double)))) {
//	  fprintf(stderr, "Cannot allocate memory for logrho_mode\n");
//		abort();
//	}
//  double dlrho = log10(eos_rhomax/eos_rhomin)/(rho_bs-1);
//	for	(int i = 0; i<rho_bs-1; i++) {
//		logrho_mode[i] = log10(eos_rhomin) + dlrho*i;
//	}
//	logrho_mode[rho_bs-1] = log10(eos_rhomax);
  
		 
  if (!(logtemp_mode = (double*)malloc(temp_bs * sizeof(double)))) {
	 fprintf(stderr, "Cannot allocate memory for logtemp_mode\n");
	 abort();
  }
  double dlt = log10(eos_tempmax/eos_tempmin)/(temp_bs-1);
  for	(int i = 0; i<temp_bs-1; i++) {
    logtemp_mode[i] = log10(eos_tempmin) + dlt*i;
  }
  logtemp_mode[temp_bs-1] = log10(eos_tempmax);


  if (!(entr_mode = (double*)malloc(entr_bs * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for logentr_mode\n");
	 abort();
  }
  double dentr = (eos_entrmax-eos_entrmin)/(entr_bs-1);
  for	(int i = 0; i<entr_bs-1; i++) {
  	 entr_mode[i] = eos_entrmin + dentr*i;
  }
  entr_mode[entr_bs-1] = eos_entrmax;

	
//	if (!(yes_mode = (double*)malloc(ye_bs * sizeof(double)))) {
//	  fprintf(stderr, "Cannot allocate memory for yes_mode\n");
//		abort();
//	}
//	double dye = (eos_yemax-eos_yemin)/(ye_bs-1);
//	for	(int i = 0; i<ye_bs-1; i++) {
//		yes_mode[i] = eos_yemin + dye*i;
//	}
//	yes_mode[ye_bs-1] = eos_yemax;

	
  if (!(logprss_mode = (double*)malloc(prss_bs * sizeof(double)))) {
	 fprintf(stderr, "Cannot allocate memory for logpress_mode\n");
    abort();
  }
  double dlprss = log10(eos_prssmax/eos_prssmin)/(prss_bs-1);
  for	(int i = 0; i<prss_bs-1; i++) {
    logprss_mode[i] = log10(eos_prssmin) + dlprss*i;
  }
  logprss_mode[prss_bs-1] = log10(eos_prssmax);
	
  
  return;

}
