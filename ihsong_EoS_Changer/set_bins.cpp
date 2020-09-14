#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nuc_eos.h"

void set_bins(int rho_bs, int eps_bs, int ye_bs, double eos_rhomin, double eos_rhomax,
					double eos_epsmin, double eos_epsmax, double eos_yemin, double eos_yemax) {
  
  
	// set eps
	
	// Allocate memory for logrho, logeps and ye
	if (!(logrho2 = (double*)malloc(rho_bs * sizeof(double)))) {
	  fprintf(stderr, "Cannot allocate memery for logeps\n");
		abort();
	}
	printf("*******************************\n");
	printf("Setting rho\n");
	printf("rho_binsize: %d\n", rho_bs);
  printf("rho_min: %9.3E\n", eos_rhomin); 
	printf("rho_max: %9.3E\n", eos_rhomax);
  double dlrho = log10(eos_rhomax/eos_rhomin)/(rho_bs-1);
	for	(int i = 0; i<rho_bs-1; i++) {
		logrho2[i] = log10(eos_rhomin) + dlrho*i;
	}
	logrho2[rho_bs-1] = log10(eos_rhomax);
  
	//for (int i = 0; i<rho_bs; i++) {
	//  printf("logrho[%d]: %f\n", i, logrho2[i]);
	//}

		 
	if (!(logeps2 = (double*)malloc(eps_bs * sizeof(double)))) {
	  fprintf(stderr, "Cannot allocate memery for logeps\n");
		abort();
	}
	printf("*******************************\n");
	printf("Setting eps\n");
	printf("eps_binsize: %d\n", eps_bs);
  printf("eps_min: %9.3E\n", eos_epsmin); 
	printf("eps_max: %9.3E\n", eos_epsmax);
  double dleps = log10(eos_epsmax/eos_epsmin)/(eps_bs-1);
	for	(int i = 0; i<eps_bs-1; i++) {
		logeps2[i] = log10(eos_epsmin) + dleps*i;
	}
	logeps2[eps_bs-1] = log10(eos_epsmax);

	//for (int i = 0; i<eps_bs; i++) {
	//  printf("logeps[%d]: %f\n", i, logeps2[i]);
	//}
	
	
	if (!(yes2 = (double*)malloc(ye_bs * sizeof(double)))) {
	  fprintf(stderr, "Cannot allocate memery for logeps\n");
		abort();
	}
	printf("*******************************\n");
	printf("Setting ye\n");
	printf("ye_binsize: %d\n", ye_bs);
  printf("ye_min: %9.3E\n", eos_yemin); 
	printf("ye_max: %9.3E\n", eos_yemax);
	double dye = (eos_yemax-eos_yemin)/(ye_bs-1);
	for	(int i = 0; i<ye_bs-1; i++) {
		yes2[i] = eos_yemin + dye*i;
	}
	yes2[ye_bs-1] = eos_yemax;
  
	//for (int i = 0; i<ye_bs; i++) {
	//  printf("ye[%d]: %f\n", i, yes2[i]);
	//}

	
  return;
}

