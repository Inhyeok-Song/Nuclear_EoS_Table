#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include "nuc_eos.h"

double linterp2D(double* xs, double* ys, double* fs, double x, double y);

void bisection(double lr, double lt0, double ye, double leps0, double* ltout,
	       int iv, double prec, int* keyerrt);


void nuc_eos_C_findtemp_pressure( double lr, double lt0, double ye, double lprsin, double *ltout, 
                                  double prec, int *keyerrt ) {
  
   // local vars
  int itmax = 20; // use at most 20 iterations, then go to bisection
  double dlprsdlt; // derivative dlogprs/dlogT
  double ldt;
  double lprs,lprs0,lprs1; // temp vars for prs
  double ltn, lt, lt1; // temp vars for temperature
  double ltmax = logtemp[ntemp-1]; // max temp
  double ltmin = logtemp[0]; // min temp

  // setting up some vars
  *keyerrt = 0;
  lprs0 = lprsin;
  lprs1 = lprs0;
  lt = lt0;
  lt1 = lt;

  // step 1: do we already have the right temperature
  nuc_eos_C_linterp_for_prss(lr,lt,ye,&lprs,alltables,nrho,ntemp,nye,
			     logrho,logtemp,yes,&dlprsdlt);

  if(fabs(lprs-lprs0) < prec*fabs(lprs0)) {
    *ltout = lt0;
    return;
  }
  lt1 = lt;
  lprs1 = lprs;

  int it = 0;
  while(it < itmax) {

    // step 2: check if the two bounding values of the temperature
    //         give prs values that enclose the new prs.
    int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-10) * dtempi),1),ntemp-1);
    int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-10) * drhoi),1),nrho-1);
    int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-10) * dyei),1),nye-1);

    double prst1, prst2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 0 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1
      ifs = 0 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2
      ifs = 0 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 0 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
      fs[3] = alltables[ifs];

      prst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 0 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1
      ifs = 0 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2
      ifs = 0 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 0 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
      fs[3] = alltables[ifs];

      prst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }

    // Check if we are already bracketing the input internal
    // pressure. If so, interpolate for new T.
    if(lprs0 >= prst1 && lprs0 <= prst2) {

      *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (prst2 - prst1) *
	(lprs0 - prst1) + logtemp[itemp-1];

#if DEBUG
      fprintf(stderr,"it: %d, bracketed solution\n", it);
#endif
      //return;

    }

    // well, then do a Newton-Raphson step
    ldt = -(lprs - lprs0) / dlprsdlt;
    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt1 = lt;
    lt = ltn;
    lprs1 = lprs;

    nuc_eos_C_linterp_for_prss(lr,lt,ye,&lprs,alltables,nrho,ntemp,nye,
			       logrho,logtemp,yes,&dlprsdlt);
#if DEBUG
    fprintf(stderr,"findtemp it: %d, err: %15.6E \n", it, fabs((lprs-lprs0) / lprs0));
#endif

    if(fabs(lprs-lprs0) < prec*fabs(lprs0)) {
      *ltout = lt;
      return;
    }

    // if we are closer than 10^-3  to the
    // root (prs-prs0)=0, we are switching to
    // the secant method, since the table is rather coarse and the
    // derivatives may be garbage.
    if(fabs(lprs-lprs0) < 1.0e-3*fabs(lprs0)) {
      dlprsdlt = (lprs-lprs1)/(lt-lt1);
    }

    it++;
  }

  if(it >= itmax-1) {
    // try bisection
#if DEBUG
    fprintf(stderr,"trying bisection\n");
#endif
    bisection(lr, lt0, ye, lprs0, ltout, 0, prec, keyerrt);
#if DEBUG
    fprintf(stderr,"bisection keyerrt: %d\n",*keyerrt);
#endif
    return;
  }

  fprintf(stderr,"We should never reach this point! Aborting!\n");
  abort();

  return;
}


void nuc_eos_C_findtemp_entropy(double lr, double lt0, double ye, 
			double xs, double *ltout, double prec, int *keyerrt) {

    // local vars
    double dlepsdlt; // derivative dlogeps/dlogT
    int itmax = 20;
    double s, lt, ldt, tol, d1,d2,d3,s0, s1, lt1;
    double ltn, tinput;
    double ltmax = logtemp[ntemp-1]; // max temp
    double ltmin = logtemp[0]; // min temp

    *keyerrt = 0;
    tol = prec;
    lt = lt0;
    lt1 = lt;
    s0 = xs;
    s1 = s0;

    // step 1: do we already have the right temperature
    nuc_eos_C_linterp_for_entr(lr,lt,ye,&s,alltables,nrho,ntemp,nye,
			     logrho,logtemp,yes,&dlepsdlt);

    if (fabs(s-s0) < tol*fabs(s0)) 
    {
        return;
    }

    lt1 = lt;
    s1 = s;
    
    int it = 0;
    while(it < itmax) {

    // step 2: check if the two bounding values of the entropy
    //         give eps values that enclose the new entr.
    int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-10) * dtempi),1),ntemp-1); 
    int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-10) * drhoi),1),nrho-1); 
    int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-10) * dyei),1),nye-1); 

    double ss1, ss2;
    // lower entropy
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      ss1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    // upper entripy
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      ss2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    
    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if(s0 >= ss1 && s0 <= ss2) {
      
      *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (ss2 - ss1) * 
	           (s0 - ss1) + logtemp[itemp-1];
#if DEBUG
      fprintf(stderr,"it: %d, bracketed solution\n", it);
#endif
      return;

    }

    // well, then do a Newton-Raphson step
    ldt = -(s - s0) / dlepsdlt;
    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt1 = lt;
    lt = ltn;
    s1 = s;

    nuc_eos_C_linterp_for_entr(lr,lt,ye,&s,alltables,nrho,ntemp,nye,
			       logrho,logtemp,yes,&dlepsdlt);
#if DEBUG
    fprintf(stderr,"findtemp it: %d, err: %15.6E \n", it, fabs((s-s0) / s0));
#endif

    if(fabs(s-s0) < prec*fabs(s0)) {
      *ltout = lt;
      return;
    }
    
    // if we are closer than 10^-3  to the 
    // root (eps-eps0)=0, we are switching to 
    // the secant method, since the table is rather coarse and the
    // derivatives may be garbage.
    if(fabs(s-s0) < 1.0e-3*fabs(s0)) {
      dlepsdlt = (s-s1)/(lt-lt1);
    }

    it++;
  }

  if(it >= itmax - 1) {
    // try bisection
#if DEBUG
    fprintf(stderr,"trying bisection\n");
#endif
    bisection(lr, lt0, ye, s0, ltout, 2, prec, keyerrt);
#if DEBUG
    fprintf(stderr,"bisection keyerrt: %d\n",*keyerrt);
#endif
    return;
  }

  fprintf(stderr,"We should never reach this point! Aborting!\n");
  abort();
}


void nuc_eos_C_findtemp(double lr, double lt0, double ye, 
			double lepsin, double prec, double *ltout,
			int *keyerrt) {


//while (prec < 2.0E-2) {
 
  // local vars
  int itmax = 20; // use at most 10 iterations, then go to bisection
  double dlepsdlt; // derivative dlogeps/dlogT
  double ldt;
  double leps,leps0,leps1; // temp vars for eps
  double ltn, lt, lt1; // temp vars for temperature
  double ltmax = logtemp[ntemp-1]; // max temp
  double ltmin = logtemp[0]; // min temp

  // setting up some vars
  *keyerrt = 0;
  leps0 = lepsin;
  leps1 = leps0;
  lt = lt0;
  lt1 = lt;

  // step 1: do we already have the right temperature
  nuc_eos_C_linterp_for_temp(lr,lt,ye,&leps,alltables,nrho,ntemp,nye,
			     logrho,logtemp,yes,&dlepsdlt);

  if(fabs(leps-leps0) < prec*fabs(leps0)) {
    *ltout = lt0;
    return;
  }
  lt1 = lt;
  leps1 = leps;

  int it = 0;
  while(it < itmax) {

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-10) * dtempi),1),ntemp-1); 
    int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-10) * drhoi),1),nrho-1); 
    int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-10) * dyei),1),nye-1); 

    double epst1, epst2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      epst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      epst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    
    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if(leps0 >= epst1 && leps0 <= epst2) {
      
      *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (epst2 - epst1) * 
	(leps0 - epst1) + logtemp[itemp-1];
#if DEBUG
      fprintf(stderr,"it: %d, bracketed solution\n", it);
#endif
      return;

    }

    // well, then do a Newton-Raphson step
    ldt = -(leps - leps0) / dlepsdlt;
    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt1 = lt;
    lt = ltn;
    leps1 = leps;

    nuc_eos_C_linterp_for_temp(lr,lt,ye,&leps,alltables,nrho,ntemp,nye,
			       logrho,logtemp,yes,&dlepsdlt);
#if DEBUG
    fprintf(stderr,"findtemp it: %d, err: %15.6E \n", it, fabs((leps-leps0) / leps0));
#endif

    if(fabs(leps-leps0) < prec*fabs(leps0)) {
      *ltout = lt;
      return;
    }
    
    // if we are closer than 10^-3  to the 
    // root (eps-eps0)=0, we are switching to 
    // the secant method, since the table is rather coarse and the
    // derivatives may be garbage.
    if(fabs(leps-leps0) < 1.0e-3*fabs(leps0)) {
      dlepsdlt = (leps-leps1)/(lt-lt1);
    }

    it++;
  }

  if(it >= itmax - 1) {
    // try bisection
#if DEBUG
    fprintf(stderr,"trying bisection\n");
#endif
    bisection(lr, lt0, ye, leps0, ltout, 1, prec, keyerrt);
#if DEBUG
    fprintf(stderr,"bisection keyerrt: %d\n",*keyerrt);
#endif
    return;
  }
//  printf("%5.3E\n", prec);
//  prec = prec*1.2;
//}

  //fprintf(stderr,"We should never reach this point! Aborting!\n");
  //abort();

}


double linterp2D(double* xs, double* ys, double* fs, double x, double y)
{

  //  2     3 
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid
  
  double t1 = (fs[1]-fs[0])/(xs[1]-xs[0]) * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])/(xs[1]-xs[0]) * (x - xs[0]) + fs[2];

  return (t2 - t1)/(ys[1]-ys[0]) * (y-ys[0]) + t1;
}

void bisection(double lr, double lt0, double ye, double leps0, double* ltout,
	       int iv, double prec, int* keyerrt) {

  // iv is the index of the table variable we do the bisection on

  int bcount = 0; 
  int maxbcount = 30;
  int itmax = 50;
  
  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = logtemp[0];
  double ltmax = logtemp[ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double dlepsdlt;
  double f1a[3];
  double f2a[3];

  // prepare
  lt  = lt0;
  lt1 = log10( MIN(pow(10.0,ltmax),1.2*(pow(10.0,lt0))) );
  lt2 = log10( MAX(pow(10.0,ltmin),0.8*(pow(10.0,lt0))) );

  int nvars = 3;
  nuc_eos_C_linterp_some(lr, lt1, ye, f1a, alltables, 
  			 ivs_short, nrho, ntemp, nye, nvars,
  			 logrho, logtemp, yes);

  nuc_eos_C_linterp_some(lr, lt2, ye, f2a, alltables, 
			 ivs_short, nrho, ntemp, nye, nvars,
  			 logrho, logtemp, yes);


  f1=f1a[iv]-leps0;
  f2=f2a[iv]-leps0;
  
  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
#if 0
  int ifixdeg = 0;
  int ifixdeg_max = 20;
#endif
  while(f1*f2 >= 0.0) {
    
    lt1 = log10( MIN(pow(10.0,ltmax),1.2*(pow(10.0,lt1))) );
    lt2 = log10( MAX(pow(10.0,ltmin),0.8*(pow(10.0,lt2))) );
    nuc_eos_C_linterp_some(lr, lt1, ye, f1a, alltables, 
			   ivs_short, nrho, ntemp, nye, nvars,
			   logrho, logtemp, yes);
    nuc_eos_C_linterp_some(lr, lt2, ye, f2a, alltables, 
			   ivs_short, nrho, ntemp, nye, nvars,
			   logrho, logtemp, yes);

#if 0
    // special enforcement of eps(lt1)>eps(lt2)
    while(f1a < f2a && ifixdeg < ifixdeg_max) {
      lt1 = log10( MIN(pow(10.0,ltmax),1.2*(pow(10.0,lt1))) );
      nuc_eos_C_linterp_some(lr, lt1, ye, &f1a, alltables, 
			     ivs, nrho, ntemp, nye, nvars,
			     logrho, logtemp, yes);
      ifixdeg++;

    }
#endif


    f1=f1a[iv]-leps0;
    f2=f2a[iv]-leps0;



#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
	    bcount,f1,f2,lt1,lt2,f1a[iv],f2a[iv],leps0);
#endif

    bcount++;
    if(bcount >= maxbcount) {
      *keyerrt = 667;
      return;
    }

  }

  if(f1 < 0.0) {
    lt = lt1;
    dlt = log10( pow(10.0,lt2) - pow(10.0,lt1) );
  } else {
    lt = lt2;
    dlt = log10( pow(10.0,lt1) - pow(10.0,lt2) );
  }

  //  fprintf(stderr,"%15.6E %15.6E %15.6E %15.6E\n",lt,dlt,lt1,lt2);
  //fprintf(stderr,"%15.6E %15.6E \n",pow(10.0,lt2) - pow(10.0,lt1),
  //	  pow(10.0,lt1) - pow(10.0,lt2));

  int it;
  for(it=0;it<itmax;it++) {
    dlt = log10( pow(10.0,dlt) * 0.5 );
    ltmid = log10( pow(10.0,lt) + pow(10.0,dlt) );
    nuc_eos_C_linterp_some(lr, ltmid, ye, f2a, alltables, 
			   ivs_short, nrho, ntemp, nye, nvars,
			   logrho, logtemp, yes);

    fmid=f2a[iv]-leps0;
    if(fmid <= 0.0) lt=ltmid;

#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	          it,fmid,dlt,ltmid);
#endif

    if(fabs(1.0-f2a[iv]/leps0) <= prec) {
      *ltout = ltmid;
      return;
    }

  }

  if(it >= itmax-1) {
    *keyerrt = 667;
    return;
  }

}
