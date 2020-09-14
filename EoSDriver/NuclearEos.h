//#ifndef __NUCLEAREOS_H__
//#define __NUCLEAREOS_H__

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 16
//#define DEBUG 1

// variables for nuclear Eos

extern int nrho;
extern int neps;
extern int nye;
extern int nmode;

extern double *alltables;
extern double *alltables_mode;
extern double *logrho; 
extern double *logeps;
extern double *yes;
extern double energy_shift;
extern double drho, drhoi;
extern double deps, depsi;
extern double dye, dyei;

extern double *logtemp_mode;
extern double *entr_mode;
extern double *logprss_mode;

// some vectors for selecting variables for more
// efficient interpolation

// table key 
// (alltables)
// 0 logpress 
// 1 logtemp
// 2 entropy
// 3 munu
// 4 cs2
// 5 muhat
// 6 mu_e
// 7 mu_p
// 8 mu_n
// 9 Xa
// 10 Xh
// 11 Xn
// 12 Xp
// 13 Abar
// 14 Zbar
// 15 Gamma

// alltables_mode
//  0 logenergy for temp mode
//  1 lpgenergy for entr mode
//  2 logenergy for prss mode


// frontend function declarations
void nuc_eos_C_short( double xrho,  double *xenr,  double xye, 
                      double *xtemp,  double *xent,  double *xprs,  
                      double *xcs2,   double *xmunu, double energy_shift,
                      int nrho, int neps, int nye, int nmode,
                      double *alltables, double *alltables_mode,
                      double *logrho, double *logeps, double *yes,
                      double *logtemp_mode, double *entr_mode, double *logprss_mode,
                      int keymode, int *keyerr, double rfeps );


// core function declarations
void nuc_eos_C_ReadTable( char *nuceos_table_name );

void nuc_eos_C_linterp_some( double x, double y, double z,
			                    double *output_vars, double *alltables,
                             int nx, int ny, int nz, int nvars, 
                             double *xt,double *yt, double *zt );

void nuc_eos_C_cubinterp_some( double x, double y, double z, 
					                double *output_vars, double *alltables,
                               int nx, int ny, int nz, int nvars, 
                               double *xt, double *yt, double *zt );

void find_energy( double x, double y, double z,
                  double *found_leps, double *alltables_mode,
                  int nx, int ny, int nz, int neps,
                  double *xt, double *yt, double *zt, double *logeps,
                  int keymode, int *keyerr );

void nuc_eos_C_testing();

//#endif // __NUCLEAREOS_H__
