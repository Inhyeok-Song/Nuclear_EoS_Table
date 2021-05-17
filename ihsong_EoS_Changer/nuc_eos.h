#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 19
//#define DEBUG 1

extern int nrho;
extern int nrho2;
extern int ntemp;
extern int neps2;
extern int nye;
extern int nye2;

extern int nrho_mode;
extern int ntemp_mode;
extern int nentr_mode;
extern int nprss_mode;
extern int nye_mode;
extern int nmode;

extern double *alltables;
extern double *logrho;
extern double *logrho2;
extern double *logtemp;
extern double *logeps2;
extern double *yes;
extern double *yes2;
extern double energy_shift;
extern double energy_shift2;
extern double dtemp, dtempi;
extern double drho, drhoi;
extern double dye, dyei;

extern double *logrho_mode;
extern double *logtemp_mode;
extern double *entr_mode;
extern double *logprss_mode;
extern double *yes_mode;

// min and max values

extern double eos_rhomax, eos_rhomin;
extern double eos_tempmin, eos_tempmax;
extern double eos_yemin, eos_yemax;
extern double eos_epsmin, eos_epsmax;


// outvars

extern double xrho, xtemp, xeps, xye;
extern double xent, xprs, xcs2;
extern double xmunu, xmuhat, xmu_e, xmu_p, xmu_n;
extern double xXa, xXh, xXn, xXp;
extern double xAbar, xZbar;
extern double xgamma;
extern double dpde;


// table key
// 0 logpress 
// 1 logenergy
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

// some vectors for selecting variables for more
// efficient interpolation
extern int ivs_short[19];

// frontend function declarations
void nuc_eos_C_short(double xrho, double *xtemp, double xye,
         double *xenr, double *xprs, double *xent,
         double *xcs2, double *xmunu, double *xmuhat,
         double *xmu_e, double *xmu_p, double *xmu_n,
         double *xXa, double *xXh, double *xXn,
	       double *xXp, double *xAbar, double *xZbar,
	       double *xgamma, double *dpde, int keytemp, int *keyerr,
		     double rfeps); 

// core function declarations

void nuc_eos_C_ReadTable(char *nuceos_table_name);

void nuc_eos_C_linterp_some(double x, double y, double z,
			    double *f, double *ft, 
			    int *ivs,
			    int nx, int ny, int nz, int nvars,
			    double *xt,double *yt, double *zt);

void nuc_eos_C_cubinterp_some(double x, double y, double z,
			    double *f, double *ft, 
			    int *ivs,
			    int nx, int ny, int nz, int nvars,
			    double *xt,double *yt, double *zt);

void nuc_eos_C_linterp_for_temp(double x, double y, double z,
				double *f, double *ft, 
				int nx, int ny, int nz, 
				double *xt, double *yt, double *zt,
				double *linterp_for_temp);

void nuc_eos_C_linterp_for_entr(double x, double y, double z,
        double *f, double *ft,
        int nx, int ny, int nz,
        double *xt,double *yt, double *zt,
        double *dlepsdlt);

void nuc_eos_C_linterp_for_prss(double x, double y, double z,
        double *f, double *ft,
        int nx, int ny, int nz,
        double *xt,double *yt, double *zt,
        double *dlepsdlt);


void nuc_eos_C_findtemp(double lr, double lt0, double ye, 
			double lepsin, double prec, double *ltout, int *keyerrt);

void nuc_eos_C_findtemp_entropy(double lr, double lt0, double ye,
      double xs, double *ltout, double prec, int *keyerrt);

void nuc_eos_C_findtemp_pressure(double lr, double lt0, double ye,
      double lprsin, double *ltout, double prec, int *keyerrt);

void set_bins(int rho_bs, int eps_bs, int ye_bs, double eos_rhomin, double eos_rhomax,
		  double eos_epsmin, double eos_epsmax, double eos_yemin, double eos_yemax);

void set_bins2(int rho_bs, int temp_bs, int entr_bs, int prss_bs, int ye_bs,
               double eos_rhomin,  double eos_rhomax,
               double eos_tempmin, double eos_tempmax,
               double eos_entrmin, double eos_entrmax,
               double eos_prssmode, double eos_prssmax,
               double eos_yemin,   double eos_yemax);

void write_eos_table(char *nuceos_eps_table_name);
