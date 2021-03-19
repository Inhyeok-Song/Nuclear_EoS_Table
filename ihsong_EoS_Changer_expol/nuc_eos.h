#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 16
//#define DEBUG 1

extern int nrho;
extern int neps;
extern int nye;

extern int nrho_mode;
extern int nmode;
extern int nye_mode;

extern double *alltables;
extern double *alltables_mode;
extern double *logrho;
extern double *logeps;
extern double *yes;
extern double energy_shift;

extern double *logrho_mode;
extern double *logtemp_mode;
extern double *entr_mode;
extern double *logprss_mode;
extern double *yes_mode;


// outvars

extern double xrho, xtemp, xeps, xye;
extern double xent, xprs, xcs2;
extern double xmunu, xmuhat, xmu_e, xmu_p, xmu_n;
extern double xXa, xXh, xXn, xXp;
extern double xAbar, xZbar;
extern double xgamma;


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
extern int ivs_short[16];

// core function declarations

void nuc_eos_C_ReadTable(char *nuceos_table_name);


void write_eos_table(char *nuceos_eps_table_name);