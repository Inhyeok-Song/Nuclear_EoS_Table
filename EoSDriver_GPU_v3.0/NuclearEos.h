//#ifndef __NUCLEAREOS_H__
//#define __NUCLEAREOS_H__
//#define DEBUG 1

#define npoints 1024*1024

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define EOS_NTABLE_MAX       20
#define EOS_NAUX_MAX         20
#define NUC_TABLE_NVAR       16     // number of variables in the EoS table lookup
#define NUC_TABLE_NPTR        8     // number of table pointers to be sent to GPU


extern int    g_nrho;
extern int    g_neps;
extern int    g_nye;
extern int    g_nmode;
extern double g_energy_shift;

extern double *g_alltables;
extern double *g_alltables_mode;
extern double *g_logrho;
extern double *g_logeps;
extern double *g_yes;
extern double *g_logtemp_mode;
extern double *g_logprss_mode;
extern double *g_entr_mode;

extern double *h_EoS_Table[EOS_NTABLE_MAX];

extern int    *d_nrho;
extern int    *d_neps;
extern int    *d_nye;
extern int    *d_nmode;
extern double *d_energy_shift;

extern double *d_EoS_Table[EOS_NTABLE_MAX];

extern double *d_Rand_Vars;


// auxiliary array indices
#define NUC_AUX_ESHIFT        0     // AuxArray_Flt: energy_shift
#define NUC_AUX_DENS2CGS      1     // AuxArray_Flt: convert density    to cgs
#define NUC_AUX_PRES2CGS      2     // AuxArray_Flt: convert pressure   to cgs
#define NUC_AUX_VSQR2CGS      3     // AuxArray_Flt: convert velocity^2 to cgs
#define NUC_AUX_PRES2CODE     4     // AuxArray_Flt: convert pressure   to code unit
#define NUC_AUX_VSQR2CODE     5     // AuxArray_Flt: convert velocity^2 to code unit

#define NUC_AUX_NRHO          0     // AuxArray_Int: nrho
#define NUC_AUX_NEPS          1     // AuxArray_Int: neps
#define NUC_AUX_NYE           2     // AuxArray_Int: nye
#define NUC_AUX_NMODE         3     // AuxArray_Int: nmode


// table indices
#define NUC_TAB_ALL           0     // alltables
#define NUC_TAB_ALL_MODE      1     // alltables_mode
#define NUC_TAB_RHO           2     // logrho
#define NUC_TAB_EPS           3     // logeps
#define NUC_TAB_YE            4     // yes
#define NUC_TAB_TEMP_MODE     5     // logtemp_mode
#define NUC_TAB_ENTR_MODE     6     // entr_mode
#define NUC_TAB_PRES_MODE     7     // logprss_mode


// EoS modes
#define NUC_MODE_ENGY         0     // energy mode
#define NUC_MODE_TEMP         1     // temperature mode
#define NUC_MODE_ENTR         2     // entropy mode
#define NUC_MODE_PRES         3     // pressure mode

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
//  1 logenergy for entr mode
//  2 logenergy for prss mode


// core function declarations
void nuc_eos_C_ReadTable( char *nuceos_table_name );
//

//void NuclearEoS_Init();

//void nuc_eos_C_linterp_some( double x, double y, double z,
//			                    double *output_vars, double *alltables,
//                             int nx, int ny, int nz, int nvars, 
//                             double *xt,double *yt, double *zt );
//
//void nuc_eos_C_cubinterp_some( double x, double y, double z, 
//					                double *output_vars, double *alltables,
//                               int nx, int ny, int nz, int nvars, 
//                               double *xt, double *yt, double *zt );
//
//void find_energy( double x, double y, double z,
//                  double *found_leps, double *alltables_mode,
//                  int nx, int ny, int nz, int neps,
//                  double *xt, double *yt, double *zt, double *logeps,
//                  int keymode, int *keyerr );
