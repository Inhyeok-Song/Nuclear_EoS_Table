//#ifndef __NUCLEAREOS_H__
//#define __NUCLEAREOS_H__

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//#define DEBUG 1

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


// frontend function declaration
//__device__ static
//void nuc_eos_C_short( double xrho,  double *xenr,  double xye, 
//                      double *xtemp,  double *xent,  double *xprs,  
//                      double *xcs2,   double *xmunu, double energy_shift,
//                      int nrho, int neps, int nye, int nmode,
//                      double *alltables, double *alltables_mode,
//                      double *logrho, double *logeps, double *yes,
//                      double *logtemp_mode, double *entr_mode, double *logprss_mode,
//                      int keymode, int *keyerr, double rfeps );

//void nuc_eos_C_testing();

//void PassRandPoints2GPU();
//void PassNuclearEoSTable2GPU();
//void MemFree_NuclearEoS();


#define CUDA_CHECK_ERROR( Call )   CUDA_Check_Error( Call, __FILE__, __LINE__, __FUNCTION__ )

inline void CUDA_Check_Error( cudaError Return, const char *File, const int Line, const char *Func )
{
   //if ( Return != cudaSuccess ) 
   //   Aux_Error( ERROR_INFO, "CUDA ERROR : %s !!\n", cudaGetErrorString( Return ) );
}


 // __NUCLEAREOS_H__
