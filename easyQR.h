#ifndef EASY_QR
#define EASY_QR

#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <omp.h>
#include <string.h>

#ifndef _OPENMP
	#include <time.h>
#endif	

#define BUFF_SIZE 512

///////////////////////////////////////// public structs ////////////////////////////////////////////////////////////////////////

typedef enum SHIFTSTRATEGY{ZERO,LASTDIAG,WILKINSON}shiftStrategy;			/**< Shift strategy to apply for the QR algorithm. See https://de.wikipedia.org/wiki/QR-Algorithmus for details. */

//////////////////////////////////////// public functions ///////////////////////////////////////////////////////////////////////

void testQR(int n,int threads);
void QRsetVerboseLevel(int level);
void QRiterations(void* A,double* eigenValues,double** eigenVectors,int eig_num,double tol,int max_iter,shiftStrategy st);
void getEigen2D(double* A,double* eig1_real,double* eig2_real,double* eig_im);
void getRealEigenVectors2D(double* denseH,double* eig,double* denseQ,double tol);
double denseRayleighQuotient(double* A,double* x,int dim);

#endif
