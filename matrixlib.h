#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <strings.h>
#include <cilk/cilk.h>

extern void matprint(int m,int n,double A[m][n]);
extern void vecprint(int m,double X[m]);
extern void multAx(int m,int n,double A[m][n],double X[n],double B[m]);
extern void multATx(int m,int n,double A[m][n],double X[m],double B[n]);
extern void LUfact(int n,double A[n][n]);
extern void LUsolve(int n,double LU[n][n],double x[n],double b[n]);

// Using these functions
extern void PLUsolve(int n,double LU[n][n],double *P[n],double x[n],double b[n]);
extern void PLUfact(int n,double A[n][n],double *P[n]);
//

extern double vecnorm1(int n,double x[n]);
extern double matnorm1(int n,double A[n][n]);
extern double vecnorm2(int n,double x[n]);
