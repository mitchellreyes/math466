#include <stdio.h>
#include <strings.h>
#include <complex.h>
#include <stdlib.h>
void matprint(int m,int n,double A[m][n]);
void vecprint(int m,double X[m]);
void multAx(int m,int n,
    double A[m][n],double X[n],double B[m]);
void multATx(int m,int n,
    double A[n][m],double X[n],double B[m]);
void algob(int m,int n,
    double A[m][n],double X[n],double B[m]);
void cmatprint(int m,int n,complex A[m][n]);
void cvecprint(int m,complex X[m]);
void cmultAx(int m,int n,
    complex A[m][n],complex X[n],complex B[m]);
void calgob(int m,int n,
    complex A[m][n],complex X[n],complex B[m]);
void lufact(int n,double A[n][n]);
void backsub(int n, double U[n][n], double x[n], double b[n]);
void lusolve(int n, double LU[n][n], double x[n], double b[n]);
void plufact(int n,double A[n][n], double *P[n]);
void plusolve(int n, double LU[n][n], double *P[n], double x[n], double b[n]);
double vecnorm2(int n, double x[n]);
double matnorm2(int n, double A[n][n]);
double vecnorm1(int n, double x[n]);
double matnorm1(int n, double A[n][n]);
