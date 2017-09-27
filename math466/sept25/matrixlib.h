#include <stdio.h>
#include <strings.h>
#include <complex.h>

void matprint(int m,int n,double A[m][n]);
void vecprint(int m,double X[m]);
void multAx(int m,int n,
    double A[m][n],double X[n],double B[m]);
void algob(int m,int n,
    double A[m][n],double X[n],double B[m]);
void cmatprint(int m,int n,complex A[m][n]);
void cvecprint(int m,complex X[m]);
void cmultAx(int m,int n,
    complex A[m][n],complex X[n],complex B[m]);
void calgob(int m,int n,
    complex A[m][n],complex X[n],complex B[m]);
