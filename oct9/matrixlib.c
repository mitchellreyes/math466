#include <stdio.h>
#include <strings.h>
#include <complex.h>
//#include <cilk/cilk.h>
#include "matrixlib.h"

void matprint(int m,int n,double A[m][n]){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("%-15g\t",A[i][j]);
        }
        printf("\n");
    }
}
void vecprint(int m,double X[m]){
    matprint(m,1,(double (*)[1])X);
}
void multAx(int m,int n,
    double A[m][n],double X[n],double B[m]){
    bzero(B,sizeof(double)*m);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]+=A[i][j]*X[j];
        }
    }
}
void algob(int m,int n,
    double A[m][n],double X[n],double B[m]){
    bzero(B,8*m);
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++){
            B[i]+=A[i][j]*X[j];
        }
    }
}
void cmatprint(int m,int n,complex A[m][n]){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("(%g %g) ",A[i][j]);
        }
        printf("\n");
    }
}
void cvecprint(int m,complex X[m]){
    cmatprint(m,1,(complex (*)[1])X);
}
void cmultAx(int m,int n,
    complex A[m][n],complex X[n],complex B[m]){
    bzero(B,sizeof(complex)*m);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]+=A[i][j]*X[j];
        }
    }
}
void calgob(int m,int n,
    complex A[m][n],complex X[n],complex B[m]){
    bzero(B,8*m);
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++){
            B[i]+=A[i][j]*X[j];
        }
    }
}
