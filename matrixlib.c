#include <stdio.h>
#include <strings.h>
#include <complex.h>
#include <cilk/cilk.h>
#include <math.h>
#include "matrixlib.h"

void matprint(int m,int n,double A[m][n]){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("%g ",A[i][j]);
        }
        printf("\n");
    }
}
void vecprint(int m,double X[m]){
    matprint(m,1,(double (*)[1])X);
}

void multAx(int m,int n,double A[m][n],
    double X[n],double B[m]){
    bzero(B,m*sizeof(double));
    cilk_for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]+=A[i][j]*X[j];
        }
    }
}
void multATx(int m,int n,double A[m][n],
    double X[m],double B[n]){
    bzero(B,n*sizeof(double));
    for(int j=0;j<m;j++){
        cilk_for(int i=0;i<n;i++){
            B[i]+=A[j][i]*X[j];
        }
    }
}
void LUfact(int n,double A[n][n]){
    for(int j=0;j<n-1;j++){
        for(int i=j+1;i<n;i++){
            double alpha=A[i][j]/A[j][j];
            for(int k=j+1;k<n;k++){
                A[i][k]-=alpha*A[j][k];
            }
            A[i][j]=alpha;
        }
    }
}
void LUsolve(int n,double LU[n][n],
    double x[n],double b[n]){
    double y[n];
    for(int i=0;i<n;i++){ // Solve Ly=b
        y[i]=b[i];
        for(int j=0;j<i;j++){
            y[i]-=LU[i][j]*y[j];
        }
    }
    for(int i=n-1;i>=0;i--){ // Solve Ux=y
        x[i]=y[i];
        for(int j=i+1;j<n;j++){
            x[i]-=LU[i][j]*x[j];
        }
        x[i]/=LU[i][i];
    }
}

void PLUsolve(int n,double LU[n][n],double *P[n], double x[n],double b[n]){
    double y[n];
    for(int i=0;i<n;i++){ // Solve Ly=b
        y[i]=b[(P[i]-&LU[0][0])/n];
        for(int j=0;j<i;j++){
            y[i]-=P[i][j]*y[j];
        }
    }
    for(int i=n-1;i>=0;i--){ // Solve Ux=y
        x[i]=y[i];
        for(int j=i+1;j<n;j++){
            x[i]-=P[i][j]*x[j];
        }
        x[i]/=P[i][i];
    }
}

void PLUfact(int n,double A[n][n],double *P[n]){
    for(int i=0;i<n;i++){
        P[i]=&A[i][0];
    }
    for(int j=0;j<n-1;j++){
        for(int i=j+1;i<n;i++){
            if(fabs(P[j][j])<fabs(P[i][j])){
                double *t=P[j]; 
                P[j]=P[i];
                P[i]=t;
            }
        }
        cilk_for(int i=j+1;i<n;i++){
            double alpha=P[i][j]/P[j][j];
            for(int k=j+1;k<n;k++){
                P[i][k]-=alpha*P[j][k];
            }
            P[i][j]=alpha;
        }
    }
}

double vecnorm1(int n,double x[n]){
    double r=0;
    for(int i=0;i<n;i++){
        r+=fabs(x[i]);
    }
    return r;
}
double matnorm1(int n,double A[n][n]){
    double r=0;
    for(int j=0;j<n;j++){
        double s=0;
        for(int i=0;i<n;i++){
            s+=fabs(A[i][j]);
        }
        if(s>r) r=s;
    }
    return r;double B[n][n];
}
double vecnorm2(int n,double x[n]){
    double r=0;
    for(int i=0;i<n;i++){
        r+=x[i]*x[i];
    }
    return sqrt(r);
}

