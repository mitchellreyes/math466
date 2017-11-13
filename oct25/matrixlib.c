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
void multAx(int m,int n,
    double A[m][n],double X[n],double B[m]){
    bzero(B,sizeof(double)*m);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]+=A[i][j]*X[j];
        }
    }
}

void multATx(int m,int n,
    double A[n][m],double X[n],double B[m]){
    bzero(B,sizeof(double)*n);
    for(int j=0;j<m;j++){
		for(int i=0;i<n;i++)
		{
		  B[i]+=A[j][i]*X[j];
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

void lufact(int n,double A[n][n]){
	for(int i=0;i<n-1;i++){
		for(int j=i+1;j<n;j++){
			double alpha=A[j][i]/A[i][i];
			for(int k=i+1;k<n;k++){
				A[j][k]-=alpha*A[i][k];
			}
			A[j][i]=alpha;
		}
	}
}

void lusolve(int n, double LU[n][n], double x[n], double b[n]){ 
	double y[n];
	for(int i = 0; i < n; i++){ //Ly = b
		y[i] = b[i];
		for(int j = 0; j < i; j++){
			y[i] -= LU[i][j]* y[j];
		}
	}
	//Ux = y
	for(int i = n-1; i >= 0; i--){
		x[i] = y[i];
		for(int j = i + 1; j < n; j++){
			x[i] -= LU[i][j] * x[j];
		}
		x[i]/=LU[i][i];
	}
}

void plufact(int n,double A[n][n], double *P[n]){
	for(int i = 0; i < n; i++){
		P[i] = &A[i][0];
	}
	for(int i=0;i<n-1;i++){
		//switch rows
		for(int j = i + 1; j < n; j++){
			if(fabs(P[i][i]) < fabs(P[j][i])){
				//swap
				double *t = P[i];
				P[i] = P[j];
				P[j] = t;
			}
		}
		cilk_for(int j=i+1;j<n;j++){
			double alpha=P[j][i]/P[i][i];
			for(int k=i+1;k<n;k++){
				P[j][k]-=alpha*P[i][k];
			}
			P[j][i]=alpha;
		}
	}
}

void plusolve(int n, double LU[n][n], double *P[n], double x[n], double b[n]){
	double y[n];
	for(int i = 0; i < n; i++){ //Ly = b
		y[i] = b[(P[i] - &LU[0][0])/n];
		for(int j = 0; j < i; j++){
			y[i] -= P[i][j]* y[j];
		}
	}
	//Ux = y
	for(int i = n-1; i >= 0; i--){
		x[i] = y[i];
		for(int j = i + 1; j < n; j++){
			x[i] -= P[i][j] * x[j];
		}
		x[i]/=P[i][i];
	}
}

