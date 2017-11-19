#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"
#include "tictoc.h"

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
		for(int j=i+1;j<n;j++){
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

#define N 3
double A[3][3] = {{1,2,-1}, {2,5,-3}, {-1,2,0}};
double B[3], X[3] = {1, 2, 3}, XX[3];
double *P[3];
int main(){
	//printf("A=\n"); matprint(N,N,A);
	tic();
	multAx(N, N, A, X, B); //B = AX
	//printf("B=\n"); vecprint(N, B);
	plufact(N,A, P); //now A is overwritten with LU
	printf("LU=\n"); matprint(N,N,A);
	plusolve(N, A, P, XX, B);
	double t = toc();
	printf("it took %g seconds!\n", t);
	printf("XX=\n"); vecprint(N, XX);
	return 0;
}
