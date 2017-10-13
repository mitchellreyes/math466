#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cilk/cilk.h>
#include <lapacke.h>
#include <string.h>
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

#define N 2000
double A[N][N];
double B[N], X[N], XX[N];
double *P[N];
int IP[N];
int main(){
	printf("N = %d\n", N);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			A[i][j] = 2.0 * random()/RAND_MAX-1.0;
		}
		X[i] = 1;
	}
	multAx(N, N, A, X, B); //B = AX
	tic();
	//plufact(N,A, P); //now A is overwritten with LU
	//plusolve(N, A, P, XX, B);
	memcpy(XX, B, sizeof(double)*N);
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, 1, &A[0][0], N, IP, XX, 1);
	
	double t = toc();
	printf("it took %g seconds!\n", t);
	//printf("XX=\n"); vecprint(N, XX);
	double r = 0;
	for(int i = 0; i < N; i++){
		double dx = XX[i] - 1;
		r += dx*dx;
	}
	r = sqrt(r);
	printf("The error is %g.\n", r);
	printf("Speed of %g double-percision FLOPS.\n", 
	2.0 * (N / 3.0 + 1) * N * N/t/1e6/1000);
	return 0;
}
