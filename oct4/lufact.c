#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

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

#define N 3
double A[3][3] = {{1,2,-1}, {2,5,-3}, {-1,2,0}};
double B[3], X[3] = {1, 2, 3}, XX[3];
int main(){
	printf("A=\n"); matprint(N,N,A);
	multAx(N, N, A, X, B); //B = AX
	printf("B=\n"); vecprint(N, B);
	lufact(N,A); //now A is overwritten with LU
	printf("LU=\n"); matprint(N,N,A);
	lusolve(N, A, XX, B);
	printf("XX=\n"); vecprint(N, XX);
	return 0;
}
