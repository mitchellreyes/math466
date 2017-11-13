#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cilk/cilk.h>
#include <lapacke.h>
#include <string.h>
#include <strings.h>
#include<sys/time.h>
#include <sys/resource.h>
#include "matrixlib.h"
#include "tictoc.h"

double dotprod(int n, double x[n], double y[n]){
	double s = 0.0;
	for(int i = 0; i < n; i++){
		s += x[i] * y[i];
	}
	return s;
}

void gramschmidt(int m, int n, double A[m][n], double Q[m][n]){
	double vt[m];
	for(int j = 0; j < n; j++){	
		for(int i = 0; i < m; i++){
			vt[i] = A[i][j];
		}
		for(int k = 0; k < j; k++){
			double dp = dotprod(m, vt, Q[k]);
			for(int i = 0; i < m; i++){
				vt[i] -= dp * Q[k][i]; 
			}
		}
		double vtnorm = vecnorm2(m, vt);
		for(int i = 0; i < m; i++){
			Q[j][i] = vt[i] / vtnorm;
		}
	}
}


#define N 5
double A[N][N], QT[N][N], R[N][N];
int main(){
	//setrlimit(RLIMIT_STACK, 
	//&(const struct rlimit){RLIM_INFINITY, RLIM_INFINITY});
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			//A[i][j] = 2.0 * random()/RAND_MAX-1.0;
			A[i][j] = 1.0 / (i+j+1);
		}
	}
	printf("N = %d\n", N);
	printf("A = \n"); matprint(N, N, A);
	gramschmidt(N, N, A, QT);
	bzero(R, sizeof(double)*N*N);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			for(int k = 0; k < N; k++){
				R[i][j] += QT[i][k] * A[k][j];
			}
		}
	}
	printf("QT = \n"); matprint(N, N, QT);
	printf("R = \n"); matprint(N, N, R);
	double I1[N][N], I2[N][N];
	bzero(I1, sizeof(double)*N*N);
	bzero(I2, sizeof(double)*N*N);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			for(int k = 0; k < N; k++){
				I1[i][j] += QT[i][k] * QT[j][k];
				I2[i][j] += QT[k][i] * QT[k][j];
			}
		}
	}
	printf("I1 = \n"); matprint(N, N, I1);
	printf("I2 = \n"); matprint(N, N, I2);
	return 0;
}
