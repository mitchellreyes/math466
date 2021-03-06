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

void hheliminate(int m, int n, double v[m], double A[m][n]){
	for(int j = 0; j < n; j++){
		double a[m];
		for(int i = 0; i < m; i++){
			a[i] = A[i][j];
		}
		double vdota2 = 2 * dotprod(m, v, a);
		for(int i = 0; i < m; i++){
			a[i] -= v[i]*vdota2;
		}
		for(int i = 0; i < m; i++){
			A[i][j] = a[i];
		}
	}
}

void householder(int m, int n, double A[m][n], 
double Q[m][n], double R[m][n]){
	double vt[m];
	double H[n][m];
	memcpy(R, A, sizeof(double)*m*n);
	for(int j = 0; j < n; j++){
		for(int i = 0; i < m; i++){
			if(i < j){
				vt[i] = 0;
			}
			else
			{
				vt[i] = R[i][j];
			}
		}
		double c = vecnorm2(m, vt);
		if(vt[j] < 0){
			vt[j] -= c;
		}else{
			vt[j] += c;
		}
	
		double vnorm = vecnorm2(m, vt);
		for(int i = 0; i < m; i++){
			vt[i] /= vnorm;
		}
		hheliminate(m, n, vt, R);
		memcpy(H[j], vt, sizeof(double)*m);
		//printf("R = \n"); matprint(m, n, R);
	}
	bzero(Q, sizeof(double)*m*n);
	for(int j = 0; j < n; j++){
		Q[j][j] = 1;
	}
	for(int j = n-1; j >= 0; j--){
		hheliminate(m, n, H[j], Q);
	}
	for(int i = 1; i < n; i++){
		for(int j = 0; j < i; j++){
			R[i][j] = 0;
		}
	}
	
}

#define N 4
#define M 6
double A[M][N], QT[N][M], R[M][N], Q[M][N];
int main(){
	//setrlimit(RLIMIT_STACK, 
	//&(const struct rlimit){RLIM_INFINITY, RLIM_INFINITY});
	for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++){
			//A[i][j] = 2.0 * random()/RAND_MAX-1.0;
			A[i][j] = 1.0 / (i+j+1);
		}
	}
	printf("N = %d, M = %d\n", N, M);
	printf("A = \n"); matprint(M, N, A);
	householder(M, N, A, Q, R);
	printf("R = \n"); matprint(N, N, R);
	printf("Q = \n"); matprint(M, N, Q);
	double QR[M][N];
	bzero(QR, sizeof(double)*M*N);
	for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++){
			for(int k = 0; k < N; k++){
				QR[i][j] += Q[i][k] * R[k][j];
			}
		}
	}
	printf("QR = \n"); matprint(M, N, QR);
	double I1[N][N];
	bzero(I1, sizeof(double)*N*N);

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			for(int k = 0; k < M; k++){
				I1[i][j] += Q[k][i] * Q[k][j];
			}
		}
	}
	printf("I1 = \n"); matprint(N, N, I1);
	return 0;
	
}
