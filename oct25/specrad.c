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

double vecnorm2(int n, double x[n]){
	double r = 0;
	for(int i = 0; i < n; i++){
		double t = fabs(x[i]);
		r += t*t;
	}
	return sqrt(r);
}

double matnorm2(int n, double A[n][n]){
	double x[n], y[n], yk[n];
	/*
	double B[n][n], ;
	bzero(B, sizeof(double)*n*n);
	printf("starting\n");
	for(int k = 0; k < n; k++){
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
					B[i][j] += A[k][i]*A[k][j];
			}
		}	
	}
	printf("done\n");
	*/
	for(int i = 0; i < n; i++){
		x[i] = 2.0 * random()/RAND_MAX-1.0;
	}
	//multAx(n, n, B, x, y);
	double z[n];
	multAx(n, n, A, x, z);
	multATx(n, n, A, z, y);
	double t, tk,qold=0, qk;
	//yk = B^i x
	for(int i = 2; i < 10000; i++){
		//multAx(n, n, B, y, yk);
		multAx(n, n, A, x, z);
		multATx(n, n, A, z, yk);
		t = vecnorm2(n, y);
		tk = vecnorm2(n, yk);
		qk = tk/t;
		if(fabs(qk-qold)<=qk*5e-14){
			return sqrt(qk);		
		}
		printf("%g, %g tk/t = %.15e\n",t, tk, tk/t);
		for(int j = 0; j < n; j++){
			y[j] = yk[j]/tk;
		}
		qold = qk;
	}
	printf("loop ended :(\n");
	return(sqrt(tk/t));
}

double vecnorm1(int n, double x[n]){
	double r = 0;
	for(int i = 0; i < n; i++){
		r += fabs(x[i]);
	}
	return r;
}

double matnorm1(int n, double A[n][n]){
	double r = 0;
	for(int j = 0; j < n; j++){
		double s = 0;
		for(int i = 0; i < n; i++){
			s += fabs(A[i][j]);
		}
		if(s > r) r = s;
	}
	return r;
}


#define N 2000
double A[N][N];
double B[N], X[N], XX[N];
double *P[N];
int IP[N];
int main(){
	//setrlimit(RLIMIT_STACK, 
	//&(const struct rlimit){RLIM_INFINITY, RLIM_INFINITY});
	printf("N = %d\n", N);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			A[i][j] = 2.0 * random()/RAND_MAX-1.0;
		}
		X[i] = 1;
	}
	printf("|A|_1 = %g\n", matnorm1(N, A));
	printf("|A|_2 = %g\n", matnorm2(N, A));
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
	printf("Speed of %g double-percision GFLOPS.\n", 
	2.0 * (N / 3.0 + 1) * N * N/t/1e9);
	return 0;
}
