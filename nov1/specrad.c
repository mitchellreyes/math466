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
