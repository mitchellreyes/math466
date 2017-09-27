#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <complex.h>
#include <math.h>
#include <cilk/cilk.h>
#include "tictoc.h"
#include "matrixlib.h"
#define N 262144
complex X[N];
void doinit(){
	for(int k=0;k<N;k++){
		X[k]=1.0/(1+k)*(2+k)+I*1.0/(N-k);
	}
}

//only works when n is a power of 2
void fft(int n, int s, complex x[n], complex b[n]){
	if(n==1){
		b[0] = x[0];
		return;
	}
	if(n%2){
		printf("Error: n was not a power of 2!!!! >:(\n");
		exit(1);
	}
	int K = n/2;
	fft(K, 2*s, &x[s], &b[K]);
	fft(K, 2*s, &x[0], &b[0]);
	for(int k = 0; k < K; k++){
		complex w = cexp(-I*2*M_PI*k/n);
		complex be = b[k], bo=b[k+K];
		b[k] = be + w*bo;
		b[k+K]= be - w*bo;
	}
}

//only works when n is a power of 2
void pfft(int n, int s, complex x[n], complex b[n]){
	if(s > 8){
		fft(n, s, x, b);
		return;
	}
	if(n==1){
		b[0] = x[0];
		return;
	}
	if(n%2){
		printf("Error: n was not a power of 2!!!! >:(\n");
		exit(1);
	}
	int K = n/2;
	cilk_spawn pfft(K, 2*s, &x[s], &b[K]);
	pfft(K, 2*s, &x[0], &b[0]);
	cilk_sync;
	cilk_for(int k = 0; k < K; k++){
		complex w = cexp(-I*2*M_PI*k/n);
		complex be = b[k], bo=b[k+K];
		b[k] = be + w*bo;
		b[k+K]= be - w*bo;
	}
}

int main(){
	printf("N=%d\n", N);
	tic();
	doinit();
	double t=toc();
	printf("Initialization time was %g seconds.\n",t);
    complex C[N];
	//printf("X=\n"); cvecprint(N,X);
	//printf("B=\n"); cvecprint(N,B);
	tic();
	fft(N, 1, X, C);
	t = toc();
	printf("Elapsed time was %g seconds.\n",t);
	
	tic();
	pfft(N, 1, X, C);
	t = toc();
	printf("Parallel time was %g seconds.\n",t);
	//printf("C=\n"); cvecprint(N,C);
    return 0;
}
