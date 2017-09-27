#include <stdio.h>
#include <strings.h>
#include <complex.h>
#include <math.h>
#include "matrixlib.h"

#define N 4

complex A[N][N], AC[N][N];
complex X[N];

void doinit(){
	for(int k = 0; k < N; k++){
		for(int l = 0; l < N; l++){
			A[k][l] = cexp(-I*2*M_PI*k*l/N);
			AC[k][l] = cexp(I*2*M_PI*k*l/N)/N;
		}
	}
	for(int k = 0; k < N; k++){
		X[k] = 1.0 / (1+k) *(2+k);
	}
}


int main(){
	doinit();
    complex B[N], C[N];
    printf("A = \n");
    cmatprint(N, N, A);
    printf("X = \n");
    cvecprint(N, X);
    cmultAx(N, N, A, X, B);//B = A*X
    printf("B = \n");
    cvecprint(N, B);
    cmultAx(N, N, AC, B, C); //C = (AC)*B
    printf("C = \n");
    cvecprint(N, C);
    return 0;
}
