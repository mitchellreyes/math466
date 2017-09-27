#include <stdio.h>
#include <strings.h>
#include <complex.h>
#include "matrixlib.h"
void matprint(int numRows, int numCols, double A[numRows][numCols]){
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            printf("%g %lu ", A[i][j], (long unsigned)&A[i][j]);
        }
        printf("\n");
    }
}

//this is for rows
void multAx(int numRows, int numCols, double A[numRows][numCols], double X[numCols], double B[numRows]){
        bzero(B, 8*numCols);
        for(int i = 0; i < numRows; i++){
            for(int j = 0; j < numCols; j++){
                B[i] += A[i][j]*X[j];
            }
        }
}

void vecprint(int length, double X[length]){
	matprint(length, 1, (double (*)[1])X);
}

//this is for columns
void algoB(int numRows, int numCols, double A[numRows][numCols], double X[numCols], double B[numRows]){
        bzero(B, sizeof(double)*numCols);
        for(int j = 0; j < numCols; j++){
            for(int i = 0; i < numRows; i++){
                B[i] += A[i][j]*X[j];
            }
        }
}

void cmatprint(int numRows, int numCols, complex A[numRows][numCols]){
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            printf("(%g %g) ", A[i][j]);
        }
        printf("\n");
    }
}

//this is for rows
void cmultAx(int numRows, int numCols, complex A[numRows][numCols], complex X[numCols], complex B[numRows]){
        bzero(B, sizeof(complex)*numCols);
        for(int i = 0; i < numRows; i++){
            for(int j = 0; j < numCols; j++){
                B[i] += A[i][j]*X[j];
            }
        }
}

void cvecprint(int length, complex X[length]){
	cmatprint(length, 1, (complex (*)[1])X);
}

//this is for columns
void calgoB(int numRows, int numCols, complex A[numRows][numCols], complex X[numCols], complex B[numRows]){
        bzero(B, sizeof(complex)*numCols);
        for(int j = 0; j < numCols; j++){
            for(int i = 0; i < numRows; i++){
                B[i] += A[i][j]*X[j];
            }
        }
}

