#include <stdio.h>
#include <strings.h>
#include <complex.h>
void matprint(int numRows, int numCols, double A[numRows][numCols]);

//this is for rows
void multAx(int numRows, int numCols, double A[numRows][numCols], double X[numCols], double B[numRows]);

void vecprint(int length, double X[length]);

//this is for columns
void algoB(int numRows, int numCols, double A[numRows][numCols], double X[numCols], double B[numRows]);

void cmatprint(int numRows, int numCols, complex A[numRows][numCols]);

//this is for rows
void cmultAx(int numRows, int numCols, complex A[numRows][numCols], complex X[numCols], complex B[numRows]);

void cvecprint(int length, complex X[length]);

//this is for columns
void calgoB(int numRows, int numCols, complex A[numRows][numCols], complex X[numCols], complex B[numRows]);










