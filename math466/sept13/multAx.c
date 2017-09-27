#include <stdio.h>
#include <strings.h>

double A[3][2] = {{1, 2}, {3, 4}, {5, 6}};
double X[2] = {7, 8};

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

int main(){
    double B[3];
    printf("A = \n");
    matprint(3, 2, A);
    printf("X = \n");
    vecprint(2, X);
    multAx(3, 2, A, X, B);
    printf("B = \n");
    vecprint(2, B);
    return 0;
}
