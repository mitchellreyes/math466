#include <stdio.h>
#include <strings.h>

double A[3][2] = {{1, 2}, {3, 4}, {5, 6}};
double X[2] = {7, 8};

void matprint(int numRows, int numCols, double A[numRows][numCols]){
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            printf("%g %lu ", A[i][j], &A[i][j]);
        }
        printf("\n");
    }
}

void multAx(int numRows, int numCols, double A[numRows][numCols],
	X[numCols], double B[numCols]){
        bzero(B, 8*numCols);
}

int main(){
    double B[2];
    matprint(3, 2, A);
    multAx(3, 2, A, X, B);
    return 0;
}
