#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

void lufact(int size, double A[size][size]){
  //for each col
  for(int i = 0; i < size-1; i++){
    //for each row
    for(int j = i + 1; j < size; j++){
      double alpha = A[j][i]/A[i][i];
      for(int k = i + 1; k < size; k++){
        A[j][k] -= alpha*A[i][k];
      }
      A[j][i] = alpha;
    }
  }
}
#define N 3
double A[N][N] = {1, 2, -1, 2, 5, -3, -1, 2, 0};
int main(){
  printf("A = \n");
  matprint(N, N, A);
  lufact(N, A);
  printf("LU = \n");
  matprint(N, N, A);
  return 0;
}
