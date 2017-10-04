#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixlib.h"

void lufact(int n,double A[n][n]){
	for(int i=0;i<n-1;i++){
		for(int j=i+1;j<n;j++){
			double alpha=A[j][i]/A[i][i];
			for(int k=i+1;k<n;k++){
				A[j][k]-=alpha*A[i][k];
			}
			A[j][i]=alpha;
		}
	}
}
#define N 3
double A[3][3] = {{1,2,-1}, {2,5,-3}, {-1,2,0}};
int main(){
	printf("A=\n"); matprint(N,N,A);
	lufact(N,A);
	printf("LU=\n"); matprint(N,N,A);
	return 0;
}
