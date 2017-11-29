#include <stdio.h>
#include "matrixlib.h"

#define N 5

int main(){
	complex A[5][5] = {
    { -4, -3, -7,  5,  7 },
    {  0,  2, -7,  6, -6 },
    { -2, -6, -2,  4, -1 },
    { -1, -7, -5, -7, -1 },
    { -1, -2, -2, -8,  1 }};
	complex alpha = 0.5 + 4i;
	complex test = shiftinvpower(N, A, alpha);
	printf("The number is %f + i%f\n", creal(test), cimag(test)); 
	return 0;
}
