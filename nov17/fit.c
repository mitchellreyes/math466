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

double dotprod(int n, double x[n], double y[n]){
	double s = 0.0;
	for(int i = 0; i < n; i++){
		s += x[i] * y[i];
	}
	return s;
}

void gramschmidt(int m, int n, double A[m][n], double Q[m][n]){
	double vt[m];
	for(int j = 0; j < n; j++){	
		for(int i = 0; i < m; i++){
			vt[i] = A[i][j];
		}
		for(int k = 0; k < j; k++){
			double dp = dotprod(m, vt, Q[k]);
			for(int i = 0; i < m; i++){
				vt[i] -= dp * Q[k][i]; 
			}
		}
		double vtnorm = vecnorm2(m, vt);
		for(int i = 0; i < m; i++){
			Q[j][i] = vt[i] / vtnorm;
		}
	}
}

void hheliminate(int m, int n, double v[m], double A[m][n]){
	for(int j = 0; j < n; j++){
		double a[m];
		for(int i = 0; i < m; i++){
			a[i] = A[i][j];
		}
		double vdota2 = 2 * dotprod(m, v, a);
		for(int i = 0; i < m; i++){
			a[i] -= v[i]*vdota2;
		}
		for(int i = 0; i < m; i++){
			A[i][j] = a[i];
		}
	}
}

void householder(int m, int n, double A[m][n], 
double Q[m][n], double R[m][n]){
	double vt[m];
	double H[n][m];
	memcpy(R, A, sizeof(double)*m*n);
	for(int j = 0; j < n; j++){
		for(int i = 0; i < m; i++){
			if(i < j){
				vt[i] = 0;
			}
			else
			{
				vt[i] = R[i][j];
			}
		}
		double c = vecnorm2(m, vt);
		if(vt[j] < 0){
			vt[j] -= c;
		}else{
			vt[j] += c;
		}
	
		double vnorm = vecnorm2(m, vt);
		for(int i = 0; i < m; i++){
			vt[i] /= vnorm;
		}
		hheliminate(m, n, vt, R);
		memcpy(H[j], vt, sizeof(double)*m);
		//printf("R = \n"); matprint(m, n, R);
	}
	bzero(Q, sizeof(double)*m*n);
	for(int j = 0; j < n; j++){
		Q[j][j] = 1;
	}
	for(int j = n-1; j >= 0; j--){
		hheliminate(m, n, H[j], Q);
	}
	for(int i = 1; i < n; i++){
		for(int j = 0; j < i; j++){
			R[i][j] = 0;
		}
	}
	
}
/*
#define N 4

double phi0(double x){ return 1.0; }
double phi1(double x){ return x; }
double phi2(double x){ return x*x; }
double phi3(double x){ return x*x*x; }
*/
#define N 4

double phi0(double x){ return 1.0; }
double phi1(double x){ return cos(x); }
double phi2(double x){ return cos(sqrt(2.0)*x); }
double phi3(double x){ return cos(sqrt(3.0)*x); }

typedef double func(double x);
func *phi[N] = {phi0, phi1, phi2, phi3};

double phi_easy(double base, int power){
	return pow(base, power);
}	

#define N 4


int main(){
	setrlimit(RLIMIT_STACK, 
	&(const struct rlimit){RLIM_INFINITY, RLIM_INFINITY});
	int m;
	//double A[M][N], QT[N][M], R[M][N], Q[M][N];
	FILE *fp = fopen("file05b.dat", "r");
	fscanf(fp, "%d", &m);
	printf("m, N = %d, %d\n", m, N);
	double X[m], Y[m], A[m][N];
	for(int i = 0; i < m; i++){
		fscanf(fp, "%lf %lf", &X[i], &Y[i]);
	}
	for(int i = 0; i < m; i++){
		printf("%lf %lf\n", X[i], Y[i]);
	}
	for(int i = 0; i < m; i++){
		for(int j = 0; j < N; j++){
			A[i][j] = phi[j](X[i]);
			//A[i][j] = phi_easy(X[i], j);
		}
	}
	double R[m][N], Q[m][N];
	householder(m, N, A, Q, R);
	printf("R = \n"); matprint(N, N, R);
	double b[N];
	bzero(b, sizeof(double)*N);
	for(int i = 0; i < m; i++){
		for(int j = 0; j < N; j++){
			b[j] += Q[i][j]*Y[i];
		}
	}
	double c[N];
	backsub(N, R, c, b);
	printf("c = \n"); vecprint(N, c);
	return 0;
	
}
