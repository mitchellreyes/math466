#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <strings.h>
#include <cilk/cilk.h>
#include <lapacke.h>
#include "matrixlib.h"

double dotprod(int n,double x[n],double y[n]){
    double s=0;
    for(int i=0;i<n;i++){
        s+=x[i]*y[i];
    }
    return s;
}

double phi0(double x){ return 1.0; }
double phi1(double x){ return x; }
double phi2(double x){ return x*x; }
double phi3(double x){ return x*x*x; }

typedef double func(double);
func *phi[4] = { phi0, phi1, phi2, phi3 };

void gramschmidt(int m,int n,double A[m][n],
    double QT[n][m]){
    double vtilde[m];
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++){
            vtilde[i]=A[i][j];
        }
        for(int k=0;k<j;k++){
            double adotv=dotprod(m,vtilde,QT[k]);
            for(int l=0;l<m;l++){
                vtilde[l]-=adotv*QT[k][l];
            }
        }
        double vnorm=vecnorm2(m,vtilde);
        for(int i=0;i<m;i++){
            QT[j][i]=vtilde[i]/vnorm;
        }
    }
}

void hheliminate(int m,int n,double w[m],
    double A[m][n]){
    for(int j=0;j<n;j++){
        double at[m];
        for(int i=0;i<m;i++) at[i]=A[i][j];
        double wdpat=dotprod(m,w,at);
        for(int i=0;i<m;i++) at[i]-=2*w[i]*wdpat;
        for(int i=0;i<m;i++) A[i][j]=at[i];
    }
}
void householder(int m,int n,double A[m][n],
    double R[m][n], double Q[m][n]){
    double Hs[n][m];
    memcpy(R,A,sizeof(double)*m*n);
    for(int j=0;j<n;j++){
        double vtilde[m];
        for(int i=0;i<m;i++){
            if(i<j) vtilde[i]=0;
            else vtilde[i]=R[i][j];
        }
        double c=vecnorm2(m,vtilde);
        if(vtilde[j]<0) vtilde[j]-=c;
        else vtilde[j]+=c;
        double vnorm=vecnorm2(m,vtilde);
        for(int i=0;i<m;i++){
            vtilde[i]/=vnorm;
        }
        hheliminate(m,n,vtilde,R);
 //       printf("R=\n"); matprint(m,n,R);
        memcpy(Hs[j],vtilde,sizeof(double)*m);
    }
    bzero(Q,sizeof(double)*m*n);
    for(int j=0;j<n;j++) Q[j][j]=1;
    for(int j=n-1;j>=0;j--){
        hheliminate(m,n,Hs[j],Q);
    }
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            if(i>j) R[i][j]=0;
        }
    }
}

void matrixMult(int m, int n, int q, double A[m][n], double B[n][q], double result[m][q]) {
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < q; j++) {
      result[i][j] = 0;
      for(int k = 0; k < n; k++) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

void transposeMatrix(int m, int n, double A[m][n], double result[n][m]) {
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      result[j][i] = A[i][j];
    }
  }
}

#define N 5

int main(){
    setrlimit(RLIMIT_STACK,
        &(const struct rlimit)
        {RLIM_INFINITY,RLIM_INFINITY});
    printf("N=%d\n",N);
    FILE *fp=fopen("file05a.dat","r");
    int m = 5;
    //fscanf(fp,"%d",&m);
    //printf("m=%d\n",m);
    double Q[m][N],R[m][N];
    double A[5][5]={
    {  3.1411235,  2.9670499, -2.0696759, -0.1696742, -0.4045575 },
    {  0.3868716,  3.6300832, -0.7903713, -0.1573674,  0.6810963 },
    {  0.1909349,  1.5450310,  1.8109176, -0.6418465,  0.3332143 },
    { -0.0049263, -1.7362550,  0.7578766,  0.0782331,  2.9166110 },
    { -0.0240325, -2.1275275,  0.3860839, -1.2922110,  6.3396425 }};
    /*double X[m],Y[m];
    for(int i=0;i<m;i++){
        fscanf(fp,"%lf %lf",&X[i],&Y[i]);
        printf("%lf %lf\n",X[i],Y[i]);
    }*/
    //householder(m,N,A,R,Q);
    
    /* for Gram-schmidt */
    double QT[N][m];
    gramschmidt(m, N,A,QT);
    matprint(m, N, A);
    transposeMatrix(N, m, QT, Q);
    /* end of gram-schmidt*/
    
    /*
    double qty[N];
    bzero(qty,sizeof(double)*N);
    
    
    
    for(int i=0;i<m;i++){
        for(int j=0;j<N;j++){
            qty[j]+=Q[i][j]*Y[i];
        }
    }*/
    
    /* FOR gram-schmidt */
    matrixMult(N, m, N, QT, A, R);
    
    
    
    //double c[N];
    //backsub(N,R,c,qty);
    //printf("c=\n"); vecprint(N,c);
    return 0;
}    
