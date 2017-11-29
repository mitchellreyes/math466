#include <stdio.h>
#include <strings.h>
#include <complex.h>
#include <cilk/cilk.h>
#include <math.h>
#include "matrixlib.h"

void matprint(int m,int n,double A[m][n]){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("%g ",A[i][j]);
        }
        printf("\n");
    }
}
void cmatprint(int m,int n,complex A[m][n]){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("(%g %g) ",A[i][j]);
        }
        printf("\n");
    }
}
void vecprint(int m,double X[m]){
    matprint(m,1,(double (*)[1])X);
}
void cvecprint(int m,complex X[m]){
    cmatprint(m,1,(complex (*)[1])X);
}
// Algorithm (a) from page 13 of Numerical Algorithms by Justin Solomon
void multAx(int m,int n,double A[m][n],
    double X[n],double B[m]){
    bzero(B,m*sizeof(double));
    cilk_for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]+=A[i][j]*X[j];
        }
    }
}
void multATx(int m,int n,double A[m][n],
    double X[m],double B[n]){
    bzero(B,n*sizeof(double));
    for(int j=0;j<m;j++){
        cilk_for(int i=0;i<n;i++){
            B[i]+=A[j][i]*X[j];
        }
    }
}
void cmultAx(int m,int n,complex A[m][n],
    complex X[n],complex B[m]){
    bzero(B,m*sizeof(complex));
    cilk_for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]=B[i]+A[i][j]*X[j];
        }
    }
}
void LUfact(int n,double A[n][n]){
    for(int j=0;j<n-1;j++){
        for(int i=j+1;i<n;i++){
            double alpha=A[i][j]/A[j][j];
            for(int k=j+1;k<n;k++){
                A[i][k]-=alpha*A[j][k];
            }
            A[i][j]=alpha;
        }
    }
}
void backsub(int n,double U[n][n],
    double x[n],double b[n]){
    for(int i=n-1;i>=0;i--){ // Solve Ux=b
        x[i]=b[i];
        for(int j=i+1;j<n;j++){
            x[i]-=U[i][j]*x[j];
        }
        x[i]/=U[i][i];
    }
}
void LUsolve(int n,double LU[n][n],
    double x[n],double b[n]){
    double y[n];
    for(int i=0;i<n;i++){ // Solve Ly=b
        y[i]=b[i];
        for(int j=0;j<i;j++){
            y[i]-=LU[i][j]*y[j];
        }
    }
    backsub(n,LU,x,y);
}
void PLUsolve(int n,double LU[n][n],double *P[n],
    double x[n],double b[n]){
    double y[n];
    for(int i=0;i<n;i++){ // Solve Ly=b
        y[i]=b[(P[i]-&LU[0][0])/n];
        for(int j=0;j<i;j++){
            y[i]-=P[i][j]*y[j];
        }
    }
    for(int i=n-1;i>=0;i--){ // Solve Ux=y
        x[i]=y[i];
        for(int j=i+1;j<n;j++){
            x[i]-=P[i][j]*x[j];
        }
        x[i]/=P[i][i];
    }
}

void PLUfact(int n,double A[n][n],double *P[n]){
    for(int i=0;i<n;i++){
        P[i]=&A[i][0];
    }
    for(int j=0;j<n-1;j++){
        for(int i=j+1;i<n;i++){
            if(fabs(P[j][j])<fabs(P[i][j])){
                double *t=P[j]; P[j]=P[i]; P[i]=t;
//                printf("Swapped row %d with row %d\n",j,i);
            }
        }
        cilk_for(int i=j+1;i<n;i++){
            double alpha=P[i][j]/P[j][j];
            for(int k=j+1;k<n;k++){
                P[i][k]-=alpha*P[j][k];
            }
            P[i][j]=alpha;
        }
    }
}

double vecnorm1(int n,double x[n]){
    double r=0;
    for(int i=0;i<n;i++){
        r+=fabs(x[i]);
    }
    return r;
}
double matnorm1(int n,double A[n][n]){
    double r=0;
    for(int j=0;j<n;j++){
        double s=0;
        for(int i=0;i<n;i++){
            s+=fabs(A[i][j]);
        }
        if(s>r) r=s;
    }
    return r;
}
double vecnorm2(int n,double x[n]){
    double r=0;
    for(int i=0;i<n;i++){
        r+=x[i]*x[i];
    }
    return sqrt(r);
}
    
double matnorm2(int n,double A[n][n]){
/*    double B[n][n];
    bzero(B,sizeof(double)*n*n);
    for(int k=0;k<n;k++){
        for(int i=0;i<n;i++){   // B = A^t * A
            for(int j=0;j<n;j++){
                B[i][j]+=A[k][i]*A[k][j];
            }
        }
    }
    printf("here\n");
*/
    double x[n],y[n],yk[n],z[n];
    for(int i=0;i<n;i++){
        x[i]=2.0*random()/RAND_MAX+1.0;
    }
//    multAx(n,n,B,x,y);
    multAx(n,n,A,x,z);
     multATx(n,n,A,z,y);
    double normy,normyk,qold=0,qk;
    for(int k=2;k<20000;k++){
//        multAx(n,n,B,y,yk);
        multAx(n,n,A,y,z);
        multATx(n,n,A,z,yk);
        normy=vecnorm2(n,y);
        normyk=vecnorm2(n,yk);
 //       printf("%g/%g = %.15e\n",
 //           normyk,normy,normyk/normy);
        for(int j=0;j<n;j++){
            y[j]=yk[j]/normyk;
        }
        qk=normyk/normy;
        if(fabs(qk-qold)<=qk*5e-14){
            printf("Converged at k=%d\n",k);
            return sqrt(qk);
        }
        qold=qk;
//        memcpy(y,yk,sizeof(double)*n);
    }
    printf("Didn't converge!\n");
    return sqrt(normyk/normy);
}

