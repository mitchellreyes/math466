#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define size 3

void multAx(int m,int n, double A[m][n],double X[n],double B[m]);
double vecnorm2(int n, double x[n]);
double matnorm2(int n, double A[n][n]);
double determinant(double A[size][size], float k);
void cofactor(double num[size][size],float f, double qk);
void transpose(double num[size][size],
  double fac[size][size],double r, double qk);
void matcond2(double A[size][size]);

int main(){
  double A[size][size] = {{7, 1, 7}, {5, 7, 7}, {3, 7, 8}};
  matcond2(A);
  return 0;
}


double determinant(double A[size][size], float k){
  double s=1,det=0,b[size][size];
  int i,j,m,n,c;
  if (k==1)
    {
     return (A[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        n=0;
        for (i=0;i<k;i++)
          {
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=A[i][j];
                   if (n<(k-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                     }
                   }
               }
             }
          det=det + s * (A[0][c] * determinant(b,k-1));
          s=-1 * s;
          }
    }

    return (det);

}
void cofactor(double num[size][size],float f, double qk){
 double b[size][size],fac[size][size];
 int p,q,m,n,i,j;
 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     n=0;
     for (i=0;i<f;i++)
     {
       for (j=0;j<f;j++)
        {
          if (i != q && j != p)
          {
            b[m][n]=num[i][j];
            if (n<(f-2))
             n++;
            else
             {
               n=0;
               m++;
               }
            }
        }
      }
      fac[q][p]=pow(-1,q + p) * determinant(b,f-1);
    }
  }
  transpose(num,fac,f, qk);
}
void transpose(double num[size][size],
  double fac[size][size],double r, double qk){
  int i,j;
  double b[size][size],inverse[size][size],d;

  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         b[i][j]=fac[j][i];
        }
    }
  d=determinant(num,r);
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        inverse[i][j]=b[i][j] / d;
        }
    }
    printf("||A-1||2 = %g\n", matnorm2(size, inverse));
   printf("A = %g\n", qk + matnorm2(size, inverse));
}

void matcond2(double A[size][size]){
  double qk = matnorm2(size, A);
  cofactor(A, size, qk);
}


void multAx(int m,int n,
    double A[m][n],double X[n],double B[m]){
    bzero(B,sizeof(double)*m);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            B[i]+=A[i][j]*X[j];
        }
    }
}

double vecnorm2(int n, double x[n]){
	double r = 0;
	for(int i = 0; i < n; i++){
		double t = fabs(x[i]);
		r += t*t;
	}
	return sqrt(r);
}

double matnorm2(int n, double A[n][n]){
  double B[n][n], y[n], yk[n];
  bzero(B, sizeof(double)*n*n);
  for(int k = 0; k < n; k++){
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        B[i][j] += A[k][i]*A[k][j];
      }
    }
  }
  for(int i = 0; i < n; i++){
    y[i] = 2.0*random()/RAND_MAX+1.0;
  }
  double q = 0, qk;
  for(int k = 1; k < 100*n; k++){
    multAx(n, n, B, y, yk);
    qk = vecnorm2(n, yk);
    for(int j = 0; j < n; j++){
      y[j] = yk[j] / qk;
    }
    if(fabs(qk - q) < 5e-15*qk){
      return sqrt(qk);
    }
    q = qk;
  }
  printf("matnorm2: Failed to converge!\n");
  return sqrt(qk);

}
