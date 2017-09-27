//compile using gcc sn0p.c -o sn0p -lm -std=c99 -fcilkplus

#include <stdio.h>
#include <math.h>
#include <cilk/cilk.h>

double Hser(int p,int q){
    double hn=0.0;
    for(int k=q;k>=p;k--) hn+=1.0/k;
    return hn;
}
double Hpar(int p,int q){
    if(q-p<20000000) return Hser(p,q);
    int c=p/2+q/2;
    //puts the work on another thread
    double r1=cilk_spawn Hpar(p,c);
    //continues to work on this current processor
    double r2=Hpar(c+1,q);
    //waits for the work to be finished
    cilk_sync;
    return r1+r2;
}

int main() {
    int n=2000000000;
    printf("gamma(%d) = %.14f\n",n,Hpar(1,n)-log(n));
    return 0;
}
