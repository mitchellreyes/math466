#include <stdio.h>
#include <math.h>

int main() {
    int n=2000000000;
    double hn=0.0;
    for(int k=n;k>=1;k--) hn+=1.0/k;
    printf("gamma(%d) = %.14f\n",n,hn-log(n));
    return 0;
}
