#include <stdio.h>

int main(){
  double a = 1.0, b = 2.0, c = 3.0;
  printf("a = %.14e\nb = %.14e\nc = %.14e\n", a, b, c);
  double r1 = (a+b) + c;
  double r2 = a + (b + c);

  if(r1 == r2){
    printf("In this case a+b+c is associative.\n");
  }
  else{
    printf("In this case a+b+c is NOT associative\n");
  }

  return 0;
}
