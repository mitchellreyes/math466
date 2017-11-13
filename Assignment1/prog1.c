#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main(){
  double nums[3] = {1.0, 1.0, 1.0};
  double r1, r2;
  int i, rows = 0;
  srand(time(NULL));
  for(i = 0; i < 3; i++){
    while(nums[i] < 100.0){
      if(nums[i] + 0.1 > 100.0){
        break;
      }
      nums[i] += 0.1;
      r1 = (nums[0] + nums[1]) + nums[2];
      r2 = nums[0] + (nums[1] + nums[2]);

      if(r1 != r2){
        printf("a = %.14lf\nb = %.14lf\nc = %.14lf\n", nums[0], nums[1], nums[2]);
        printf("In this case a+b+c is NOT associative\n\n");
      }
    }
  }

  return 0;
}
