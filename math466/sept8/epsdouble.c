#include <stdio.h>
#include <stdlib.h>
#include <readline/readline.h>
#include <readline/history.h>

double epsilon(double x){
    double a = 0, b = x;
    for(;;){
        double y = (a + b) / 2;
        double r = x + y;
        if(r == x){
            if(a == y){
                return a;
            }
            a = y;    
        }
        else{
            if(b == y){
                return a;
            }
            b = y;
        }
    }
}

int main(){
    for(;;){
        char *ep,*xp=readline("Please enter x: ");
        if(!xp) {
            printf("\nEnd of file detected--exiting!\n");
            exit(0);
        }
        if(*xp==0) continue;
        add_history(xp);
        double x=strtod(xp,&ep);
        if(xp==ep||*ep!=0) {
            printf("Not a number, try again!\n");
            continue;
        }
        if(x==0.0) {
            printf("Zero detected--exiting!\n");
            exit(0);
        }
        printf("The value of %g was entered for x.\n",x);
        double e = epsilon(x);
        printf("The value of epsilon is %g. \n", e);
        printf("The ratio epsilon/x = %g. \n", e/x);
    }
}
