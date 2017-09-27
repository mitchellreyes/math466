#include <stdio.h>
#include <stdlib.h>
#include <readline/readline.h>
#include <readline/history.h>

long double epsilon(long double x){
    long double a = 0, b = x;
    for(;;){
        long double y = (a + b) / 2;
        long double r = x + y;
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
        long double x=strtod(xp,&ep);
        if(xp==ep||*ep!=0) {
            printf("Not a number, try again!\n");
            continue;
        }
        if(x==0.0) {
            printf("Zero detected--exiting!\n");
            exit(0);
        }
        printf("The value of %Lg was entered for x.\n",x);
        long double e = epsilon(x);
        printf("The value of epsilon is %Lg. \n", e);
        printf("The ratio epsilon/x = %Lg. \n", e/x);
    }
}
