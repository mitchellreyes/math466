#include <stdio.h>
#include <stdlib.h>
#include <readline/readline.h>
#include <readline/history.h>

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
    }
}
