all: lufact

lufact: lufact.c matrixlib.o tictoc.o
		-Wall -std=gnu99 -o lufact lufact.c matrixlib.o tictoc.o -lm -lcilkrts

matrixlib.o: matrixlib.c matrixlib.h
	-std=gnu99 -c matrixlib.c

tictoc.o: tictoc.c tictoc.h
	-std=gnu99 -c tictoc.c

clean:
	rm *.o lufact
