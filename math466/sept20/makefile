CC=/usr/local/gcc-7.1/bin/gcc
multAx: multAx.c matrixlib.o tictoc.o
	gcc -Wall -std=gnu99 -o multAx multAx.c matrixlib.o tictoc.o -lm -lcilkrts

matrixlib.o: matrixlib.c matrixlib.h
	gcc -fcilkplus -std=gnu99 -c matrixlib.c

tictoc.o: tictoc.c tictoc.h
	gcc -std=gnu99 -c tictoc.c
