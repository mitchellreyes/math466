CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math
all: linpack

linpack: linpack.c matrixlib.o tictoc.o
	$(CC) -fcilkplus -Wall -std=gnu99 -o linpack \
		linpack.c matrixlib.o tictoc.o -lm -lcilkrts

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -fcilkplus -std=gnu99 -c matrixlib.c
	
tictoc.o: tictoc.c tictoc.h
	$(CC) -fcilkplus -std=gnu99 -c tictoc.c

