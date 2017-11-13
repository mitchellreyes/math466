CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math
all: lufact

lufact: lufact.c matrixlib.o
	$(CC) -fcilkplus -Wall -std=gnu99 -o lufact \
		lufact.c matrixlib.o -lm -lcilkrts

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -fcilkplus -std=gnu99 -c matrixlib.c

