CC=/usr/local/gcc-7.1/bin/gcc
all: fft multAx

fft: fft.c matrixlib.o tictoc.o
	$(CC) -Wall -std=gnu99 -o fft \
		fft.c matrixlib.o tictoc.o -lm -lcilkrts

multAx: multAx.c matrixlib.o tictoc.o
	$(CC) -Wall -std=gnu99 -o multAx \
		multAx.c matrixlib.o tictoc.o -lm -lcilkrts

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -fcilkplus -std=gnu99 -c matrixlib.c

tictoc.o: tictoc.c tictoc.h
	$(CC) -std=gnu99 -c tictoc.c
