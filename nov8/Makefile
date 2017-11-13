CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math
all: household

household: household.c matrixlib.o tictoc.o
	$(CC) -fcilkplus -Wall -std=gnu99 -o household \
		household.c matrixlib.o tictoc.o -lm -lcilkrts -llapacke

specrad: specrad.c matrixlib.o tictoc.o
	$(CC) -fcilkplus -Wall -std=gnu99 -o specrad \
		specrad.c matrixlib.o tictoc.o -lm -lcilkrts -llapacke

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -fcilkplus -std=gnu99 -c matrixlib.c
	
tictoc.o: tictoc.c tictoc.h
	$(CC) -fcilkplus -std=gnu99 -c tictoc.c

clean:
	rm *.o specrad household
