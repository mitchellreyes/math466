CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math -fcilkplus
all: householder

clean:
	rm *.o householder

householder: householder.c matrixlib.o
	$(CC) -o householder \
		householder.c matrixlib.o -lm -llapacke

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -c matrixlib.c

