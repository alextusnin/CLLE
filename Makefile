CC=g++
CFLAGS_debug = -lfftw3 -lm -O0 
CFLAGS_run = -lfftw3 -lm -O3 


debug: main.c 
	$(CC) main.c $(CFLAGS_debug) -o main_LLE.o
run: main.c 
	$(CC) main.c $(CFLAGS_run) -o main_LLE.o
