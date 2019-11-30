# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7
MPICC = mpicc
MPIRUN = mpirun -np 4

test_sequential:
	tar -xvzf code.tar.gz
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/knnring.h ../; cd ..
	$(CC) tester.c knnring_sequential.a -o $@ -lm -lopenblas
	./test_sequential
 
