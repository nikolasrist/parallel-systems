#CC	= icc -openmp
CC	= gcc -fopenmp
CFLAGS	= -O2 -g -std=c99 -Wall -I/usr/local/include
LDFLAGS	= -L/usr/local/lib -lFHBRS -lX11 -lm

default:: heat.exe

run:: heat.exe
	./heat.exe 2000 5000 0

clean::
	-rm -f heat.exe heat.o

heat.exe: heat.o
	$(CC) -o $@ $^ $(LDFLAGS)


