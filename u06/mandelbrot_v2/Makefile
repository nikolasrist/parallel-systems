MPICC	= mpicc
CFLAGS	= -O2 -std=c99 -Wall -I/usr/local/include
LDFLAGS	= -lmetis -L/usr/local/lib -lFHBRS -L/usr/X11R6/lib -lX11

default:: mandelbrot.exe

run:: mandelbrot.exe
	mpirun -np 3 ./mandelbrot.exe

clean::
	-rm -f *.exe *.o


mandelbrot.exe: mandelbrot.o
	$(MPICC) -o $@ $< $(LDFLAGS)

mandelbrot.o: mandelbrot.c
	$(MPICC) $(CFLAGS) -c $<

