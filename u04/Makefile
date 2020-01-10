# use either standard or optimized version

# standard version (needs any BLAS compatible library)
#CC	= cc
#CFLAGS	= -O -std=c99 -I/usr/local/include
#LINKER	= cc
#LDLIBS	= -L/usr/local/lib -lFHBRS -L /usr/lib64/atlas/ -lcblas -lblas -lm

# optimized version using Intel MKL (needs MKL installed and a module load intel-cpompiler done to set the MKL search path)
CC	= cc
CFLAGS	= -DMKL -O -std=c99 -I$(MKLROOT)/include -I/usr/local/include
LINKER	= cc
LDLIBS	= -L /usr/local/lib -lFHBRS -Wl,--no-as-needed -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

default:: matmul_c.exe

run:: matmul_c.exe
	./matmul_c.exe

clean::
	-rm -f *.o *.exe

matmul_c.exe: matmul_c.o
	$(LINKER) -o $@ $^ $(LDLIBS)

