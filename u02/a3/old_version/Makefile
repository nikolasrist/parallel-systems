########################################################################
# do before module load gcc openmpi/gnu

MPICC		= mpicc
CFLAGS		= -O2 -std=c99 -Wall
LDLIBS		= -lFHBRS -lm

########################################################################

compile:: ring.exe

clean::
	-rm -f *.exe *.o

run:: ring.exe
# only as an example. possible only inside a batch script
	mpirun -np 8 ring.exe


########################################################################

%.o: %.c
	$(MPICC) $(CFLAGS) -c $<
%.exe: %.o
	$(MPICC) -o $@ $< $(LDLIBS)

########################################################################
