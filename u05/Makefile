CC	= icc -qopenmp
CFLAGS	= -std=c99 -O2 -g -I/usr/local/include
LDLIBS	= -L/usr/local/lib -lFHBRS -lpthread -lm

default:: sched_test.exe

run::
	./sched_test.exe

clean::
	-rm -f *.exe *.o


sched_test.exe: sched_test.o sched.o
	$(CC) -o $@ $^ $(LDLIBS)


sched.o: sched.c sched.h
sched_test.o: sched_test.c sched.h
