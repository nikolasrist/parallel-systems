CFLAGS	= -std=c99 -Wall -O -g
LDLIBS	= -lmetis -lFHBRS -lm


# compile everything
default:: part.exe

run:: part.exe
	./part.exe

# clean
clean::
	-rm -f part.exe part.o

part.exe: part.o
	$(CC) -o $@ $^ $(LDLIBS)


