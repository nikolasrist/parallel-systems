CC	= gcc
CFLAGS	= -std=c99 -Wall -O2

default: bubblesort.exe

run:bubblesort.exe
	./bubblesort.exe 10000

clean:
	-rm -f *.o *.exe

%.exe: %.o
	$(CC) -o $@ $< $(LDLIBS)
