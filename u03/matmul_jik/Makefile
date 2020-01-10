CC	= gcc
CFLAGS	= -std=c99 -Wall -O2

default: matmul.exe

run: matmul.exe
	./matmul.exe 1000

clean:
	-rm -f *.o *.exe

%.exe: %.o
	$(CC) -o $@ $< $(LDLIBS)
