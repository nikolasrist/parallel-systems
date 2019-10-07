#===============================================================================

# choose compiler

# GNU compiler
CC	= gcc
CFLAGS	= -Wall -std=c99

# Intel compiler
#CC	= icc
#CFLAGS	= -Wall -std=c99

# PGI compiler
#CC	= pgcc
#CFLAGS	=

# link flags
LDLIBS	= -lFHBRS

#===============================================================================
#
#       implicit rules
#

# Compiling
%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $<


#===============================================================================

DEFAULT:: prefix.exe

run:: prefix.exe
	./prefix.exe

clean::
	-rm -f *.o *.exe

#===============================================================================

prefix.exe: main.o prefix.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

main.o: main.c prefix.h

prefix.o: prefix.c prefix.h

#===============================================================================

