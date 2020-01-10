#!/usr/bin/python

import sys

port=0
program=""

def usage():
    print "Usage: mkjobs.py port program arggroup1 [arggoup2 ...]"
    print "port:      smallest port to use"
    print "arggroupX: program arguments to use"
    print "the script will generate all job.sh calls for each arggroup"
    print ""
    print "Example2: mkjobs.py 12000 bubblesort.exe 1000 > runjobs.sh"
    print "Example1: mkjobs.py 13000 bubblesort.exe 100 200 > runjobs.sh"
    print ""
    print "This program is originally written by me, Jochen Wierum."
    print "You can do whatever you want with this code, but I won't take"
    print "responsibility for anything."

def mkDineroArgs(
        l1dsize="32k",  l1dassoc='8',
        l2usize="256k", l2uassoc='8',
        l3usize="32m",   l3uassoc='128',
        lsize='64'):
    return ['-l1-dsize '  + l1dsize,  '-l2-usize '  + l2usize, \
            '-l3-usize '  + l3usize,  '-l1-dbsize ' + lsize, \
            '-l2-ubsize ' + lsize,    '-l3-ubsize ' + lsize, \
            '-l1-dassoc ' + l1dassoc, '-l2-uassoc ' + l2uassoc, \
            '-l3-uassoc ' + l3uassoc]

def printCallString(pargs, name, value):
    global port
    args = {name: value}
    print 'sbatch --export=PORT=' + str(port) \
        + ',PROGRAM=' + program \
        + ',PROGRAM_ARGS=' + pargs \
        + ',DIR=' + pargs + '/' + name \
        + ',FILE=' + value + '.txt,DINERO_ARGS="' + ' '.join(mkDineroArgs(**args)) \
        + '",ALL job.sh'
    port = port + 1

def printCallStrings(args):
    sizes = ['1k', '2k', '4k', '8k', '16k', '64k', '128k', '256k', '512k']
    print "# L1 cache size: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'l1dsize', size)

    sizes = ['64k', '128k', '256k', '512k', '1m', '2m', '4m', '8m', '16m']
    print "# L2 cache size: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'l2usize', size)

    sizes = ['512k', '1m', '2m', '4m', '8m', '16m', '32m', '64m', '128m']
    print "# L3 cache size: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'l3usize', size)

    sizes = ['1', '2', '3', '4', '5', '6', '7', '8']
    print "# L1 associativity: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'l1dassoc', size)

    sizes = ['1', '2', '3', '4', '8', '16', '32', '64', '128']
    print "# L2 associativity: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'l2uassoc', size)

    sizes = ['1', '2', '3', '4', '8', '16', '32', '64', '128']
    print "# L3 associativity: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'l3uassoc', size)

    sizes = ['8', '16', '32', '64', '128', '256', '512', '1024']
    print "# line size: " + ', '.join(sizes)
    for size in sizes:
        printCallString(args, 'lsize', size)

print len(sys.argv)
if len(sys.argv) < 4:
    usage()
    sys.exit(1)
else:
    port = int(sys.argv[1])
    program = sys.argv[2]
    for arg in sys.argv[3:]:
        printCallStrings(arg)
