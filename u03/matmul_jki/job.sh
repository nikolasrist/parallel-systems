#!/bin/sh
##############################################################
#SBATCH --output slurm.%N.%j.out # STDOUT (%N: nodename, %j: job-ID)
#SBATCH --error slurm.%N.%j.err  # STDERR
#SBATCH --partition hpc          # partition (queue)
#SBATCH --ntasks 4               # number of tasks/cores
#SBATCH --mem 2G                 # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time 72:00:00          # total runtime of job allocation (format D-HH:MM:SS)
##############################################################
#
# generate a data memory reference trace for a program,
# run the dinero cache simulator with this memory reference trace
# trace data is send through a socket
#
##############################################################
# log everything
set -x

module load gcc pin dinero4

# which tool to use
TOOL=/usr/local/pin/pin-3.7/source/tools/memtrace/obj-intel64/memtrace.so

##############################################################

# log arguments submitted to script
echo "Running with PORT=$PORT,PROGRAM=$PROGRAM,PROGRAM_ARGS=$PROGRAM_ARGS,DIR=$DIR,FILE=$FILE,DINERO_ARGS=$DINERO_ARGS" >&2

# create output directory for this job (if it doesn't exist) and change to it
mkdir -p output/$DIR
cd output/$DIR

#!!!!!!! CHANGE port !!!!!!!!
# free socket port to use for communication between server and client
export PIN_ON_SOCKET=$PORT
#e.g. 14423

#!!!!!!! CHANGE program name !!!!!!!!
# program to instrument
#export PROGRAM="./bubblesort.exe"
#!!!!!!! CHANGE program arguments !!!!!!!!
#export PROGRAM_ARGS="30000"

# increase stack size for some programs
ulimit -s 10000000

##############################################################
# run server and client

# note: valgrind has also a (less accurate) memory tracer: "valgrind --tool=lackey --trace-mem=yes" with an appropriate filter
# start trace generator with application to instrument (including program arguments)
pin -t $TOOL -- $SLURM_SUBMIT_DIR/$PROGRAM $PROGRAM_ARGS &

# otherwise the server isn't ready to accept client requests on port
sleep 10

#!!!!!!! CHANGE dinero parameters !!!!!!!!
# memtrace client with cache simulation (change dinero parameters to your needs)
time memtrace_client.exe localhost $PIN_ON_SOCKET | \
	 dineroIV $DINERO_ARGS > $FILE

##############################################################

exit
