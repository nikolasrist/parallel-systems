#!/bin/sh
#SBATCH --partition=wr43         # partition (queue)
#SBATCH --tasks=96               # number of tasks
#SBATCH --mem=16G                # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time=20:00             # total runtime of job allocation (format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID)
#SBATCH --error=slurm.%j.err     # filename for STDERR

# load latest GNU compiler and library with utility functions
module load gcc libFHBRS

# start application through Makefile
make run

# finished
exit
