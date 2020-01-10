#!/bin/bash
#SBATCH --output slurm.%N.%j.out # STDOUT
#SBATCH --error slurm.%N.%j.err  # STDERR
#SBATCH --partition wr43         # partition (queue)
#SBATCH --ntasks=96              # per node parallelism degree
#SBATCH --mem 10G                # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time 20:00             # total runtime of job allocation (format D-HH:MM:SS)

# for MKL library
module load intel-compiler

# execute the program
./matmul_c.exe

exit
