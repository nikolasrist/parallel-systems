#!/bin/bash
#SBATCH --output slurm.%N.%j.out # STDOUT
#SBATCH --error slurm.%N.%j.err  # STDERR
#SBATCH --partition wr43         # partition (queue)
#SBATCH --ntasks=96              # per node parallelism degree
#SBATCH --mem 16G                # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time 30:00             # total runtime of job allocation (format D-HH:MM:SS)

module load intel-compiler

./sched_test.exe

exit
