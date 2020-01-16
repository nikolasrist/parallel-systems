#!/bin/bash
#SBATCH --output slurm.%N.%j.out # STDOUT
#SBATCH --error slurm.%N.%j.err  # STDERR
#SBATCH --partition hpc3         # partition (queue)
#SBATCH --nodes=1                # number of tasks/cores
#SBATCH --ntasks=48              # per node parallelism degree
#SBATCH --mem 16G                # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time 30:00             # total runtime of job allocation (format D-HH:MM:SS)

# load modules necessary
module load gcc openmpi/gnu metis/5.1.0-32

# degree of parallelism to be tested
nprocs="2 3 5 9 17 33"

# loop over all combinations of processor numbers
NNODES=3

# start here your program, for example an MPI program
echo "running on " $NNODES " processors"
mpirun -np $NNODES mandelbrot.exe

exit
