#!/bin/sh
#SBATCH --partition hpc3         # partition (queue)
#SBATCH --ntasks=64              # number of tasks
#SBATCH --ntasks-per-core=1      # don't use HyperThreading
#SBATCH --mem=4G                 # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time 10:00             # total runtime of job allocation ((format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID)
#SBATCH --error=slurm.%j.err     # filename for STDERR

# load modules
module load gcc openmpi/gnu libFHBRS

# start here your MPI program (example here: ring.exe)
mpirun -np 64 ./ring.exe

exit
