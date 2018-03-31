#!/bin/bash
#SBATCH --time=0:1:00
#SBATCH --job-name=fractal
#SBATCH --nodes=3
#SBATCH --cpus-per-task=16
#SBATCH --partition=YOUR-PARTITION

module load mpi/openmpi-x86_64
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

time mpirun -n 3 ./fract.exe {insert input file}
