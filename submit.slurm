#!/bin/bash --login

#SBATCH --job-name=Hayward_Rupture
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

#SBATCH --account=n03-tnm
#SBATCH --partition=highmem
#SBATCH --qos=highmem

cd Hayward_Rupture

export OMP_NUM_THREADS=1

module load PrgEnv-gnu
module load cray-fftw/3.3.10.3
module load metis/5.1.0
module load cray-hdf5/1.12.0.3
module load cray-parallel-netcdf/1.12.3.1
module load cmake/3.21.3

srun --distribution=block:block --hint=nomultithread ./axisem3d > axisem.out
