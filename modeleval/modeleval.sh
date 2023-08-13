#!/bin/bash
#SBATCH -J modeleval
#SBATCH -o modeleval_%j.out
#SBATCH -e modeleval_%j.err
#SBATCH -N 1
#SBATCH -c 12
#SBATCH -t 0-03:00
#SBATCH -p huce_ice
#SBATCH --mem=10000
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2020b-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "runtests_23_06_04.m;"
