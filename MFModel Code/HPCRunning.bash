#!/bin/bash
#
##SBATCH --nodes=1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --time=1:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=myMatlabTest
#SBATCH --mail-type=END
##SBATCH --mail-user=zx555@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2020a

cd /scratch/$USER/NYU-Vision
cat HPCTrail.m | srun matlab -nodisplay


