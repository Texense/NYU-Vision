#!/bin/bash
##SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=20:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=Figure1HRes
#SBATCH --mail-type=END
#SBATCH --mail-user=zx555@nyu.edu
#SBATCH --output=Figure1.out

module purge
module load matlab/2020b

cd /scratch/$USER/NYU-Vision
matlab -r "Figure1_SuperFineContour"








