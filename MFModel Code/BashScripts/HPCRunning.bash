#!/bin/bash
#
##SBATCH --nodes=1
##SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=100GB
#SBATCH --job-name=VisionMF7D
#SBATCH --mail-type=END
#SBATCH --mail-user=zx555@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2020b

cd /scratch/$USER/NYU-Vision
##cat HPCTrail.m | srun matlab -nodisplay
cat MF_7DContour_ScriptServer.m | srun matlab -nodisplay


