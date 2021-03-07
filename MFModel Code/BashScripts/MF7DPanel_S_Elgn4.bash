#!/bin/bash
#
##SBATCH --nodes=1
##SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=5:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=Run20Panel
##SBATCH --mail-type=END
##SBATCH --mail-user=zx555@nyu.edu
#SBATCH --output=%j_Run4_Panel%a.out

module purge
module load matlab/2020b

cd /scratch/$USER/NYU-Vision
echo ${SLURM_ARRAY_TASK_ID}
matlab -r "MF_7DContour_HPC(${SLURM_ARRAY_TASK_ID},4,5,4)"

