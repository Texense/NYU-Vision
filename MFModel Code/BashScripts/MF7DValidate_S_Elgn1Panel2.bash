#!/bin/bash
#
##SBATCH --nodes=1
##SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=Valid20Panel
##SBATCH --mail-type=END
##SBATCH --mail-user=zx555@nyu.edu
#SBATCH --output=%j_Validate1_Panel2.out

module purge
module load matlab/2020b

cd /scratch/$USER/NYU-Vision
echo ${SLURM_ARRAY_TASK_ID}
matlab -nodisplay -r "MF_7D_NWValidation(1,2)"


