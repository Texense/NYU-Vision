#!/bin/bash
##SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=64GB
##SBATCH --job-name=Run20Panel
##SBATCH --mail-type=END
##SBATCH --mail-user=zx555@nyu.edu
##SBATCH --output=%j_Run1_Panel%a.out

module purge
module load matlab/2020b

cd /scratch/$USER/NYU-Vision
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
matlab -r "Figure6_HPCSpotCheck(${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII})"








