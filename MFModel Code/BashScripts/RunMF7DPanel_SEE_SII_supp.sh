declare -a SElgnMat SIlgnMat rIl6Mat SEEMat SIIMat
SElgnMat=(4 4)
SIlgnMat=(2 2)
rIl6Mat=(3 6)
SEEMat=(2 1)
SIIMat=(2 4)
#
for MatInd in $(seq 0 1); do
#
SELGN=${SElgnMat[${MatInd}]}
SILGN=${SIlgnMat[${MatInd}]}
rIl6=${rIl6Mat[${MatInd}]}
SEE=${SEEMat[${MatInd}]}
SII=${SIIMat[${MatInd}]}
#
echo "${SELGN}, ${SILGN}, ${rIl6}, ${SEE}, ${SII}"
export SELGN SILGN rIl6 SEE SII
#
sbatch -o out_SElgn${SELGN}_SIlgn${SILGN}_rIl6${rIl6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e out_SElgn${SELGN}_SIlgn${SILGN}_rIl6${rIl6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=MF7D_${SELGN}_${SILGN}_${rIL6} \
       MF7DPanel_SEE_SII_Final.bash
#
sleep 1 # pause to be kind to the scheduler
done