for SELGN in 2 3; do
for SILGN in 2 3; do
for rIL6  in 6; do
for SEE in $(seq 1 3); do
for SII in $(seq 1 5); do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o out_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e out_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=MF7D_${SELGN}_${SILGN}_${rIL6} \
       MF7DPanel_S_EE_SII_Final.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done