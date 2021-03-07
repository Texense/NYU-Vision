for SELGN in 2 3; do
for SILGN in 2 3; do
for SEE in $(seq 1 3); do
for SII in $(seq 1 5); do
#
echo "${SELGN}, ${SILGN}, ${SEE}, ${SII}"
export SELGN SILGN SEE SII
#
sbatch -o out_SElgn${SELGN}_SIlgn${SILGN}_SEE${SEE}_SII${SII}.stdout.txt \
-e out_SElgn${SELGN}_SIlgn${SILGN}_SEE${SEE}_SII${SII}.stdout.txt \
--job-name=MF7D_selgn${SELGN} \
MF7DPanel_S_EE_SII.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done