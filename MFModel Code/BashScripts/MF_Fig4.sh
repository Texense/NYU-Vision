for SELGN in 3; do
for SILGN in 2; do
for rIL6  in 3; do
for SEE in $(seq 1 5); do
for SII in $(seq 1 4); do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o Fig4_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e Fig4_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=Fig4_${SEE}_${SII} \
       Figure234Panel.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done