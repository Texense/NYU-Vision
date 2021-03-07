for SELGN in 2; do
for SILGN in $(seq 1 4); do
for rIL6  in $(seq 1 4); do
for SEE in 3; do
for SII in 2; do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o Fig2_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e Fig2_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=Fig2_${SILGN}_${rIL6} \
       Figure234Panel.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done