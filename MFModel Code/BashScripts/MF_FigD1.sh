for SELGN in 2; do
for SILGN in 2; do
for rIL6  in 2; do
for SEE in 3; do
for SII in $(seq 2 3); do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o FigD1_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e FigD1_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=FigD1_${SEE}_${SII} \
       FigureD1Panel.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done