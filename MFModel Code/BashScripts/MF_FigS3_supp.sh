for SELGN in 1; do
for SILGN in 1; do
for rIL6  in 2; do
for SEE in 4; do
for SII in 2; do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o FigS3_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e FigS3_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=FS_${SEE}${SII}${SELGN}${SILGN}${rIL6} \
       Figure234Panel.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done