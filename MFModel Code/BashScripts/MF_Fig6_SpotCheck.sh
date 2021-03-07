## For Fig2
for SELGN in 2; do
for SILGN in $(seq 1 4); do
for rIL6  in $(seq 1 4); do
for SEE in 3; do
for SII in 2; do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o Fig2SC_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e Fig2SC_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=Fig2_${SILGN}_${rIL6} \
       MF_Fig6_SpotCheck.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done

## For Fig3
for SELGN in $(seq 1 4); do
for SILGN in $(seq 1 4); do
for rIL6  in 2 3; do
for SEE in 3; do
for SII in 2; do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o Fig3SC_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e Fig3SC_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=Fig3_${SELGN}_${SILGN}_${rIL6} \
       Figure6SpotCheck.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done

## For Fig 4
for SELGN in 3; do
for SILGN in 2; do
for rIL6  in 3; do
for SEE in $(seq 1 5); do
for SII in $(seq 1 4); do
#
echo "${SELGN}, ${SILGN}, ${rIL6}, ${SEE}, ${SII}"
export SELGN SILGN rIL6 SEE SII
#
sbatch -o Fig4SC_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       -e Fig4SC_SElgn${SELGN}_SIlgn${SILGN}_rIL6${rIL6}_SEE${SEE}_SII${SII}.stdout.txt \
       --job-name=Fig4_${SEE}_${SII} \
       Figure6SpotCheck.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done