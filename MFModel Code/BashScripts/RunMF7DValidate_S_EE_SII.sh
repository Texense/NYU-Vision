declare -A arr1=([1]=1 [2]=1 [3]=2 [4]=2 [5]=3 [6]=3)
declare -A arr2=([1]=1 [2]=5 [3]=3 [4]=4 [5]=1 [6]=4) 

for SELGN in 2 3; do
for SILGN in 2 3; do
for MatInd in $(seq 1 6); do
SEE=${arr1[${MatInd}]}
SII=${arr2[${MatInd}]}
#
echo "${SELGN}, ${SILGN}, ${SEE}, ${SII}"
export SELGN SILGN SEE SII
#
sbatch -o Validate_SElgn${SELGN}_SIlgn${SILGN}_SEE${SEE}_SII${SII}.stdout.txt \
-e Validate_SElgn${SELGN}_SIlgn${SILGN}_SEE${SEE}_SII${SII}.stdout.txt \
--job-name=MF7D_selgn${SELGN} \
MF7DValidate_SEE_SII.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
done