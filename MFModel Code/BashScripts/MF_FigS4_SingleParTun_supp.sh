## ParRatIn divided by 100
for ParInd in 9 10 11 12 15; do
for ParRatIn in $(seq 40 79) $(seq 121 160); do
#
echo "${ParInd}, ${ParRatIn}"
export ParInd ParRatIn
#
sbatch -o FigS4_ParInd${ParInd}_ParRatIn${ParRatIn}.stdout.txt \
       -e FigS4_ParInd${ParInd}_ParRatIn${ParRatIn}.stdout.txt \
       --job-name=PT_${ParInd}_${ParRatIn} \
       MF_FigS4_Run.bash
#
sleep 1 # pause to be kind to the scheduler
done
done
