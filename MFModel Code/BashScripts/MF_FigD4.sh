##We are doing different axes... 
declare -A arr1=([1]=1 [2]=3 [3]=5 [4]=1 [5]=3 [6]=6)
declare -A arr2=([1]=2 [2]=4 [3]=6 [4]=5 [5]=6 [6]=7) 
##declare -A arr1=([1]=1 [2]=1)
##declare -A arr2=([1]=2 [2]=5) 

for SeqInd in $(seq 1 6); do 
ParInd1=${arr1[${SeqInd}]}
ParInd2=${arr2[${SeqInd}]}
#
echo "${ParInd1}, ${ParInd2}"
export ParInd1 ParInd2
#
sbatch -o FigD4_${ParInd1}_${ParInd2}.stdout.txt \
       -e FigD4_${ParInd1}_${ParInd2}.stdout.txt \
       --job-name=D4_${ParInd1}_${ParInd2} \
       FigureD4Panel.bash
#
sleep 1 # pause to be kind to the scheduler
done
