for SEI in 14 20; do #$(seq 1 21)
#
echo "${SEI}"
export SEI
#
sbatch -o FigD1b_SEI${SEI}.stdout.txt \
       -e FigD1b_SEI${SEI}.stdout.txt \
       --job-name=D1b_${SEI} \
       FigureD1bPanel.bash
#
sleep 1 # pause to be kind to the scheduler
done
