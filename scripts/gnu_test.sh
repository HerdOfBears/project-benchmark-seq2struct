#!/bin/bash
#SBATCH --job-name=esmfold_predns_mdruns
#SBATCH --account=ctb-rmansbac
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=2
#SBATCH --mem=24G

# Assumes running from project-benchmark-seq2struct directory

# load modules/venv
source venv/bin/activate
module load StdEnv/2020 cuda/11.4 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0

pwd
WDIR=$(pwd)
# "outputs/" is deep learning model outputs, but MD sim inputs
start=`date +%s.%N` 
parallel -j 2 --joblog parallel_test.log < ./scripts/gnu_test_sim.txt
end=`date +%s.%N`

runtime=$( echo "$end - $start" | bc -l )
echo "duration: $runtime"
echo Done.
echo Done.
