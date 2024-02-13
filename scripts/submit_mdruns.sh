#!/bin/bash
#SBATCH --job-name=small_test_mdruns
#SBATCH --account=ctb-rmansbac
#SBATCH --time=0-01:00:00
#SBATCH --nodes=2
#SBATCH --gpus-per-node=4        # 6 nodes, 4 GPUs each, 24 total
#SBATCH --mem=15G                # mem per node
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=jyler.menard@mail.concordia.ca

# Assumes running from project-benchmark-seq2struct directory

pwd
WDIR=$(pwd)
JOBS_PER_NODE=4 # 4 GPUs per node

# load modules/venv
source $WDIR"/venv/bin/activate"
module load StdEnv/2020 cuda/11.4 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0

scontrol show hostname > ./node_list_${SLURM_JOB_ID} # save node list for parallel

# "outputs/" is deep learning model outputs, but MD sim inputs
start=`date +%s.%N` 
parallel -j $JOBS_PER_NODE \
        --joblog parallel_test.log \
        --sshloginfile $WDIR/node_list_${SLURM_JOB_ID} \
        --sshdelay 30 \
        --workdir $WDIR \
        < ./scripts/gnu_cmds.txt
end=`date +%s.%N`

runtime=$( echo "$end - $start" | bc -l )
echo "duration: $runtime"
echo Done.
