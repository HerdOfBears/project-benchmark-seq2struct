#!/bin/bash
#SBATCH --job-name=not_gnu_esmfold_predns_mdruns
#SBATCH --account=ctb-rmansbac
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=24G

# Assumes running from project-benchmark-seq2struct directory

# load modules/venv
source venv/bin/activate
module load StdEnv/2020 cuda/11.4 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0


# "outputs/" is deep learning model outputs, but MD sim inputs
start=`date +%s.%N` 
for spep in starPep_00218 starPep_43458 starPep_12300 starPep_09923; do
    python scripts/simulate_protein_in_water.py --input_dir outputs/esmfold/ --pdb_file ${spep}_esmfold_prediction_pdbfixed.pdb --output_dir outputs/gnu_test/ --slurm_id $SLURM_JOB_ID --prefix non_gnu_test
done
end=`date +%s.%N`

runtime=$( echo "$end - $start" | bc -l )
echo "duration: $runtime"
echo Done.