#!/bin/bash
#SBATCH --job-name=run_pipeline
#SBATCH --account=ctb-rmansbac
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=12G

# load modules/venv
module load pymol/2.5
module load dssp/3.1.4
source venv/bin/activate

python pepstructure/pipeline.py --model1 esmfold --model2 omegafold
