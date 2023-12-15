#!/bin/bash
#SBATCH --job-name=ESMFold_batch
#SBATCH --account=ctb-rmansbac
#SBATCH --time=05:30:00
#SBATCH --nodes=1

# load modules/venv
source venv/bin/activate

python scripts/1_perform_esm_predictions.py