#!/bin/bash
#SBATCH --account=ctb-rmansbac
#SBATCH --time=00:30:00
#SBATCH --job-name="esmfold test"
#SBATCH --nodes=1
#SBATCH --mem=32G

# Load modules/virtual environments
source venv/bin/activate

python pepstructure/esm_example.py GKIIKLKASLKLL