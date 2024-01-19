#!/bin/bash
#SBATCH --job-name=check_device
#SBATCH --account=ctb-rmansbac
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --mem=32G

# move from scripts/ to project-benchmarking/
# cd ../
# pwd

# load modules/venv
source venv/bin/activate
module load StdEnv/2020 cuda/11.4 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0


python scripts/simulate_protein_in_water.py --pdb_dir outputs/esmfold/ --pdb_file starPep_03523_esmfold_prediction_pdbfixed.pdb