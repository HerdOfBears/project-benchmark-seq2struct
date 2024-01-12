#!/bin/bash
#SBATCH --job-name=fix_pdbs
#SBATCH --account=ctb-rmansbac
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=32G

# load modules/venv
source venv/bin/activate
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0
pip install --no-index pdbfixer

# move from scripts/ to project-benchmarking/
cd ../
pwd

python scripts/fix_pdb_files.py --pdb_dir outputs/esmfold/