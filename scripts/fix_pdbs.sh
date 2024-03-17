#!/bin/bash
#SBATCH --job-name=fix_pdbs
#SBATCH --account=ctb-rmansbac
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=12G

# load modules/venv
source venv/bin/activate
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0

# run script
python scripts/fix_pdb_files.py --pdb_dir outputs/af2/
