#!/bin/bash
#SBATCH --job-name=fix_pbc_image
#SBATCH --account=ctb-rmansbac
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

# Assumes running from project-benchmark-seq2struct directory

# load modules/venv
source venv/bin/activate
#module load StdEnv/2020 cuda/11.4 gcc/9.3.0 openmpi/4.0.3
#module load openmm/8.0.0

# "outputs/" is deep learning model outputs, but MD sim inputs
python scripts/fix_pbc_image.py --input outputs/esmfold/starPep_00385/starPep_00385_esmfold_prediction_coordinates_peptideOnly.pdb --output outputs/esmfold/starPep_00385/starPep_00385_esmfold_predn_coords_peptideOnly_pbcfixed.pdb
echo Done.
