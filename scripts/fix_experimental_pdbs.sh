#!/bin/bash
#SBATCH --job-name=fix_expt_pdbs
#SBATCH --account=def-rmansbac
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=12G

# load modules/venv
source venv/bin/activate
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load openmm/8.0.0

# run script
DIR="../../shared/Samith/Benchmark_proj/PDB_exp"
for dir in $DIR/*; do
    if [ -d "$dir" ]; then
        echo "Processing $dir"
        python scripts/fix_pdb_files.py --pdb_dir $dir"/Best_Structure/" --experimental y
    fi
done