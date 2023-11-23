#!/bin/bash
# This file clones the RoseTTAFold repo, 
# then sets up a virtual environment for RoseTTAFold
# it assumes we do NOT want RoseTTAFold's pyRosetta structure refinement method

##########################################
##########################################
# get repo and set up virtual environment
##########################################
##########################################

cd .. # file assumes we begin in the project-benchmark-seq2struct/scripts directory
deactivate # deactivate any existing virtual environments

# clone RoseTTAFold
git clone https://github.com/RosettaCommons/RoseTTAFold

# move the RoseTTAFold pip requirements file into the RoseTTAFold directory
mv requirements_rosettafold.txt RoseTTAFold/

# move into the RoseTTAFold directory
cd RoseTTAFold 

# RoseTTAFold assumes python 3.8.10
module load python/3.8.10
python3.8 -m venv venv_RoseTTAFold

# activate the virtual environment
source venv_RoseTTAFold/bin/activate

# install the pip requirements and use compute canada mirrors (--no-index)
pip install -r requirements_rosettatfold.txt --no-index

##########################################
##########################################
# Download the RoseTTAFold model weights
##########################################
##########################################
wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
tar xfz weights.tar.gz

##########################################
##########################################
# Download three large databases: UniRef30_2020_06, BFD, and pdb100_2021Mar03
##########################################
##########################################

# uniref30 [46G]
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
mkdir -p UniRef30_2020_06
tar xfz UniRef30_2020_06_hhsuite.tar.gz -C ./UniRef30_2020_06

# BFD [272G]
wget https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz
mkdir -p bfd
tar xfz bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz -C ./bfd

# structure templates [10G]
wget https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2021Mar03.tar.gz
tar xfz pdb100_2021Mar03.tar.gz