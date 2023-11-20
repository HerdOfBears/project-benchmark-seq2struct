#!/bin/bash
cd .. # file assumes you have moved into scripts/ directory, but commands assume project-benchmark-seq2struct/ directory
module load python/3.9 # load python 3.9, required python version for esmfold
python3.9 -m venv venv # make a virtual environment under project-benchmark-seq2struct/

echo "activating virtual environment"
source venv/bin/activate # activate the virtual environment

# check if requirements.txt file in directory
if test -f "requirements.txt"; then
	sed -i '/openfold/d' requirements.txt # remove openfold from requirements.txt otherwise errors will occur
	echo "continuing with package installations"
	echo "commencing pip install of requirements.txt"
	pip install -r requirements.txt --no-index # install pacakges using Compute Canada's wheels

	# OpenFold requires nvcc availability 
	echo "openfold requires loading cuda"
	module load cuda # openfold requires nvcc installed
	pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
fi
