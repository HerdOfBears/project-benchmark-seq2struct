#!/bin/bash
cd .. # file assumes we are located in project-benchmark-seq2struct directory
module load python/3.9
python3.9 -m venv venv

echo "activating virtual environment"
source venv/bin/activate

# check if requirements.txt file in directory
if test -f "requirements.txt"; then
	echo "continuing with package installations"
	echo "commencing pip install of requirements.txt"
	pip install -r requirements.txt --no-index

	# OpenFold requires nvcc availability 
	echo "openfold requires loading cuda"
	module load cuda
	pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
fi
