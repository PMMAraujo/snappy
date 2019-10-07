#!/bin/bash

# warning
read -p "This scrip in going to install the conda distribution miniconda. If you already have conda installed on your machine this is not recomended.`echo $'\n> '`  Press [Enter] if you want to proceed or CTRL+C if you want to cancel"  response


# download miniconda3 version 4.6.14 from repository
wget https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ./miniconda.sh
# run miniconda instalation script
bash ./miniconda.sh -b -p $HOME/miniconda
rm miniconda.sh

# add conda path to available in this session
source ~/miniconda/etc/profile.d/conda.sh

# initiate conda and add the path to bash
conda init

# a warning to not close the terminal
echo "SNAPPy IS STILL INSTALLING PLEASE DO NOT CLOSE THE TERMINAL!"

# create snappy environment
conda env create -f environment.yaml

# activate snappy environment
conda activate snappy

# run tests to ensure instalation is correct
py.test

# instruction in how to future use snappy
echo ""
echo "Please open a new terminal and activate de environment: 'conda activate snappy'"
