#!/bin/bash

wget https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ./miniconda.sh
bash ./miniconda.sh -b -p $HOME/miniconda

export PATH="$HOME/miniconda/bin:$PATH"

echo -e '\n# added by SNAPPy to run miniconda' >> ~/.bashrc

conda init
conda-env create -f environment.yaml
#conda activate test_subtyper

echo "Please open a new terminal and activate de environment: 'conda activate snappy'"
