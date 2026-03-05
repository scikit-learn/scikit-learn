#!/bin/bash

set -e
set -x

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "${HOME}/conda"
source "${HOME}/conda/etc/profile.d/conda.sh"


# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh
conda activate base

# Run these debug commands before installing our specific conda environment.
# We want to see what is available on the runner before we make changes. But
# we need to install miniforge before being able to look at the output of the
# conda commands.
conda info --json | python -c "import sys, json; print('Conda virtual packages versions:', json.load(sys.stdin).get('virtual_pkgs', []));"
nvidia-smi

CONDA_ENV_NAME=sklearn
LOCK_FILE=build_tools/github/pylatest_conda_forge_cuda_array-api_linux-64_conda.lock
create_conda_environment_from_lock_file $CONDA_ENV_NAME $LOCK_FILE

conda activate $CONDA_ENV_NAME
conda list
