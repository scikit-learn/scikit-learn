#!/bin/bash

set -e
set -x

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "${HOME}/conda"
source "${HOME}/conda/etc/profile.d/conda.sh"


# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh
conda activate base

CONDA_ENV_NAME=sklearn
LOCK_FILE=build_tools/github/pylatest_conda_forge_cuda_array-api_linux-64_conda.lock
create_conda_environment_from_lock_file $CONDA_ENV_NAME $LOCK_FILE

conda activate $CONDA_ENV_NAME
conda list
