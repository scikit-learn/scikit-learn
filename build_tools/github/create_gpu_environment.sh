#!/bin/bash

set -e
set -x

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "${HOME}/conda"
source "${HOME}/conda/etc/profile.d/conda.sh"


# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh
conda activate base

# XXX switch once https://github.com/scikit-learn/scikit-learn/pull/29176 is merged
conda install -c conda-forge "$(get_dep conda-lock min)" -y
conda-lock install --name sklearn build_tools/github/pylatest_conda_forge_cuda_array-api_linux-64_conda.lock
