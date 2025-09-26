#!/bin/bash

set -e

"${SHELL}" <(curl -Ls micro.mamba.pm/install.sh) < /dev/null

export MAMBA_EXE='/root/.local/bin/micromamba';
$MAMBA_EXE shell init -s bash

micromamba env create -f build_tools/circle/doc_environment.yml -n sklearn-dev --yes
# Auto-activate sklearn-dev in terminal
echo "micromamba activate sklearn-dev" >> $HOME/.bashrc

# Enables users to activate environment without having to specify the full path
# echo "envs_dirs:
#  - /home/codespace/micromamba/envs" > /opt/conda/.condarc
