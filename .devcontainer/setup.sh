#!/bin/bash

set -e

"${SHELL}" <(curl -Ls micro.mamba.pm/install.sh) < /dev/null
# .bashrc has been updated by the mamba install one-liner above.
# 'source $HOME/.bashrc' sets up micromamba for later use
source $HOME/.bashrc

micromamba env create -f build_tools/circle/doc_environment.yml -n sklearn-dev --yes
# install ipykernel because it is handy to execute the Jupyter notebook VS Code
micromamba install ipykernel -n sklearn-dev --yes
# install pre-commit since it will avoid potential failures with the CI
micromamba install pre-commit -n sklearn-dev --yes
# Auto-activate sklearn-dev in terminal
echo "micromamba activate sklearn-dev" >> $HOME/.bashrc
