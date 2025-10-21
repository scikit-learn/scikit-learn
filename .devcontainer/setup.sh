#!/bin/bash

set -e

"${SHELL}" <(curl -Ls micro.mamba.pm/install.sh) < /dev/null
# .bashrc has been updated by the mamba install one-liner above.
# 'source $HOME/.bashrc' sets up micromamba for later use
source $HOME/.bashrc

micromamba env create -f build_tools/circle/doc_environment.yml -n sklearn-dev --yes
# Install additional packages:
# - ipykernel: to be able to use the VS Code Jupyter integration
# - pre-commit: avoid linting issues
micromamba install pre-commit ipykernel -n sklearn-dev --yes
# install pre-commit hooks
micromamba activate sklearn-dev
pre-commit install

# Auto-activate sklearn-dev in terminal
echo "micromamba activate sklearn-dev" >> $HOME/.bashrc
