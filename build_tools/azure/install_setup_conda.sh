#!/bin/bash

set -e
set -x

if [[ -z "${CONDA}" ]]; then
    # By default, it uses conda-forge channels and init conda
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)

    # Set CONDA environment variable to point to micromamba installation
    echo "##vso[task.setvariable variable=CONDA]$HOME/micromamba"
fi

# Add conda to PATH
echo "##vso[task.prependpath]$CONDA/bin"

if [[ -n "${CONDA}" ]]; then
    # Take ownership of conda installation
    sudo chown -R $USER $CONDA
fi
