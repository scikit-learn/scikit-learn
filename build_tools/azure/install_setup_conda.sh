#!/bin/bash

set -e
set -x

if [[ -z "${CONDA}" ]]; then
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)

    export CONDA="$HOME/micromamba"
    echo "##vso[task.setvariable variable=CONDA]$CONDA"
    micromamba shell init
else
    sudo chown -R $USER $CONDA
fi

echo "##vso[task.prependpath]$CONDA/bin"
