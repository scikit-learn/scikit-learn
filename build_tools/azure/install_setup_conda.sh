#!/bin/bash

set -e
set -x

if [[ -z "${CONDA}" ]]; then
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)

    export CONDA=$MAMBA_ROOT_PREFIX
    echo "##vso[task.setvariable variable=CONDA]$CONDA"
else
    sudo chown -R $USER $CONDA
fi

echo "##vso[task.prependpath]$CONDA/bin"
