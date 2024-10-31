#!/bin/bash

set -e
set -x

if [[ -z "${CONDA}" ]]; then
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)

    export CONDA=$MAMBA_ROOT_PREFIX
    set +x
    echo "##vso[task.setvariable variable=CONDA]$CONDA"
    set -x
else
    sudo chown -R $USER $CONDA
fi

set +x
echo "##vso[task.prependpath]$CONDA/bin"
set -x
