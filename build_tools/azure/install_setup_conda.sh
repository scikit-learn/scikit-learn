#!/bin/bash

set -e
# make sure to not use `set -x` together with `##vso` because it will append
# single quotes to the output
set +x

if [[ -z "${CONDA}" ]]; then
    # Determine architecture and OS
    ARCH=$(uname -m)
    OS=$(uname -s)

    if [[ "$ARCH" == "arm64" ]]; then
        MINIFORGE_ARCH="arm64"
    else
        MINIFORGE_ARCH="x86_64"
    fi

    if [[ "$OS" == "Darwin" ]]; then
        MINIFORGE_OS="MacOSX"
    else
        MINIFORGE_OS="Linux"
    fi

    # Download and install miniforge
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${MINIFORGE_OS}-${MINIFORGE_ARCH}.sh"
    wget ${MINIFORGE_URL} -O miniforge.sh
    chmod +x miniforge.sh
    ./miniforge.sh -b -u -p $HOME/miniforge3

    # Set CONDA environment variable to point to miniforge installation
    # We need to export because `##vso` does not create the variable in the
    # current shell
    export CONDA="$HOME/miniforge3"
    echo "##vso[task.setvariable variable=CONDA]$CONDA"
else
    sudo chown -R $USER $CONDA
fi

echo "##vso[task.prependpath]$CONDA/bin"
