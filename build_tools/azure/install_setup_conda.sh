#!/bin/bash

set -e
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
    echo "##vso[task.setvariable variable=CONDA]$HOME/miniforge3"
else
    sudo chown -R $USER $CONDA
fi

echo "##vso[task.prependpath]$CONDA/bin"
