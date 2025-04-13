#!/bin/bash

set -e
set -x

if [[ -z "${CONDA}" ]]; then
    # In some runners (macOS-13 and macOS-14 in October 2024) conda is not
    # installed so we install it ourselves
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    wget ${MINIFORGE_URL} -O miniforge.sh
    bash miniforge.sh -b -u -p $HOME/miniforge3
    CONDA="$HOME/miniforge3"
else
    # In most runners (in October 2024) conda is installed,
    # but in a system folder and we want it user writable
    sudo chown -R $USER $CONDA
fi

# Add conda to the PATH so that it can be used in further Azure CI steps.
# Need set +x for ##vso Azure magic otherwise it may add a quote in the PATH.
# For more details, see https://github.com/microsoft/azure-pipelines-tasks/issues/10331
set +x
echo "##vso[task.prependpath]$CONDA/bin"
set -x
