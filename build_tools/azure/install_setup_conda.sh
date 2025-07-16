#!/bin/bash

set -e
set -x

PLATFORM=$(uname)
if [[ "$PLATFORM" =~ MINGW|MSYS ]]; then
    PLATFORM=Windows
fi
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$PLATFORM-$(uname -m).sh"
curl -L ${MINIFORGE_URL} -o miniforge.sh
bash miniforge.sh -b -u -p $HOME/miniforge3
CONDA="$HOME/miniforge3"

# Add conda to the PATH so that it can be used in further Azure CI steps.
# Need set +x for ##vso Azure magic otherwise it may add a quote in the PATH.
# For more details, see https://github.com/microsoft/azure-pipelines-tasks/issues/10331
set +x
if [[ "$PLATFORM" == "Windows" ]]; then
   echo "##vso[task.prependpath]$CONDA/Scripts"
else
   echo "##vso[task.prependpath]$CONDA/bin"
fi
set -x
