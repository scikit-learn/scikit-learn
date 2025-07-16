#!/bin/bash

set -e
set -x

PLATFORM=$(uname)
if [[ "$PLATFORM" =~ MINGW|MSYS ]]; then
    PLATFORM=Windows
fi
if [[ "$PLATFORM" == "Windows" ]]; then
    EXTENSION="exe"
else
    EXTENSION="sh"
fi
INSTALLER="miniforge.$EXTENSION"
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$PLATFORM-$(uname -m).$EXTENSION"
curl -L ${MINIFORGE_URL} -o "$INSTALLER"

MINICONDA_DIR="$HOME/miniforge3"
if [[ "$PLATFORM" != "Windows" ]]; then
    bash "$INSTALLER" -b -u -p $MINICONDA_DIR
else
    WIN_MINICONDA_DIR=$(cygpath -w "$MINICONDA_DIR")
    cmd "/C $INSTALLER /InstallationType=JustMe /RegisterPython=0 /S /D=$WIN_MINICONDA_DIR"
fi

# Add conda to the PATH so that it can be used in further Azure CI steps.
# Need set +x for ##vso Azure magic otherwise it may add a quote in the PATH.
# For more details, see https://github.com/microsoft/azure-pipelines-tasks/issues/10331
set +x
if [[ "$PLATFORM" == "Windows" ]]; then
   echo "##vso[task.prependpath]$MINICONDA_DIR/Scripts"
else
   echo "##vso[task.prependpath]$MINICONDA_DIR/bin"
fi
set -x
