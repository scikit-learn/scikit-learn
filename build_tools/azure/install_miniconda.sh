#!/bin/bash

set -v -e

# Install Miniconda
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    if [[ "$BITS32" == "yes" ]]; then
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda.sh
    else
        wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    fi
elif [[ "$unamestr" == 'Darwin' ]]; then
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
else
  echo Error
fi
chmod +x miniconda.sh
./miniconda.sh -b
