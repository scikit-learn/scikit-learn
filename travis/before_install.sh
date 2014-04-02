#!/bin/bash

set -e

if [[ "$ENV" -eq "conda" ]]; then
    # Deactivate the travis-provided virtual environment and setup a
    # conda-based environment instead
    deactivate
    sudo apt-get install python-pip
    rm -rf /home/travis/condaenv
    sudo pip install conda
    sudo conda init
    if [[ $NUMPY_VERSION != '' ]]; then
        # Oldest version of numpy and scipy supported by conda and sklearn
        conda create -p /home/travis/condaenv --yes \
            python=$PYTHON_VERSION numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION \
            pip nose
    else
        # Latest supported version by conda
        conda create -p /home/travis/condaenv --yes \
            python=$PYTHON_VERSION numpy scipy pip nose
    fi
    export PATH=/home/travis/condaenv/bin:$PATH

elif [[ "$ENV" -eq "ubuntu" ]]; then
    # Test with standard ubuntu packages
    sudo apt-get update -qq
    sudo apt-get install -qq python-scipy python-nose
fi
