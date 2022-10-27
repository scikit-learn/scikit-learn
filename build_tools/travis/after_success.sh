#!/bin/bash

# This script is meant to be called by the "after_success" step
# defined in ".travis.yml". In particular, we upload the wheels
# of the ARM64 architecture for the continuous deployment jobs.

set -e

# The wheels cannot be uploaded on PRs
if [[ $BUILD_WHEEL == true && $TRAVIS_EVENT_TYPE != pull_request ]]; then
    # Nightly upload token and staging upload token are set in
    # Travis settings (originally generated at Anaconda cloud)
    if [[ $TRAVIS_EVENT_TYPE == cron ]]; then
        ANACONDA_ORG="scipy-wheels-nightly"
        ANACONDA_TOKEN="$SCIKIT_LEARN_NIGHTLY_UPLOAD_TOKEN"
    else
        ANACONDA_ORG="scikit-learn-wheels-staging"
        ANACONDA_TOKEN="$SCIKIT_LEARN_STAGING_UPLOAD_TOKEN"
    fi

    MINICONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"
    wget $MINICONDA_URL -O miniconda.sh
    MINICONDA_PATH=$HOME/miniconda
    chmod +x miniconda.sh && ./miniconda.sh -b -p $MINICONDA_PATH

    # Install Python 3.8 because of a bug with Python 3.9
    export PATH=$MINICONDA_PATH/bin:$PATH
    conda create -n upload -y python=3.8
    source activate upload
    conda install -y anaconda-client

    # Force a replacement if the remote file already exists
    anaconda -t $ANACONDA_TOKEN upload --force -u $ANACONDA_ORG wheelhouse/*.whl
    echo "Index: https://pypi.anaconda.org/$ANACONDA_ORG/simple"
fi
