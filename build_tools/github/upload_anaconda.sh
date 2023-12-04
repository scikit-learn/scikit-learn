#!/bin/bash

set -e
set -x

# Note: build_wheels.sh has the same branch (only for NumPy 2.0 transition)
if [[ "$GITHUB_EVENT_NAME" == "schedule" || "$CIRRUS_CRON" == "nightly" ]]; then
    ANACONDA_ORG="scientific-python-nightly-wheels"
    ANACONDA_TOKEN="$SCIKIT_LEARN_NIGHTLY_UPLOAD_TOKEN"
else
    ANACONDA_ORG="scikit-learn-wheels-staging"
    ANACONDA_TOKEN="$SCIKIT_LEARN_STAGING_UPLOAD_TOKEN"
fi

# Install Python 3.8 because of a bug with Python 3.9
export PATH=$CONDA/bin:$PATH
conda create -n upload -y python=3.8
source activate upload
conda install -y anaconda-client

# Force a replacement if the remote file already exists
anaconda -t $ANACONDA_TOKEN upload --force -u $ANACONDA_ORG $ARTIFACTS_PATH/*
echo "Index: https://pypi.anaconda.org/$ANACONDA_ORG/simple"
