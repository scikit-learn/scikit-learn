#!/bin/bash

set -e
set -x

# Install Python 3.8 because of a bug with Python 3.9
export PATH=$CONDA/bin:$PATH
conda create -n upload -y python=3.8
source activate upload
conda install -y anaconda-client
which anaconda

# TODO: how to detect nightly builds for ARM on Circle CI?
if [ "$GITHUB_EVENT_NAME" == "schedule" ]; then
    ANACONDA_ORG="scipy-wheels-nightly"
    ANACONDA_TOKEN="$SCIKIT_LEARN_NIGHTLY_UPLOAD_TOKEN"
else
    ANACONDA_ORG="scikit-learn-wheels-staging"
    ANACONDA_TOKEN="$SCIKIT_LEARN_STAGING_UPLOAD_TOKEN"
fi

if [ "$ANACONDA_TOKEN" == "" ];
    echo "No upload token defined to upload to anaconda.org/$ANACONDA_ORG: skipping."
else
    # Force a replacement if the remote file already exists
    anaconda -t $ANACONDA_TOKEN upload --force -u $ANACONDA_ORG dist/artifact/*
    echo "Index: https://pypi.anaconda.org/$ANACONDA_ORG/simple"
fi
