#!/bin/bash

set -e
set -x

if [[ "$GITHUB_EVENT_NAME" == "schedule" \
          || "$GITHUB_EVENT_NAME" == "workflow_dispatch" ]]; then
    ANACONDA_ORG="scientific-python-nightly-wheels"
    ANACONDA_TOKEN="$SCIKIT_LEARN_NIGHTLY_UPLOAD_TOKEN"
else
    ANACONDA_ORG="scikit-learn-wheels-staging"
    ANACONDA_TOKEN="$SCIKIT_LEARN_STAGING_UPLOAD_TOKEN"
fi

export PATH=$CONDA/bin:$PATH
conda create -n upload -y anaconda-client
source activate upload

# Force a replacement if the remote file already exists
anaconda -t $ANACONDA_TOKEN upload --force -u $ANACONDA_ORG $ARTIFACTS_PATH/*
echo "Index: https://pypi.anaconda.org/$ANACONDA_ORG/simple"
