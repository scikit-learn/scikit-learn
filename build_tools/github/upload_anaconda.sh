#!/bin/bash

set -e
set -x

if [ "$GITHUB_EVENT_NAME" == "schedule" ]; then
    ANACONDA_ORG="scipy-wheels-nightly"
    ANACONDA_TOKEN="$SCIKIT_LEARN_NIGHTLY_UPLOAD_TOKEN"
else
    ANACONDA_ORG="scikit-learn-wheels-staging"
    ANACONDA_TOKEN="$SCIKIT_LEARN_STAGING_UPLOAD_TOKEN"
fi

pip install git+https://github.com/Anaconda-Server/anaconda-client

# Force a replacement if the remote file already exists
anaconda -t $ANACONDA_TOKEN upload --force -u $ANACONDA_ORG dist/artifact/*
echo "Index: https://pypi.anaconda.org/$ANACONDA_ORG/simple"
