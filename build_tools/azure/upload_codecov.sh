#!/bin/bash

set -e

# called when COVERAGE=="true" and DISTRIB=="conda"
export PATH=$HOME/miniconda3/bin:$PATH
source activate $VIRTUALENV

# Need to run codecov from a git checkout, so we copy .coverage
# from TEST_DIR where pytest has been run
pushd $TEST_DIR
coverage combine --append
popd
cp $TEST_DIR/.coverage $BUILD_REPOSITORY_LOCALPATH

if [[ $OSTYPE == *"linux"* ]]; then
    curl -Os https://uploader.codecov.io/v0.3.5/linux/codecov
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os https://uploader.codecov.io/v0.3.5/macos/codecov
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
else
    curl -Os https://uploader.codecov.io/v0.3.5/windows/codecov.exe
    ./codecov.exe -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
fi
