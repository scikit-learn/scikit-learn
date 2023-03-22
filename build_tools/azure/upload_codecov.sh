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

# When we update the codecov uploader version, we need to update the checksums.
# The checksum for each codecov binary is available at
# https://uploader.codecov.io e.g. for linux
# https://uploader.codecov.io/v0.4.0/linux/codecov.SHA256SUM. In principle we
# need to check the signatures with the codecov gpg key as well, see
# https://docs.codecov.com/docs/codecov-uploader#integrity-checking-the-uploader
# for more details
CODECOV_UPLOADER_VERSION=0.4.0
CODECOV_BASE_URL="https://uploader.codecov.io/v$CODECOV_UPLOADER_VERSION"
if [[ $OSTYPE == *"linux"* ]]; then
    curl -Os "$CODECOV_BASE_URL/linux/codecov"
    SHA256SUM="671cf0d89d1c149f57e1a9a31f3fb567ab4209e4d5829f13ff7b8c104db7131f  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os "$CODECOV_BASE_URL/macos/codecov"
    SHA256SUM="7549819f0fe115e113ec3538e259d748e87d84f68afa5deadc798967ec716b8d  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
else
    curl -Os "$CODECOV_BASE_URL/windows/codecov.exe"
    SHA256SUM="15fb34be4eb9949ad4e964a0e21c4efc79657de05b2c799e041d7293dccf60eb  codecov.exe"
    echo "$SHA256SUM" | sha256sum -c
    ./codecov.exe -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
fi
