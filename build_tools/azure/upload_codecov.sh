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

# The checksum for each codecov binary is available at
# https://uploader.codecov.io
# When updating the binary, we also need to update the checksum
if [[ $OSTYPE == *"linux"* ]]; then
    curl -Os https://uploader.codecov.io/v0.3.5/linux/codecov
    SHA256SUM="080b43eaec3434326bb0f61653a82d27aba15c311ddde9d3f68cb364314f7aae  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os https://uploader.codecov.io/v0.3.5/macos/codecov
    SHA256SUM="dfd7b0e3b245967477933c7f0c2b9f8f11dce775bc121a44de89b2b8b04251cd  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
else
    curl -Os https://uploader.codecov.io/v0.3.5/windows/codecov.exe
    SHA256SUM="b3a5971ac1d5c2b50476bb9f1079e6e6fc967c4f0d1ede9a9d0ae2311d5d3dad codecov.exe"
    echo "$SHA256SUM" | shasum -a256 -c
    ./codecov.exe -t ${CODECOV_TOKEN} --rootDir $BUILD_REPOSITORY_LOCALPATH
fi
