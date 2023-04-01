#!/bin/bash

set -e

# When we update the codecov uploader version, we need to update the checksums.
# The checksum for each codecov binary is available at
# https://uploader.codecov.io e.g. for linux
# https://uploader.codecov.io/v0.4.0/linux/codecov.SHA256SUM.

# Instead of hardcoding a specific version and signature in this script, it
# would be possible to use the "latest" symlink URL but then we need to
# download both the codecov.SHA256SUM files each time and check the signatures
# with the codecov gpg key as well, see:
# https://docs.codecov.com/docs/codecov-uploader#integrity-checking-the-uploader
# However this approach would yield a larger number of downloads from
# codecov.io and keybase.io, therefore increasing the risk of running into
# network failures.
CODECOV_UPLOADER_VERSION=0.4.1
CODECOV_BASE_URL="https://uploader.codecov.io/v$CODECOV_UPLOADER_VERSION"
echo "About to run codecov for repo at $BUILD_REPOSITORY_LOCALPATH and searching for coverage files in $TEST_DIR:"
if [[ $OSTYPE == *"linux"* ]]; then
    curl -Os "$CODECOV_BASE_URL/linux/codecov"
    SHA256SUM="32cb14b5f3aaacd67f4c1ff55d82f037d3cd10c8e7b69c051f27391d2e66e15c  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} --r $BUILD_REPOSITORY_LOCALPATH -s $TEST_DIR -Z
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os "$CODECOV_BASE_URL/macos/codecov"
    SHA256SUM="4ab0f06f06e9c4d25464f155b0aff36bfc1e8dbcdb19bfffd586beed1269f3af  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} -r $BUILD_REPOSITORY_LOCALPATH -s $TEST_DIR -Z
else
    curl -Os "$CODECOV_BASE_URL/windows/codecov.exe"
    SHA256SUM="e0cda212aeaebe695509ce8fa2d608760ff70bc932003f544f1ad368ac5450a8 codecov.exe"
    echo "$SHA256SUM" | sha256sum -c
    ./codecov.exe -t ${CODECOV_TOKEN} -r $BUILD_REPOSITORY_LOCALPATH -s $TEST_DIR -Z
fi
