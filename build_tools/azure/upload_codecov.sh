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
CODECOV_UPLOADER_VERSION=0.4.0
CODECOV_BASE_URL="https://uploader.codecov.io/v$CODECOV_UPLOADER_VERSION"


# XXX: debug
echo "ls -la $BUILD_REPOSITORY_LOCALPATH"
ls -la $BUILD_REPOSITORY_LOCALPATH

# Check that the git repo is located at the expected location:
if [[ ! -d "$BUILD_REPOSITORY_LOCALPATH/.git" ]]; then
    echo "Could not find the git checkout at $BUILD_REPOSITORY_LOCALPATH"
    exit 1
fi
# Check that the combined coverage file exists at the expected location:
if [[ ! -f "$BUILD_REPOSITORY_LOCALPATH/.coverage" ]]; then
    echo "Could not find the combined coverage file at $BUILD_REPOSITORY_LOCALPATH/.coverage"
    exit 1
fi
if [[ $OSTYPE == *"linux"* ]]; then
    curl -Os "$CODECOV_BASE_URL/linux/codecov"
    SHA256SUM="671cf0d89d1c149f57e1a9a31f3fb567ab4209e4d5829f13ff7b8c104db7131f  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} -r $BUILD_REPOSITORY_LOCALPATH -f $BUILD_REPOSITORY_LOCALPATH/.coverage -Z
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os "$CODECOV_BASE_URL/macos/codecov"
    SHA256SUM="7549819f0fe115e113ec3538e259d748e87d84f68afa5deadc798967ec716b8d  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} -r $BUILD_REPOSITORY_LOCALPATH -f $BUILD_REPOSITORY_LOCALPATH/.coverage -Z
else
    curl -Os "$CODECOV_BASE_URL/windows/codecov.exe"
    SHA256SUM="15fb34be4eb9949ad4e964a0e21c4efc79657de05b2c799e041d7293dccf60eb codecov.exe"
    echo "$SHA256SUM" | sha256sum -c
    ./codecov.exe -t ${CODECOV_TOKEN} -r $BUILD_REPOSITORY_LOCALPATH -f $BUILD_REPOSITORY_LOCALPATH/.coverage -Z
fi
