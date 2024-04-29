#!/bin/bash

set -e

# Do not upload to codecov on forks
if [[ "$BUILD_REPOSITORY_NAME" != "scikit-learn/scikit-learn" ]]; then
    exit 0
fi

# When we update the codecov uploader version, we need to update the checksums.
# The checksum for each codecov binary is available at
# https://uploader.codecov.io e.g. for linux
# https://uploader.codecov.io/v0.7.1/linux/codecov.SHA256SUM.

# Instead of hardcoding a specific version and signature in this script, it
# would be possible to use the "latest" symlink URL but then we need to
# download both the codecov.SHA256SUM files each time and check the signatures
# with the codecov gpg key as well, see:
# https://docs.codecov.com/docs/codecov-uploader#integrity-checking-the-uploader
# However this approach would yield a larger number of downloads from
# codecov.io and keybase.io, therefore increasing the risk of running into
# network failures.
CODECOV_UPLOADER_VERSION=0.7.1
CODECOV_BASE_URL="https://uploader.codecov.io/v$CODECOV_UPLOADER_VERSION"


# Check that the git repo is located at the expected location:
if [[ ! -d "$BUILD_REPOSITORY_LOCALPATH/.git" ]]; then
    echo "Could not find the git checkout at $BUILD_REPOSITORY_LOCALPATH"
    exit 1
fi

# Check that the combined coverage file exists at the expected location:
export COVERAGE_XML="$BUILD_REPOSITORY_LOCALPATH/coverage.xml"
if [[ ! -f "$COVERAGE_XML" ]]; then
    echo "Could not find the combined coverage file at $COVERAGE_XML"
    exit 1
fi

if [[ $OSTYPE == *"linux"* ]]; then
    curl -Os "$CODECOV_BASE_URL/linux/codecov"
    SHA256SUM="b9282b8b43eef83f722646d8992c4dd36563046afe0806722184e7e9923a6d7b  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} -R $BUILD_REPOSITORY_LOCALPATH -f coverage.xml -Z --verbose
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os "$CODECOV_BASE_URL/macos/codecov"
    SHA256SUM="e4ce34c144d3195eccb7f8b9ca8de092d2a4be114d927ca942500f3a6326225c  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov -t ${CODECOV_TOKEN} -R $BUILD_REPOSITORY_LOCALPATH -f coverage.xml -Z --verbose
else
    curl -Os "$CODECOV_BASE_URL/windows/codecov.exe"
    SHA256SUM="f5de88026f061ff08b88a5895f9c11855523924ceb8174e027403dd20fa5e4d6  codecov.exe"
    echo "$SHA256SUM" | sha256sum -c
    ./codecov.exe -t ${CODECOV_TOKEN} -R $BUILD_REPOSITORY_LOCALPATH -f coverage.xml -Z --verbose
fi
