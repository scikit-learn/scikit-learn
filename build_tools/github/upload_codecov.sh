#!/bin/bash

set -e

# Do not upload to codecov on forks
if [[ "$BUILD_REPOSITORY_NAME" != "scikit-learn/scikit-learn" ]]; then
    exit 0
fi

# When we update the codecov uploader version, we need to update the checksums.
# The checksum for each codecov binary is available at
# https://cli.codecov.io e.g. for linux
# https://cli.codecov.io/v10.2.1/linux/codecov.SHA256SUM.

# Instead of hardcoding a specific version and signature in this script, it
# would be possible to use the "latest" symlink URL but then we need to
# download both the codecov.SHA256SUM files each time and check the signatures
# with the codecov gpg key as well, see:
# https://docs.codecov.com/docs/codecov-uploader#integrity-checking-the-uploader
# However this approach would yield a larger number of downloads from
# codecov.io and keybase.io, therefore increasing the risk of running into
# network failures.
CODECOV_CLI_VERSION=10.2.1
CODECOV_BASE_URL="https://cli.codecov.io/v$CODECOV_CLI_VERSION"

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
    SHA256SUM="39dd112393680356daf701c07f375303aef5de62f06fc80b466b5c3571336014  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov upload-coverage -t ${CODECOV_TOKEN} -f coverage.xml -Z
    ./codecov do-upload --disable-search --report-type test_results --file $JUNIT_FILE
elif [[ $OSTYPE == *"darwin"* ]]; then
    curl -Os "$CODECOV_BASE_URL/macos/codecov"
    SHA256SUM="01183f6367c7baff4947cce389eaa511b7a6d938e37ae579b08a86b51f769fd9  codecov"
    echo "$SHA256SUM" | shasum -a256 -c
    chmod +x codecov
    ./codecov upload-coverage -t ${CODECOV_TOKEN} -f coverage.xml -Z
    ./codecov do-upload --disable-search --report-type test_results --file $JUNIT_FILE
else
    curl -Os "$CODECOV_BASE_URL/windows/codecov.exe"
    SHA256SUM="e54e9520428701a510ef451001db56b56fb17f9b0484a266f184b73dd27b77e7  codecov.exe"
    echo "$SHA256SUM" | sha256sum -c
    ./codecov.exe upload-coverage -t ${CODECOV_TOKEN} -f coverage.xml -Z
    ./codecov.exe do-upload --disable-search --report-type test_results --file $JUNIT_FILE
fi
