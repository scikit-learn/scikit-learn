#!/bin/bash

set -e
set -x

wget $GITHUB_ARTIFACT_URL
mkdir -p doc/_build/html/stable
unzip doc*.zip -d doc/_build/html/stable
