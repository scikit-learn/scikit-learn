#!/bin/bash

set -e
set -x

wget $GITHUB_ARTIFACT_URL
unzip doc*.zip -d doc/_build/html/stable
