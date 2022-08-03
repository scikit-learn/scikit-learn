#!/bin/bash

# This script is meant to be called by the "script" step defined
# in the ".travis.yml" file. While this step is forbidden by the
# continuous deployment jobs, we have to execute the scripts for
# testing the continuous integration jobs.

if [[ $BUILD_WHEEL != true ]]; then
    # This trick will make Travis terminate the continuation of the pipeline
    bash build_tools/travis/test_script.sh || travis_terminate 1
    bash build_tools/travis/test_docs.sh || travis_terminate 1
fi
