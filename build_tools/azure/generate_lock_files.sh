#!/bin/bash

set -e

build_name=pylatest_conda_forge_mkl
platforms="linux-64 osx-64"

for platform in ${platforms}; do
    conda-lock \
        --platform ${platform} \
        --file ${build_name}_${platform}_environment.yml \
        --filename-template "${build_name}_${platform}.lock"
done
