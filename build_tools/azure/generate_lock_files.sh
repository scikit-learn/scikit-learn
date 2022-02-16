#!/bin/bash

set -e

function conda_lock_func() {
    build_name=$1
    platform=$2
    environment_base_name=${3:-${build_name}}
    conda-lock lock \
        --platform ${platform} \
        --file ${environment_base_name}_environment.yml \
        --filename-template "${build_name}_${platform}.lock"
}

build_name=pylatest_conda_forge_mkl
platforms="linux-64 osx-64"

for platform in ${platforms}; do
    conda_lock_func ${build_name} ${platform} ${build_name}_${platform}
done

build_name=py37_conda_defaults_openblas
conda_lock_func ${build_name} linux-64
