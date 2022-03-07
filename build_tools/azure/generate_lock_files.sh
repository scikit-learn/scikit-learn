#!/bin/bash

set -e

function conda_lock_func() {
    build_name=$1
    platform=$2
    environment_base_name=${3:-${build_name}}
    conda-lock lock \
        --mamba \
        --platform ${platform} \
        --file ${environment_base_name}_environment.yml \
        --lockfile "${build_name}_${platform}_conda-lock.yml"
}

function pip_lock_func() {
    build_name=$1
    pip-compile "${build_name}_requirements.txt" -o "${build_name}_lock.txt"
}

build_name=pylatest_conda_forge_mkl
platforms="linux-64 osx-64"

for platform in ${platforms}; do
    conda_lock_func ${build_name} ${platform} ${build_name}_${platform}
done

build_name=py38_conda_defaults_openblas
conda_lock_func ${build_name} linux-64

build_name=pylatest_pip_openblas_pandas
conda_lock_func ${build_name} linux-64

# TODO: how to make sure that python 3.8 (Ubuntu 20.04 has python 3.8.10)is
# used for this command here
build_name=ubuntu_atlas
pip_lock_func ${build_name}
