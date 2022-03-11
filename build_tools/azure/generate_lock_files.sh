#!/bin/bash

set -e

dir_name=./build_tools/azure

function conda_lock_func() {
    build_name=$1
    platform=$2
    environment_base_name=${3:-${build_name}}
    conda-lock lock \
        --mamba \
        --platform ${platform} \
        --file "${dir_name}/${environment_base_name}_environment.yml" \
        --lockfile "${dir_name}/${build_name}_${platform}_conda-lock.yml"
}

function pip_lock_func() {
    build_name=$1
    pip-compile "${dir_name}/${build_name}_requirements.txt" -o "${dir_name}/${build_name}_lock.txt"
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

build_name=pylatest_pip_openblas_pandas
conda_lock_func ${build_name} linux-64

build_name=pylatest_conda_forge_mkl_no_coverage
conda_lock_func ${build_name} linux-64

build_name=pypy3
conda_lock_func ${build_name} linux-64

build_name=pylatest_pip_scipy_dev
conda_lock_func ${build_name} linux-64

build_name=pylatest_conda_mkl_no_openmp
conda_lock_func ${build_name} osx-64

# TODO: how to make sure that python 3.8 (Ubuntu 20.04 has python 3.8.10) is
# used for this command here
build_name=ubuntu_atlas
pip_lock_func ${build_name}

# TODO: how to make sure that python 3.9 (debian-32 docker image has python
# 3.9) is used for this command here
build_name=debian_atlas_32bit
pip_lock_func ${build_name}
