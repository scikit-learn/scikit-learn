#!/bin/bash

set -e

conda-lock \
    --platform linux-64 \
    --platform osx-64 \
    --file pylatest_conda_forge_mkl_environment.yml \
    --filename-template "pylatest_conda_forge_mkl_{platform}.lock"
