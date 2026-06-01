#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <blis|openblas|newaccelerate>"
  exit 2
fi

blas_impl="$1"
env_name="blis-repro-${blas_impl}"

if command -v mamba >/dev/null 2>&1; then
  mamba_bin="mamba"
elif command -v micromamba >/dev/null 2>&1; then
  mamba_bin="micromamba"
else
  echo "Neither mamba nor micromamba is available." >&2
  exit 1
fi

"$mamba_bin" create -y -n "$env_name" numpy "libblas=*=*_${blas_impl}" >/tmp/mamba-${env_name}.log 2>&1

# Make runs deterministic and stress threaded BLAS behavior.
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

if [[ "$blas_impl" == "blis" ]]; then
  export BLIS_NUM_THREADS=8
  unset OPENBLAS_NUM_THREADS
else
  export OPENBLAS_NUM_THREADS=8
  unset BLIS_NUM_THREADS
fi

"$mamba_bin" run -n "$env_name" python build_tools/github/blis_reproducer/numpy_blis_reproducer.py
