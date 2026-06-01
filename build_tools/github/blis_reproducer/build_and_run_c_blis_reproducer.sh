#!/usr/bin/env bash
set -euo pipefail

workdir="${RUNNER_TEMP:-/tmp}/blis-c-repro"
rm -rf "$workdir"
mkdir -p "$workdir"
cd "$workdir"

curl -L -o blis.tar.gz https://github.com/flame/blis/archive/refs/tags/1.1.tar.gz
tar -xzf blis.tar.gz
cd blis-1.1

./configure auto --enable-cblas
make -j"$(sysctl -n hw.ncpu)"

lib_dir="$(find "$workdir/blis-1.1/lib" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
if [[ -z "${lib_dir}" ]]; then
  echo "Unable to locate BLIS library directory" >&2
  exit 1
fi

arch_name="$(basename "$lib_dir")"
include_dir="$workdir/blis-1.1/include/$arch_name"
if [[ ! -d "$include_dir" ]]; then
  include_dir="$workdir/blis-1.1/include"
fi

cd /tmp/workspace/scikit-learn/scikit-learn
cc \
  -I"$include_dir" \
  -I"$workdir/blis-1.1/include" \
  build_tools/github/blis_reproducer/blis_gemm_reproducer.c \
  -L"$lib_dir" \
  -lblis -lm \
  -Wl,-rpath,"$lib_dir" \
  -o "$workdir/blis_gemm_reproducer"

export BLIS_NUM_THREADS=8
"$workdir/blis_gemm_reproducer"
