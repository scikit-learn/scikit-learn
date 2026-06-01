# BLIS macOS arm64 reproducers

This folder contains two independent reproducers for the macOS arm64 BLAS regressions:

- `numpy_blis_reproducer.py`: numpy-only reproducer that compares a GEMM-based pairwise distance computation against a pure numpy reference implementation.
- `blis_gemm_reproducer.c`: pure C reproducer that calls `cblas_dgemm` and validates results against a scalar long-double reference loop.

## Run the numpy reproducer in dedicated BLAS envs

```bash
./build_tools/github/blis_reproducer/run_numpy_blas_reproducer.sh blis
./build_tools/github/blis_reproducer/run_numpy_blas_reproducer.sh openblas
./build_tools/github/blis_reproducer/run_numpy_blas_reproducer.sh newaccelerate
```

## Run the pure C reproducer with locally built BLIS

```bash
./build_tools/github/blis_reproducer/build_and_run_c_blis_reproducer.sh
```

The CI workflow `.github/workflows/blis-macos-arm64-reproducer.yml` runs both reproducers on macOS arm64.
