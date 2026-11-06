"""
Benchmark of `sandwich_dot(X, W) = X.T @ diag(W) @ X` for dense X.

Compare the previous GEMM implementation with the current one in
`sklearn.linear_model._linear_loss`, which uses BLAS SYRK via `A.T @ A`.

Usage:
    python bench_sandwich_dot.py [--n-threads N] [--dtype float64 float32]
"""

import argparse
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from sklearn.linear_model._linear_loss import sandwich_dot


def sandwich_dot_gemm(X, W):
    """Previous implementation, a full GEMM."""
    return X.T @ (W[:, None] * X)


def best_time(func, *args, n_repeat=5):
    func(*args)  # warm-up
    times = []
    for _ in range(n_repeat):
        start = perf_counter()
        func(*args)
        times.append(perf_counter() - start)
    return min(times)


def run(dtypes, shapes):
    info = threadpool_info()
    if info:
        print(f"BLAS: {info[0]['internal_api']}, threads: {info[0]['num_threads']}")
    rng = np.random.default_rng(0)
    header = f"{'dtype':8} {'n_samples':>10} {'n_features':>10} "
    header += f"{'gemm [ms]':>10} {'new [ms]':>10} {'speedup':>8} {'rel. diff':>10}"
    print(header)
    for dtype in dtypes:
        for n_samples, n_features in shapes:
            X = rng.standard_normal((n_samples, n_features)).astype(dtype)
            W = rng.uniform(0.1, 1, n_samples).astype(dtype)
            expected = sandwich_dot_gemm(X, W)
            diff = np.max(np.abs(sandwich_dot(X, W) - expected)) / np.max(
                np.abs(expected)
            )
            t_gemm = best_time(sandwich_dot_gemm, X, W)
            t_new = best_time(sandwich_dot, X, W)
            print(
                f"{dtype:8} {n_samples:>10} {n_features:>10} "
                f"{t_gemm * 1e3:>10.1f} {t_new * 1e3:>10.1f} "
                f"{t_gemm / t_new:>8.2f} {diff:>10.1e}",
                flush=True,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-threads", type=int, default=None, help="Limit BLAS threads."
    )
    parser.add_argument(
        "--dtype", nargs="+", default=["float64", "float32"], help="dtypes to use."
    )
    args = parser.parse_args()
    shapes = [
        (100_000, 20),
        (100_000, 100),
        (100_000, 500),
        (1_000_000, 50),
        (1_000_000, 200),
        (20_000, 1000),
    ]
    with threadpool_limits(limits=args.n_threads):
        run(args.dtype, shapes)
