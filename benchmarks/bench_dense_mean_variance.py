"""Benchmark the single-pass dense per-feature mean and variance.

Compares `_dense_mean_variance_axis0` against `np.mean` alone (a lower bound:
one pass, no variance), `np.mean` + `np.var` and `_incremental_mean_and_var`,
in terms of runtime and numerical error (w.r.t. a long double two-pass
reference computed on well-conditioned and large-offset data).
"""

import argparse
from time import perf_counter

import numpy as np
from sklearn.utils._mean_variance import _dense_mean_variance_axis0

from sklearn.utils.extmath import _incremental_mean_and_var


def np_mean(X):
    return np.mean(X, axis=0), None


def np_mean_var(X):
    return np.mean(X, axis=0), np.var(X, axis=0)


def incremental(X):
    mean, var, _ = _incremental_mean_and_var(X, 0.0, 0.0, 0)
    return mean, var


METHODS = {
    "cython (this)": _dense_mean_variance_axis0,
    "np.mean only": np_mean,
    "np.mean+np.var": np_mean_var,
    "_incremental_mean_and_var": incremental,
}


def best_time(func, X, repeat):
    times = []
    for _ in range(repeat):
        tic = perf_counter()
        func(X)
        times.append(perf_counter() - tic)
    return min(times)


def reference(X):
    X = X.astype(np.longdouble)
    mean = X.mean(axis=0)
    return mean, ((X - mean) ** 2).mean(axis=0)


def max_rel_err(value, ref, scale=None):
    scale = ref if scale is None else scale
    return float(np.max(np.abs((value - ref) / scale)))


def bench_runtime(shapes, repeat):
    print("\n## Runtime (best of %d, GB/s = size of X / time)\n" % repeat)
    print(f"| {'shape':>16} | dtype   | order | " + " | ".join(METHODS) + " |")
    print("|" + "---|" * (3 + len(METHODS)))
    rng = np.random.default_rng(0)
    for shape in shapes:
        X64 = rng.standard_normal(shape)
        for dtype in (np.float32, np.float64):
            for order in "CF":
                X = np.asarray(X64, dtype=dtype, order=order)
                cells = []
                for func in METHODS.values():
                    t = best_time(func, X, repeat)
                    cells.append(f"{t * 1e3:8.1f}ms ({X.nbytes / t / 1e9:4.1f})")
                print(
                    f"| {shape!s:>16} | {np.dtype(dtype).name} | {order}     | "
                    + " | ".join(cells)
                    + " |"
                )
        del X64, X


def bench_accuracy(n_samples, n_features):
    print(
        "\n## Max error w.r.t. a long double reference: |mean - ref| / std and"
        " |var - ref| / var\n"
    )
    rng = np.random.default_rng(0)
    noise = rng.standard_normal((n_samples, n_features))
    methods = {k: v for k, v in METHODS.items() if k != "np.mean only"}
    print(
        "| data | dtype | order | " + " | ".join(f"{m} (mean / var)" for m in methods)
    )
    print("|" + "---|" * (3 + len(methods)))
    cases = [
        ("N(0, 1)", np.float64, 0.0),
        ("N(0, 1)", np.float32, 0.0),
        ("1e3 + N(0, 1)", np.float32, 1e3),
        ("1e8 + N(0, 1)", np.float64, 1e8),
        ("1e10 + N(0, 1)", np.float64, 1e10),
    ]
    for name, dtype, offset in cases:
        for order in "CF":
            X = np.asarray(offset + noise, dtype=dtype, order=order)
            ref_mean, ref_var = reference(X)
            cells = []
            for func in methods.values():
                mean, var = func(X)
                cells.append(
                    f"{max_rel_err(mean, ref_mean, np.sqrt(ref_var)):.1e} / "
                    f"{max_rel_err(var, ref_var):.1e}"
                )
            print(
                f"| {name} | {np.dtype(dtype).name} | {order} | "
                + " | ".join(cells)
                + " |"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    if args.quick:
        shapes = [(100_000, 50), (2_000, 5_000)]
    else:
        shapes = [
            (1_000_000, 100),  # large, memory bandwidth bound
            (10_000_000, 5),  # tall and skinny
            (5_000, 20_000),  # wide
            (20_000, 50),  # fits in cache
        ]
    bench_runtime(shapes, args.repeat)
    bench_accuracy(200_000, 20)
