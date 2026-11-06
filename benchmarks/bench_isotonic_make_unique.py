"""Benchmark Cython _make_unique vs a NumPy-only replacement."""

import time

import numpy as np
from sklearn._isotonic import _make_unique as _make_unique_cy


def _make_unique_np(X, y, sample_weights):
    """Vectorized approximation: groups by adjacent diff, not by group anchor.

    Diverges from the Cython version when a chain of consecutive values is
    each < eps apart from its neighbor, but the cumulative drift from the
    first value of the chain is >= eps (the Cython version starts a new
    group there, this version does not).
    """
    X = np.asarray(X)
    y = np.asarray(y)
    w = np.asarray(sample_weights)
    eps = np.finfo(X.dtype).resolution
    starts = np.flatnonzero(np.concatenate(([True], np.diff(X) >= eps)))
    w_out = np.add.reduceat(w, starts)
    y_out = np.add.reduceat(y * w, starts) / w_out
    return X[starts], y_out.astype(X.dtype, copy=False), w_out


def _make_unique_np_exact(X, y, sample_weights):
    """NumPy replacement that exactly replicates the Cython anchor-based
    grouping (new group when x - group_anchor >= eps), via a Python loop
    over groups (not over all samples) using searchsorted to jump ahead.
    """
    X = np.asarray(X)
    y = np.asarray(y)
    w = np.asarray(sample_weights)
    eps = np.finfo(X.dtype).resolution
    n = len(X)
    starts = []
    i = 0
    while i < n:
        starts.append(i)
        threshold = X[i] + eps
        # first index j (relative, searching from i) with X[j] >= threshold
        j = i + np.searchsorted(X[i:], threshold, side="left")
        i = j if j > i else i + 1
    starts = np.asarray(starts)
    w_out = np.add.reduceat(w, starts)
    y_out = np.add.reduceat(y * w, starts) / w_out
    return X[starts], y_out.astype(X.dtype, copy=False), w_out


def make_data(n, dtype, dup_fraction, rng):
    """dup_fraction in [0, 1): fraction of 'collapsed' unique values."""
    n_unique = max(1, int(n * (1 - dup_fraction)))
    unique_x = np.sort(rng.uniform(0, 1000, size=n_unique)).astype(dtype)
    idx = np.sort(rng.integers(0, n_unique, size=n))
    X = unique_x[idx]
    y = rng.normal(size=n).astype(dtype)
    w = rng.uniform(0.5, 1.5, size=n).astype(dtype)
    return X, y, w


def time_it(func, args, repeat=7, number=3):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        for _ in range(number):
            func(*args)
        dt = (time.perf_counter() - t0) / number
        best = min(best, dt)
    return best


def check_close(X, y, w, X2, y2, w2):
    return (
        X.shape == X2.shape
        and np.allclose(X, X2)
        and np.allclose(y, y2, rtol=1e-5, atol=1e-6)
        and np.allclose(w, w2, rtol=1e-5, atol=1e-6)
    )


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    sizes = [1_000, 100_000, 10_000_000]
    dtypes = [np.float32, np.float64]
    dup_fractions = [0.0, 0.5, 0.99]

    header = (
        f"{'n':>10} {'dtype':>8} {'dup_frac':>8} "
        f"{'cy (ms)':>10} {'np (ms)':>10} {'speedup(cy/np)':>15} {'match':>6}"
    )
    print(header)
    print("-" * len(header))

    for dtype in dtypes:
        for n in sizes:
            for dup_fraction in dup_fractions:
                X, y, w = make_data(n, dtype, dup_fraction, rng)

                number = 3 if n <= 1_000_000 else 1
                repeat = 7 if n <= 1_000_000 else 3

                t_cy = time_it(_make_unique_cy, (X, y, w), repeat=repeat, number=number)
                t_np = time_it(_make_unique_np, (X, y, w), repeat=repeat, number=number)

                out_cy = _make_unique_cy(X, y, w)
                out_np = _make_unique_np(X, y, w)
                match = check_close(*out_cy, *out_np)

                print(
                    f"{n:>10} {np.dtype(dtype).name:>8} {dup_fraction:>8} "
                    f"{t_cy * 1e3:>10.4f} {t_np * 1e3:>10.4f} "
                    f"{t_cy / t_np:>15.3f} {match!s:>6}"
                )

    print()
    print(
        "Exact-equivalent NumPy version (anchor-based grouping via "
        "searchsorted loop over groups, not over all samples)"
    )
    header2 = (
        f"{'n':>10} {'dtype':>8} {'dup_frac':>8} "
        f"{'cy (ms)':>10} {'np_exact (ms)':>14} {'speedup':>10} {'match':>6}"
    )
    print(header2)
    print("-" * len(header2))

    # Skip the largest, lowest-duplication cases: the group loop runs once
    # per group, so it degenerates to ~n Python-level iterations when
    # duplicates are rare.
    sizes_exact = [1_000, 100_000, 1_000_000]
    for dtype in dtypes:
        for n in sizes_exact:
            for dup_fraction in dup_fractions:
                X, y, w = make_data(n, dtype, dup_fraction, rng)

                number = 3
                repeat = 5

                t_cy = time_it(_make_unique_cy, (X, y, w), repeat=repeat, number=number)
                t_ex = time_it(
                    _make_unique_np_exact, (X, y, w), repeat=repeat, number=number
                )

                out_cy = _make_unique_cy(X, y, w)
                out_ex = _make_unique_np_exact(X, y, w)
                match = check_close(*out_cy, *out_ex)

                print(
                    f"{n:>10} {np.dtype(dtype).name:>8} {dup_fraction:>8} "
                    f"{t_cy * 1e3:>10.4f} {t_ex * 1e3:>14.4f} "
                    f"{t_cy / t_ex:>10.3f} {match!s:>6}"
                )
