#!/usr/bin/env python
"""Numpy-only reproducer for suspected BLIS DGEMM issues on macOS arm64."""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from typing import Any

import numpy as np


def _pairwise_via_gemm(x: np.ndarray) -> np.ndarray:
    x2 = np.einsum("ij,ij->i", x, x)
    d2 = x2[:, None] - 2.0 * (x @ x.T) + x2[None, :]
    np.maximum(d2, 0.0, out=d2)
    return np.sqrt(d2, out=d2)


def _pairwise_reference(x: np.ndarray) -> np.ndarray:
    diff = x[:, None, :] - x[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2))


def _digest(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).view(np.uint8)).hexdigest()


def _detected_blas_name() -> str:
    if hasattr(np.__config__, "show"):
        try:
            cfg: dict[str, Any] = np.__config__.show(mode="dicts")
            return str(cfg.get("Build Dependencies", {}).get("blas", {}).get("name", "")).lower()
        except TypeError:
            # Older NumPy versions do not support show(mode="dicts").
            pass
    return ""


def run(seed: int, n_samples: int, n_features: int, repeats: int, atol: float, rtol: float) -> int:
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((n_samples, n_features), dtype=np.float64)
    # Use a non-trivial memory layout to stress GEMM code paths.
    x = np.asfortranarray(x[:, ::-1])

    ref = _pairwise_reference(x)

    failures = []
    digests = []
    for i in range(repeats):
        d = _pairwise_via_gemm(x)
        digests.append(_digest(d))

        if not np.isfinite(d).all():
            failures.append(f"repeat={i}: non-finite values detected")
            continue

        if not np.allclose(d, ref, atol=atol, rtol=rtol):
            abs_err = np.max(np.abs(d - ref))
            rel_err = np.max(np.abs(d - ref) / np.maximum(np.abs(ref), 1e-15))
            failures.append(
                f"repeat={i}: mismatch (max_abs_err={abs_err:.6g}, max_rel_err={rel_err:.6g})"
            )

    print(f"numpy={np.__version__}")
    print(f"detected_blas={_detected_blas_name() or 'unknown'}")
    print(f"BLIS_NUM_THREADS={os.getenv('BLIS_NUM_THREADS')}")
    print(f"OPENBLAS_NUM_THREADS={os.getenv('OPENBLAS_NUM_THREADS')}")
    print(f"ACCELERATE_LOG_LEVEL={os.getenv('ACCELERATE_LOG_LEVEL')}")
    print(f"unique_result_digests={len(set(digests))}/{len(digests)}")

    if failures:
        print("FAIL")
        for failure in failures[:10]:
            print(f"  - {failure}")
        return 1

    if len(set(digests)) != 1:
        print("FAIL")
        print("  - output changed across repeats")
        return 1

    print("PASS")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-samples", type=int, default=384)
    parser.add_argument("--n-features", type=int, default=192)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--expected-blas", choices=["blis", "openblas", "newaccelerate"], default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    detected_blas = _detected_blas_name()
    if args.expected_blas is not None:
        expected_fragment = "accelerate" if args.expected_blas == "newaccelerate" else args.expected_blas
        if expected_fragment not in detected_blas:
            print(
                f"FAIL: requested BLAS '{args.expected_blas}' but NumPy reports '{detected_blas or 'unknown'}'"
            )
            sys.exit(1)
    sys.exit(
        run(
            seed=args.seed,
            n_samples=args.n_samples,
            n_features=args.n_features,
            repeats=args.repeats,
            atol=args.atol,
            rtol=args.rtol,
        )
    )
