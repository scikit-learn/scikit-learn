"""Compare sparse HGB with its dense-input baseline on synthetic wide data.

Run from the repository root, for example::

    python benchmarks/bench_hist_gradient_boosting_sparse.py

The script runs sparse and dense modes in separate processes so their peak
resident memory measurements do not include each other's allocations. The
dense mode materializes the same input with ``toarray`` and serves as the
memory and timing baseline; it is not a competing estimator.

The implementation keeps input and binned data in canonical CSR form. A
temporary CSC representation is used to estimate per-feature thresholds by
scanning the stored column entries. Sparse bin codes reserve an implicit
default code for the ordinary bin containing zero; histogram kernels accumulate
stored entries and derive that bin from the node totals. Split partitioning and
tree prediction use binary search within each CSR row. Dense arrays continue to
use their existing histogram and predictor kernels. A generic dense/sparse
accessor and full dense materialization were rejected because they either
complicate the dense hot path or restore sample-by-feature memory scaling.
Sparse categorical features remain unsupported.
"""

import argparse
import ctypes
import json
import os
import platform
import subprocess
import sys
import time

import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits

from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score


def _peak_rss_bytes():
    """Return peak resident set size using the current platform API."""
    if os.name == "nt":
        size_t = ctypes.c_size_t

        class ProcessMemoryCountersEx(ctypes.Structure):
            _fields_ = [
                ("cb", ctypes.c_ulong),
                ("PageFaultCount", ctypes.c_ulong),
                ("PeakWorkingSetSize", size_t),
                ("WorkingSetSize", size_t),
                ("QuotaPeakPagedPoolUsage", size_t),
                ("QuotaPagedPoolUsage", size_t),
                ("QuotaPeakNonPagedPoolUsage", size_t),
                ("QuotaNonPagedPoolUsage", size_t),
                ("PagefileUsage", size_t),
                ("PeakPagefileUsage", size_t),
                ("PrivateUsage", size_t),
            ]

        counters = ProcessMemoryCountersEx()
        counters.cb = ctypes.sizeof(counters)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        get_process_memory_info = psapi.GetProcessMemoryInfo
        get_process_memory_info.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ProcessMemoryCountersEx),
            ctypes.c_ulong,
        ]
        get_process_memory_info.restype = ctypes.c_int
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        process = kernel32.GetCurrentProcess()
        if not get_process_memory_info(process, ctypes.byref(counters), counters.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        return counters.PeakWorkingSetSize

    if sys.platform.startswith("linux"):
        with open("/proc/self/status") as status_file:
            for line in status_file:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) * 1024
        raise RuntimeError("Could not read VmHWM from /proc/self/status")

    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak * (1024 if sys.platform != "darwin" else 1))


def _make_data(n_samples, n_features, density, seed):
    rng = np.random.default_rng(seed)
    n_signal_features = min(8, n_features)
    signal = sparse.random(
        n_samples,
        n_signal_features,
        density=0.3,
        format="csr",
        random_state=rng,
        data_rvs=lambda size: rng.normal(size=size),
    )
    n_noise_features = n_features - n_signal_features
    if n_noise_features:
        noise = sparse.random(
            n_samples,
            n_noise_features,
            density=density,
            format="csr",
            random_state=rng,
            data_rvs=lambda size: rng.normal(size=size),
        )
    else:
        noise = sparse.csr_matrix((n_samples, 0))
    X = sparse.hstack((signal, noise), format="csr")
    weights = np.linspace(1.0, 0.25, n_signal_features)
    scores = np.asarray(signal @ weights).ravel()
    y = (scores >= np.median(scores)).astype(np.uint8)
    return X, y


def _run_worker(args):
    X, y = _make_data(args.n_samples, args.n_features, args.density, args.seed)
    split = int(args.n_samples * 0.8)
    X_train_sparse = X[:split].tocsr()
    X_test_sparse = X[split:].tocsr()
    y_train, y_test = y[:split], y[split:]
    sparse_bytes = X_train_sparse.data.nbytes + X_train_sparse.indices.nbytes
    sparse_bytes += X_train_sparse.indptr.nbytes

    if args.worker == "dense":
        X_train = X_train_sparse.toarray()
        X_test = X_test_sparse.toarray()
        del X, X_train_sparse, X_test_sparse
        input_format = "dense"
    else:
        X_train, X_test = X_train_sparse, X_test_sparse
        del X
        input_format = "CSR"

    model = HistGradientBoostingClassifier(
        max_iter=args.n_iterations,
        max_leaf_nodes=15,
        min_samples_leaf=10,
        max_bins=64,
        early_stopping=False,
        random_state=args.seed,
    )
    with threadpool_limits(limits=args.n_threads):
        start = time.perf_counter()
        model.fit(X_train, y_train)
        training_seconds = time.perf_counter() - start

        start = time.perf_counter()
        prediction = model.predict(X_test)
        prediction_seconds = time.perf_counter() - start

    result = {
        "mode": input_format,
        "n_samples": args.n_samples,
        "n_features": args.n_features,
        "density": args.density,
        "nnz": int(X_train_sparse.nnz) if args.worker == "sparse" else None,
        "sparse_input_mib": sparse_bytes / 2**20,
        "dense_input_mib": split
        * args.n_features
        * np.dtype(np.float64).itemsize
        / 2**20,
        "training_seconds": training_seconds,
        "prediction_seconds": prediction_seconds,
        "prediction_samples_per_second": len(y_test) / prediction_seconds,
        "test_accuracy": accuracy_score(y_test, prediction),
        "peak_rss_mib": _peak_rss_bytes() / 2**20,
        "python": platform.python_version(),
        "machine": platform.machine(),
        "n_threads": args.n_threads,
    }
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=6000)
    parser.add_argument("--n-features", type=int, default=3000)
    parser.add_argument("--density", type=float, default=0.005)
    parser.add_argument("--n-iterations", type=int, default=5)
    parser.add_argument("--n-threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--worker", choices=("sparse", "dense"))
    args = parser.parse_args()

    if args.worker:
        _run_worker(args)
        return

    results = []
    for worker in ("sparse", "dense"):
        command = [
            sys.executable,
            __file__,
            "--worker",
            worker,
            "--n-samples",
            str(args.n_samples),
            "--n-features",
            str(args.n_features),
            "--density",
            str(args.density),
            "--n-iterations",
            str(args.n_iterations),
            "--n-threads",
            str(args.n_threads),
            "--seed",
            str(args.seed),
        ]
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
        results.append(json.loads(completed.stdout.strip().splitlines()[-1]))

    print(
        f"Synthetic HGB benchmark: {args.n_samples:,} rows, "
        f"{args.n_features:,} columns, noise density={args.density:.3%}, "
        f"{args.n_iterations} trees"
    )
    print(
        "mode    train(s)  predict(s)  predict/s  peak RSS(MiB)  accuracy  input(MiB)"
    )
    for result in results:
        input_mib = (
            result["sparse_input_mib"]
            if result["mode"] == "CSR"
            else result["dense_input_mib"]
        )
        print(
            f"{result['mode']:<6}  {result['training_seconds']:>8.3f}  "
            f"{result['prediction_seconds']:>10.4f}  "
            f"{result['prediction_samples_per_second']:>9.0f}  "
            f"{result['peak_rss_mib']:>13.1f}  "
            f"{result['test_accuracy']:>8.3f}  {input_mib:>10.1f}"
        )
    print(
        f"Environment: Python {results[0]['python']}, {results[0]['machine']}, "
        f"{args.n_threads} thread(s)"
    )


if __name__ == "__main__":
    main()
