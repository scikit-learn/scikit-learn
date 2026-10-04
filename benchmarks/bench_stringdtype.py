"""StringDType time/memory benchmark with isolated workers and correctness controls.

Example (requires psutil):
    python benchmarks/bench_stringdtype.py --rows 100000 --repeats 3 \
        --profiles short long_tail --output /tmp/stringdtype-bench.json

Run with each NumPy environment separately. For larger runs use --rows 1000000.
Every record runs in a fresh subprocess. Construction and estimator execution
are timed separately; reference comparisons run AFTER memory/time measurements.
RSS includes allocator overhead and retained temporary allocations, not just
live array storage. The 1 ms sampler can miss short peaks; process high-water
RSS also includes input generation. array.nbytes excludes out-of-buffer payload.
Object arrays deliberately share repeated Python string objects by default;
--object-layout fresh also tests independently allocated equal strings.
Missing data uses NaN, never the string 'nan'. Unicode+missing is excluded because
fixed-width Unicode cannot preserve that sentinel. Invalid label tasks are
excluded explicitly. There are no machine-dependent timing assertions.
"""

import argparse
import gc
import itertools
import json
import os
import platform
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np
import psutil
from scipy import sparse

import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import (
    LabelBinarizer,
    LabelEncoder,
    OneHotEncoder,
    OrdinalEncoder,
)
from sklearn.utils._encode import _encode, _unique
from sklearn.utils._missing import is_scalar_nan
from sklearn.utils._unique import _metadata_cache, attach_unique, cached_unique
from sklearn.utils.multiclass import check_classification_targets, type_of_target
from sklearn.utils.validation import check_array

OPERATIONS = (
    "storage",
    "unique",
    "discovery",
    "mapping",
    "conversion",
    "object_unique",
    "targets",
    "ordinal",
    "onehot",
    "labels",
    "cache",
    "pipeline",
)
PROFILES = ("short", "long_tail", "unicode", "unique", "missing")


def make_values(rows, profile, layout):
    count = rows if profile == "unique" else 32
    pool = [f"category_{i:08d}" for i in range(count)]
    if profile == "unicode":
        pool = ["café_東京_🙂_" + value for value in pool]
    elif profile == "long_tail":
        pool[-1] += "x" * 1000
    rng = np.random.default_rng(42)
    indices = rng.integers(0, count, size=rows)
    if profile == "unique":
        indices = rng.permutation(rows)
    elif profile == "long_tail":
        indices %= count - 1
        indices[::1000] = count - 1
    values = [pool[i] for i in indices]
    if layout == "fresh":
        values = [value.encode("utf-8").decode("utf-8") for value in values]
    if profile == "missing":
        values[::20] = [np.nan] * len(values[::20])
    return values


def execute(operation, X, timings=None):
    def measure(name, function):
        start = time.perf_counter()
        result = function()
        if timings is not None:
            timings[name] = timings.get(name, 0) + time.perf_counter() - start
        return result

    if operation == "discovery":
        return dict(categories=measure("discovery", lambda: _unique(X)))
    if operation == "conversion":
        return dict(values=measure("to_object", lambda: np.asarray(X, dtype=object)))
    if operation == "object_unique":
        values = measure("to_object", lambda: np.asarray(X, dtype=object))
        return dict(categories=measure("object_discovery", lambda: _unique(values)))
    if operation == "mapping":
        categories = measure("discovery", lambda: _unique(X))
        codes = measure("mapping", lambda: _encode(X, uniques=categories))
        return dict(categories=categories, codes=codes)
    if operation == "targets":
        measure("target_validation", lambda: check_classification_targets(X))
        return dict(valid=True)

    if operation == "storage":
        return {
            "values": check_array(
                X, dtype=None, ensure_2d=False, ensure_all_finite="allow-nan"
            )
        }
    if operation == "unique":
        categories, inverse, counts = _unique(
            X, return_inverse=True, return_counts=True
        )
        return dict(categories=categories, inverse=inverse, counts=counts)
    if operation in ("ordinal", "onehot"):
        estimator = (
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1)
            if operation == "ordinal"
            else OneHotEncoder(handle_unknown="ignore")
        )
        column = X.reshape(-1, 1)
        measure("fit", lambda: estimator.fit(column))
        encoded = measure("transform", lambda: estimator.transform(column))
        return dict(encoded=encoded, categories=estimator.categories_[0])
    if operation == "labels":
        encoder, binarizer = LabelEncoder(), LabelBinarizer(sparse_output=True)
        codes = measure("label_encoder", lambda: encoder.fit_transform(X))
        binary = measure("label_binarizer", lambda: binarizer.fit_transform(X))
        # Repeated public target detection/scoring complements the cache microbenchmark.
        scores = measure("scoring", lambda: [accuracy_score(X, X) for _ in range(5)])
        return dict(
            codes=codes,
            binary=binary,
            classes=encoder.classes_,
            reconstructed=encoder.inverse_transform(codes),
            scores=scores,
            target_type=type_of_target(X),
        )
    if operation == "cache":
        with _metadata_cache():
            cached = attach_unique(X)
            categories = [cached_unique(cached) for _ in range(10)]
        return dict(categories=categories)
    if operation == "pipeline":
        # The numeric target depends on string content, not representation or sorting.
        y = measure(
            "target_construction",
            lambda: np.array([sum(str(value).encode("utf-8")) % 2 for value in X]),
        )
        split = len(X) * 4 // 5
        pipeline = make_pipeline(
            OneHotEncoder(handle_unknown="ignore"), LogisticRegression(max_iter=100)
        )
        measure(
            "pipeline_fit", lambda: pipeline.fit(X[:split].reshape(-1, 1), y[:split])
        )
        return dict(
            prediction=measure(
                "predict", lambda: pipeline.predict(X[split:].reshape(-1, 1))
            ),
            probabilities=measure(
                "predict_proba",
                lambda: pipeline.predict_proba(X[split:].reshape(-1, 1)),
            ),
            categories=pipeline[0].categories_[0],
        )
    raise ValueError(operation)


def equivalent(actual, expected):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for key in actual:
            equivalent(actual[key], expected[key])
    elif isinstance(actual, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            equivalent(a, b)
    elif sparse.issparse(actual):
        assert actual.shape == expected.shape
        difference = actual - expected
        np.testing.assert_allclose(difference.data, 0, atol=1e-10)
    else:
        a, b = np.asarray(actual), np.asarray(expected)
        assert a.shape == b.shape
        if a.dtype.kind in "OUT" or b.dtype.kind in "OUT":
            canonical = lambda x: [
                None if is_scalar_nan(v) else str(v) for v in x.ravel()
            ]
            assert canonical(a) == canonical(b)
        else:
            np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-8, equal_nan=True)


def worker(args):
    if args.profile == "missing" and (
        args.representation == "unicode"
        or args.operation in ("labels", "cache", "targets")
    ):
        return dict(status="excluded", reason="No equivalent sentinel/label contract")
    process = psutil.Process()
    rss_baseline = process.memory_info().rss
    values = make_values(args.rows, args.profile, args.object_layout)
    dtype = {
        "object": object,
        "unicode": str,
        "stringdtype": np.dtypes.StringDType(na_object=np.nan)
        if args.profile == "missing"
        else np.dtypes.StringDType(),
    }[args.representation]
    start = time.perf_counter()
    X = np.array(values, dtype=dtype)
    construction_seconds = time.perf_counter() - start
    del values
    gc.collect()
    rss_input = process.memory_info().rss
    samples = [rss_input]
    stop = threading.Event()

    def sample():
        while not stop.wait(0.001):
            samples.append(process.memory_info().rss)

    monitor = threading.Thread(target=sample)
    monitor.start()
    try:
        start = time.perf_counter()
        stage_seconds = {}
        actual = execute(args.operation, X, stage_seconds)
        seconds = time.perf_counter() - start
    finally:
        samples.append(process.memory_info().rss)
        stop.set()
        monitor.join()
    result = dict(
        status="pass",
        stage_seconds=stage_seconds,
        construction_seconds=construction_seconds,
        seconds=seconds,
        rss_baseline=rss_baseline,
        rss_input=rss_input,
        rss_after_operation=process.memory_info().rss,
        sampled_operation_peak=max(samples),
        buffer_bytes=X.nbytes,
    )
    if sys.platform != "win32":
        import resource

        highwater = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result["process_peak_rss"] = highwater * (
            1 if sys.platform == "darwin" else 1024
        )
    # Reference allocation and comparisons must not contaminate measured peaks.
    reference = execute(args.operation, np.asarray(X, dtype=object))
    equivalent(actual, reference)
    if args.operation == "object_unique":
        equivalent(actual, execute("discovery", X))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=10000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--profiles", nargs="+", choices=PROFILES, default=list(PROFILES)
    )
    parser.add_argument(
        "--operations", nargs="+", choices=OPERATIONS, default=list(OPERATIONS)
    )
    parser.add_argument(
        "--object-layout", choices=("shared", "fresh"), default="shared"
    )
    parser.add_argument("--output", type=Path, default=Path("stringdtype-results.json"))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=PROFILES, help=argparse.SUPPRESS)
    parser.add_argument("--operation", choices=OPERATIONS, help=argparse.SUPPRESS)
    parser.add_argument(
        "--representation",
        choices=("object", "unicode", "stringdtype"),
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    if args.rows < 100 or args.repeats < 1:
        parser.error("Use at least 100 rows and one repeat")
    if args.worker:
        try:
            print(json.dumps(worker(args)))
        except Exception as exc:
            print(
                json.dumps(
                    dict(
                        status="failure",
                        error=type(exc).__name__,
                        message=str(exc),
                        traceback=traceback.format_exc(),
                    )
                )
            )
        return
    env = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        LOKY_MAX_CPU_COUNT="1",
    )
    git = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(sklearn.__file__).parent,
        capture_output=True,
        text=True,
    )
    dirty = (
        subprocess.run(
            ["git", "diff", "--quiet", "HEAD"], cwd=Path(sklearn.__file__).parent
        ).returncode
        != 0
    )
    report = dict(
        benchmark_schema=2,
        working_tree_dirty=dirty,
        numpy=np.__version__,
        sklearn=sklearn.__version__,
        sklearn_path=sklearn.__file__,
        commit=git.stdout.strip(),
        python=sys.version,
        platform=platform.platform(),
        rows=args.rows,
        object_layout=args.object_layout,
        results=[],
    )
    for repeat, profile, operation, representation in itertools.product(
        range(args.repeats),
        args.profiles,
        args.operations,
        ("object", "unicode", "stringdtype"),
    ):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            "--rows",
            str(args.rows),
            "--profile",
            profile,
            "--operation",
            operation,
            "--representation",
            representation,
            "--object-layout",
            args.object_layout,
        ]
        try:
            completed = subprocess.run(
                command,
                env=env,
                capture_output=True,
                text=True,
                timeout=300,
                check=True,
            )
            row = json.loads(completed.stdout)
            if completed.stderr:
                row["warnings"] = completed.stderr[-4000:]
        except (subprocess.SubprocessError, ValueError) as exc:
            row = dict(status="worker_failure", message=str(exc))
        row.update(
            repeat=repeat,
            profile=profile,
            operation=operation,
            representation=representation,
        )
        report["results"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(
            f"{repeat} {profile} {operation} {representation}: {row['status']}",
            flush=True,
        )
    if any(r["status"] in ("failure", "worker_failure") for r in report["results"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
