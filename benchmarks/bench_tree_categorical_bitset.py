"""
Benchmark categorical decision-tree bitset performance.

This script is intended to be run twice: once before and once after a bitset
storage/helper refactor. Example:

    python benchmarks/bench_tree_categorical_bitset.py \
        --label old --output ./tree-bitset-old.json
    python benchmarks/bench_tree_categorical_bitset.py \
        --label new --output ./tree-bitset-new.json
    python benchmarks/bench_tree_categorical_bitset.py \
        --compare ./tree-bitset-old.json ./tree-bitset-new.json
"""

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path
from statistics import median

import numpy as np


def _git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except Exception:
        return None


def make_categorical_regression(n_samples, n_features, n_categories, random_state):
    rng = np.random.default_rng(random_state)
    X = rng.integers(
        0, n_categories, size=(n_samples, n_features), dtype=np.int32
    ).astype(np.float32)

    effects = rng.normal(size=(n_features, n_categories))
    effects -= effects.mean(axis=1, keepdims=True)
    y = np.zeros(n_samples, dtype=np.float64)
    for feature_idx in range(n_features):
        y += effects[feature_idx, X[:, feature_idx].astype(np.intp)]
    y += 0.01 * rng.normal(size=n_samples)
    return X, y


def make_categorical_classification(n_samples, n_features, n_categories, random_state):
    X, score = make_categorical_regression(
        n_samples, n_features, n_categories, random_state
    )
    y = score > np.median(score)
    return X, y.astype(np.intp)


def time_call(func, *, warmups, repeats):
    for _ in range(warmups):
        func()

    timings = []
    for _ in range(repeats):
        start = time.perf_counter()
        func()
        timings.append(time.perf_counter() - start)

    return {
        "best": min(timings),
        "median": median(timings),
        "all": timings,
    }


def count_categorical_splits(estimator):
    tree = estimator.tree_
    feature = tree.feature
    internal = feature >= 0
    n_categories = estimator.n_categories_in_feature_
    return int(np.sum(internal & (n_categories[feature] > 0)))


def tree_info(estimator):
    return {
        "depth": int(estimator.tree_.max_depth),
        "node_count": int(estimator.tree_.node_count),
        "categorical_splits": count_categorical_splits(estimator),
    }


def benchmark_one_estimator(
    estimator,
    X_train,
    y_train,
    X_pred,
    *,
    categorical_features,
    warmups,
    repeats,
    grid_resolution,
    include_partial_dependence,
):
    fit_params = {
        "max_depth": estimator.max_depth,
        "min_samples_leaf": estimator.min_samples_leaf,
        "random_state": estimator.random_state,
        "categorical_features": categorical_features,
    }

    def fit_once():
        return estimator.__class__(**fit_params).fit(X_train, y_train)

    fitted = fit_once()

    target_feature = [0]
    custom_values = {
        0: np.arange(
            min(grid_resolution, int(np.max(X_train[:, 0])) + 1), dtype=np.float32
        )
    }

    timings = {
        "fit": time_call(fit_once, warmups=warmups, repeats=repeats),
        "predict": time_call(
            lambda: fitted.predict(X_pred), warmups=warmups, repeats=repeats
        ),
        "decision_path": time_call(
            lambda: fitted.decision_path(X_pred), warmups=warmups, repeats=repeats
        ),
    }

    if include_partial_dependence:
        from sklearn.inspection import partial_dependence

        timings["partial_dependence_recursion"] = time_call(
            lambda: partial_dependence(
                fitted,
                X_train,
                target_feature,
                method="recursion",
                custom_values=custom_values,
            ),
            warmups=warmups,
            repeats=repeats,
        )

    return {
        "timings": timings,
        "tree": tree_info(fitted),
    }


def run_benchmark(args):
    from sklearn import __version__ as sklearn_version
    from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

    categorical_features = np.arange(args.n_features)
    results = {
        "label": args.label,
        "git_commit": _git_commit(),
        "python": sys.version,
        "platform": platform.platform(),
        "sklearn_version": sklearn_version,
        "params": {
            "n_samples": args.n_samples,
            "n_pred_samples": args.n_pred_samples,
            "n_features": args.n_features,
            "categories": args.categories,
            "max_depth": args.max_depth,
            "min_samples_leaf": args.min_samples_leaf,
            "warmups": args.warmups,
            "repeats": args.repeats,
            "grid_resolution": args.grid_resolution,
        },
        "cases": [],
    }

    for n_categories in args.categories:
        X_reg, y_reg = make_categorical_regression(
            args.n_samples, args.n_features, n_categories, args.random_state
        )
        X_clf, y_clf = make_categorical_classification(
            args.n_samples, args.n_features, n_categories, args.random_state
        )
        X_pred = X_reg[: args.n_pred_samples]

        print(f"Benchmarking n_categories={n_categories}")
        regressor = DecisionTreeRegressor(
            max_depth=args.max_depth,
            min_samples_leaf=args.min_samples_leaf,
            random_state=args.random_state,
            categorical_features=categorical_features,
        )
        classifier = DecisionTreeClassifier(
            max_depth=args.max_depth,
            min_samples_leaf=args.min_samples_leaf,
            random_state=args.random_state,
            categorical_features=categorical_features,
        )

        case = {
            "n_categories": n_categories,
            "regressor": benchmark_one_estimator(
                regressor,
                X_reg,
                y_reg,
                X_pred,
                categorical_features=categorical_features,
                warmups=args.warmups,
                repeats=args.repeats,
                grid_resolution=args.grid_resolution,
                include_partial_dependence=True,
            ),
            "classifier": benchmark_one_estimator(
                classifier,
                X_clf,
                y_clf,
                X_pred,
                categorical_features=categorical_features,
                warmups=args.warmups,
                repeats=args.repeats,
                grid_resolution=args.grid_resolution,
                include_partial_dependence=False,
            ),
        }
        results["cases"].append(case)
        print_case_summary(case)

    if args.output is not None:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(results, indent=2, sort_keys=True))
        print(f"Wrote {output}")

    return results


def format_seconds(seconds):
    return f"{seconds:.6f}s"


def print_case_summary(case):
    print(f"n_categories={case['n_categories']}")
    for estimator_name in ("regressor", "classifier"):
        tree = case[estimator_name]["tree"]
        print(
            f"  {estimator_name}: depth={tree['depth']} "
            f"nodes={tree['node_count']} "
            f"categorical_splits={tree['categorical_splits']}"
        )
        for workload, timing in case[estimator_name]["timings"].items():
            print(
                f"    {workload}: median={format_seconds(timing['median'])} "
                f"best={format_seconds(timing['best'])}"
            )


def compare_results(args):
    baseline = json.loads(Path(args.compare[0]).read_text())
    candidate = json.loads(Path(args.compare[1]).read_text())

    print(f"baseline: {baseline.get('label')} {baseline.get('git_commit')}")
    print(f"candidate: {candidate.get('label')} {candidate.get('git_commit')}")
    print()

    baseline_cases = {case["n_categories"]: case for case in baseline["cases"]}
    candidate_cases = {case["n_categories"]: case for case in candidate["cases"]}

    for n_categories in sorted(baseline_cases):
        if n_categories not in candidate_cases:
            continue
        print(f"n_categories={n_categories}")
        for estimator_name in ("regressor", "classifier"):
            print(f"  {estimator_name}")
            base_timings = baseline_cases[n_categories][estimator_name]["timings"]
            cand_timings = candidate_cases[n_categories][estimator_name]["timings"]
            for workload in base_timings:
                base = base_timings[workload]["median"]
                cand = cand_timings[workload]["median"]
                pct = 100 * (cand - base) / base
                print(
                    f"    {workload}: {format_seconds(base)} -> "
                    f"{format_seconds(cand)} ({pct:+.2f}%)"
                )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", default="run", help="Label stored in JSON output.")
    parser.add_argument("--output", help="Path to write JSON benchmark results.")
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("BASELINE_JSON", "CANDIDATE_JSON"),
        help="Compare two JSON outputs instead of running benchmarks.",
    )
    parser.add_argument("--n-samples", type=int, default=50_000)
    parser.add_argument("--n-pred-samples", type=int, default=50_000)
    parser.add_argument("--n-features", type=int, default=12)
    parser.add_argument(
        "--categories",
        type=int,
        nargs="+",
        default=[64, 128, 256],
        help="Category cardinalities to benchmark.",
    )
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--min-samples-leaf", type=int, default=20)
    parser.add_argument("--grid-resolution", type=int, default=256)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--random-state", type=int, default=0)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.compare is not None:
        compare_results(args)
    else:
        run_benchmark(args)


if __name__ == "__main__":
    main()
