"""Instructions
1. Build this PR and run:

```bash
python bench_missing_extratrees.py bench ~/bench_results pr
```

2. On main run:

```bash
python bench_missing_extratrees.py bench ~/bench_results main
```

3. Plotting

```bash
python bench_missing_extratrees.py plot ~/bench_results pr main results_image.png
```
"""


from functools import partial
import argparse
from time import perf_counter
from statistics import mean, stdev
from itertools import product
import csv
from pathlib import Path


from sklearn.tree import ExtraTreeRegressor, ExtraTreeClassifier
from sklearn.datasets import make_classification, make_regression, make_low_rank_matrix
import numpy as np
from scipy.sparse import csc_matrix


def make_poisson_data(n_samples, n_features=50, random_state=0, has_missing=False):
    rng = np.random.RandomState(random_state)
    X = make_low_rank_matrix(
        n_samples=n_samples, n_features=n_features, random_state=rng
    )
    coef = rng.uniform(low=-2, high=2, size=n_features) / np.max(X, axis=0)
    y = rng.poisson(lam=np.exp(X @ coef))
    if has_missing:
        missing_mask = rng.choice([True, False], size=X.shape, p=[0.1, 0.9])
        X[missing_mask] = np.nan
    return X, y


def make_low_card_data(n_samples, n_features=50, random_state=0, has_missing=False):
    rng = np.random.RandomState(random_state)
    X = rng.choice([0.0, 1.0, 2.0], size=(n_samples, n_features))
    y = rng.choice([0, 1], size=n_samples)
    if has_missing:
        missing_mask = rng.choice([True, False], size=X.shape, p=[0.1, 0.9])
        X[missing_mask] = np.nan
    return X, y


def make_regression_custom(*args, has_missing=False, random_state=None, **kwargs):
    X, y = make_regression(*args, random_state=random_state, **kwargs)
    rng = np.random.RandomState()
    if has_missing:
        rng = np.random.RandomState(random_state)
        missing_mask = rng.choice([True, False], size=X.shape, p=[0.1, 0.9])
        X[missing_mask] = np.nan
    return X, y


def make_classification_custom(*args, has_missing=False, random_state=None, **kwargs):
    X, y = make_classification(*args, random_state=random_state, **kwargs)
    rng = np.random.RandomState()
    if has_missing:
        rng = np.random.RandomState(random_state)
        missing_mask = rng.choice([True, False], size=X.shape, p=[0.1, 0.9])
        X[missing_mask] = np.nan
    return X, y


N_REPEATS = 20

benchmark_config = [
    (
        ExtraTreeRegressor,
        list(
            product(
                ["squared_error"],
                [
                    partial(make_regression_custom, n_targets=2),
                    make_low_card_data,
                ],
                [40_000],
                ["dense"],
                ["random"],
                [False],
            )
        ),
    ),
    (
        ExtraTreeRegressor,
        list(
            product(
                ["poisson"],
                [make_poisson_data],
                [40_000],
                ["dense"],
                ["random"],
                [False],
            )
        ),
    ),
    (
        ExtraTreeClassifier,
        list(
            product(
                ["gini", "entropy"],
                [
                    partial(make_classification_custom, n_informative=10, n_classes=5),
                    make_low_card_data,
                ],
                [40_000],
                ["dense"],
                ["random"],
                [False],
            )
        ),
    ),
]


def bench(args):
    bench_results, branch = args.bench_results, args.branch
    results_dir = Path(bench_results)
    results_dir.mkdir(exist_ok=True)

    results_path = results_dir / f"{branch}.csv"

    with results_path.open("w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "criterion",
                "n_samples",
                "make_data",
                "container",
                "splitter",
                "has_missing",
                "n_repeat",
                "duration",
            ],
        )
        writer.writeheader()

        for Klass, items in benchmark_config:

            for config in items:
                (
                    criterion,
                    make_data,
                    n_samples,
                    container,
                    splitter,
                    has_missing,
                ) = config
                if isinstance(make_data, partial):
                    make_data_str = make_data.func.__name__
                else:
                    make_data_str = make_data.__name__

                default_config = {
                    "criterion": criterion,
                    "n_samples": n_samples,
                    "make_data": make_data_str,
                    "container": container,
                    "splitter": splitter,
                    "has_missing": has_missing,
                }
                combine_config = " ".join(f"{k}={v}" for k, v in default_config.items())

                klass_results = []
                for n_repeat in range(N_REPEATS):
                    X, y = make_data(
                        n_samples=n_samples,
                        random_state=n_repeat,
                        n_features=100,
                        has_missing=has_missing,
                    )
                    tree = Klass(
                        random_state=n_repeat, criterion=criterion, splitter=splitter
                    )

                    if container == "sparse":
                        X = csc_matrix(X, dtype=np.float32)

                    start = perf_counter()
                    tree.fit(X, y)
                    duration = perf_counter() - start
                    klass_results.append(duration)
                    writer.writerow(
                        {
                            **default_config,
                            **{
                                "n_repeat": n_repeat,
                                "duration": duration,
                            },
                        }
                    )
                results_mean, results_stdev = mean(klass_results), stdev(klass_results)
                print(
                    f"{combine_config} with {results_mean:.3f} +/- {results_stdev:.3f}"
                )


def plot(args):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    results_path = Path(args.bench_results)
    pr_path = results_path / f"{args.pr_name}.csv"
    main_path = results_path / f"{args.main_name}.csv"
    image_path = results_path / args.image_path

    df_pr = pd.read_csv(pr_path).assign(branch=args.pr_name)
    df_main = pd.read_csv(main_path).assign(branch=args.main_name)
    df_all = pd.concat((df_pr, df_main), ignore_index=True)

    df_all = df_all.assign(
        make_data=df_all["make_data"]
        .str.replace("_custom", "")
        .str.replace("make_", "")
        .str.replace("_data", "")
    )

    gb = df_all.groupby(["criterion", "make_data"])
    groups = gb.groups

    n_rows, n_cols = 2, 4
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 8), constrained_layout=True)
    axes_flat = axes.ravel()
    for i, (keys, idx) in enumerate(groups.items()):
        ax = axes_flat[i]
        ax.set_title(" | ".join(keys))
        sns.boxplot(data=df_all.loc[idx], y="duration", x="branch", ax=ax)
        if i % n_cols != 0:
            ax.set_ylabel("")

    axes_flat[-1].set_visible(False)

    fig.savefig(image_path)
    print(f"Saved image to {image_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    bench_parser = subparsers.add_parser("bench")
    bench_parser.add_argument("bench_results")
    bench_parser.add_argument("branch")
    bench_parser.set_defaults(func=bench)

    plot_parser = subparsers.add_parser("plot")
    plot_parser.add_argument("bench_results")
    plot_parser.add_argument("pr_name")
    plot_parser.add_argument("main_name")
    plot_parser.add_argument("image_path")
    plot_parser.set_defaults(func=plot)

    args = parser.parse_args()
    args.func(args)