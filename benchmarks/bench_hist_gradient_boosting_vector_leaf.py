"""Benchmark trees with vector leaves in HistGradientBoosting estimators.

Compare ``vector_leaf=True`` (a single tree per boosting iteration with one
leaf value per output) with ``vector_leaf=False`` (one tree per boosting
iteration and per output) for multiclass classification, multilabel
classification and multi-quantile regression, on a grid of hyper-parameters.

Each configuration is scored on the training set and on a held-out test set
of ``TEST_SIZE`` samples (at least 1000). Lower is better:

- multiclass: log loss;
- multilabel: binary log loss averaged over the labels;
- quantile: pinball loss averaged over the quantile levels.

Results are appended to a CSV file (already computed configurations are
skipped). Use ``--plot`` to draw the Pareto plots of the fit time vs the
training scores and vs the held-out test scores.

Example::

    python bench_hist_gradient_boosting_vector_leaf.py --output results.csv
    python bench_hist_gradient_boosting_vector_leaf.py --output results.csv --plot
"""

import argparse
import csv
import itertools
import os
from time import perf_counter

import numpy as np

from sklearn.datasets import (
    fetch_california_housing,
    fetch_covtype,
    fetch_openml,
    make_classification,
    make_multilabel_classification,
)
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from sklearn.metrics import log_loss, mean_pinball_loss
from sklearn.model_selection import train_test_split

TEST_SIZE = 1000

GRID = {
    "max_iter": [30, 100, 300],
    "max_leaf_nodes": [7, 31, 127],
    "learning_rate": [0.03, 0.1, 0.3],
    "l2_regularization": [0.0, 1.0, 10.0],
}
KEY_FIELDS = ["dataset", "task", "n_outputs", "vector_leaf", *GRID]
FIELDS = [*KEY_FIELDS, "fit_time", "train_loss", "test_loss"]


def covtype():
    X, y = fetch_covtype(return_X_y=True)
    split = train_test_split(
        X, y, train_size=30_000, test_size=TEST_SIZE, stratify=y, random_state=0
    )
    return "multiclass", *split, {}


def synthetic_multiclass():
    X, y = make_classification(
        n_samples=30_000,
        n_features=30,
        n_informative=15,
        n_classes=20,
        n_clusters_per_class=1,
        random_state=0,
    )
    split = train_test_split(X, y, test_size=TEST_SIZE, stratify=y, random_state=0)
    return "multiclass", *split, {}


def yeast():
    X, y = fetch_openml("yeast", version=4, return_X_y=True, as_frame=False)
    y = y.astype(np.float64)
    split = train_test_split(X, y, test_size=TEST_SIZE, random_state=0)
    return "multilabel", *split, {}


def synthetic_multilabel():
    X, y = make_multilabel_classification(
        n_samples=30_000, n_features=40, n_classes=20, n_labels=4, random_state=0
    )
    split = train_test_split(X, y, test_size=TEST_SIZE, random_state=0)
    return "multilabel", *split, {}


def california():
    X, y = fetch_california_housing(return_X_y=True)
    quantiles = np.linspace(0.1, 0.9, 9).round(2).tolist()
    split = train_test_split(X, y, test_size=TEST_SIZE, random_state=0)
    return "quantile", *split, {"loss": "quantile", "quantile": quantiles}


def synthetic_quantile():
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(30_000, 10))
    mean = 10 * np.sin(np.pi * X[:, 0] * X[:, 1]) + 20 * (X[:, 2] - 0.5) ** 2
    scale = 0.5 + 3 * X[:, 3] + 2 * X[:, 4] ** 2
    y = mean + scale * rng.standard_t(df=3, size=X.shape[0])
    quantiles = np.linspace(0.05, 0.95, 19).round(2).tolist()
    split = train_test_split(X, y, test_size=TEST_SIZE, random_state=0)
    return "quantile", *split, {"loss": "quantile", "quantile": quantiles}


DATASETS = {
    "covtype": covtype,
    "synthetic_multiclass": synthetic_multiclass,
    "yeast": yeast,
    "synthetic_multilabel": synthetic_multilabel,
    "california": california,
    "synthetic_quantile": synthetic_quantile,
}


def score(task, estimator, X, y, params):
    if task == "multiclass":
        proba = estimator.predict_proba(X)
        return log_loss(y, proba, labels=estimator.classes_)
    if task == "multilabel":
        proba = np.clip(estimator.predict_proba(X), 1e-15, 1 - 1e-15)
        return -np.mean(y * np.log(proba) + (1 - y) * np.log1p(-proba))
    pred = estimator.predict(X)
    return np.mean(
        [
            mean_pinball_loss(y, pred[:, j], alpha=q)
            for j, q in enumerate(params["quantile"])
        ]
    )


def run(datasets, output):
    done = set()
    if os.path.exists(output):
        with open(output) as f:
            for row in csv.DictReader(f):
                done.add(tuple(row[k] for k in KEY_FIELDS))
    write_header = not os.path.exists(output)
    with open(output, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()
        for name in datasets:
            task, X_train, X_test, y_train, y_test, params = DATASETS[name]()
            if task == "quantile":
                Estimator = HistGradientBoostingRegressor
                n_outputs = len(params["quantile"])
            else:
                Estimator = HistGradientBoostingClassifier
                if task == "multiclass":
                    n_outputs = len(np.unique(y_train))
                else:
                    n_outputs = y_train.shape[1]
            for values in itertools.product(*GRID.values()):
                grid_params = dict(zip(GRID, values))
                for vector_leaf in [False, True]:
                    row = dict(
                        dataset=name,
                        task=task,
                        n_outputs=n_outputs,
                        vector_leaf=vector_leaf,
                        **grid_params,
                    )
                    key = tuple(str(row[k]) for k in KEY_FIELDS)
                    if key in done:
                        continue
                    estimator = Estimator(
                        vector_leaf=vector_leaf,
                        early_stopping=False,
                        random_state=0,
                        **grid_params,
                        **params,
                    )
                    tic = perf_counter()
                    estimator.fit(X_train, y_train)
                    row["fit_time"] = perf_counter() - tic
                    row["train_loss"] = float(
                        score(task, estimator, X_train, y_train, params)
                    )
                    row["test_loss"] = float(
                        score(task, estimator, X_test, y_test, params)
                    )
                    writer.writerow(row)
                    f.flush()
                    print(row, flush=True)


def pareto_front(fit_time, loss):
    order = np.lexsort((loss, fit_time))
    front, best = [], np.inf
    for i in order:
        if loss[i] < best:
            front.append(i)
            best = loss[i]
    return np.asarray(front)


def plot(output, figure, loss_column, split_label):
    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.read_csv(output)
    datasets = list(dict.fromkeys(df["dataset"]))
    n_cols = 3
    n_rows = int(np.ceil(len(datasets) / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(6 * n_cols, 4.5 * n_rows), squeeze=False
    )
    task_loss = {
        "multiclass": "log loss",
        "multilabel": "binary log loss (mean over labels)",
        "quantile": "pinball loss (mean over quantiles)",
    }
    for ax, name in zip(axes.ravel(), datasets):
        sub = df[df["dataset"] == name]
        task = sub["task"].iloc[0]
        for vector_leaf, color in [(False, "tab:blue"), (True, "tab:orange")]:
            s = sub[sub["vector_leaf"] == vector_leaf]
            fit_time = s["fit_time"].to_numpy()
            loss = s[loss_column].to_numpy()
            ax.scatter(fit_time, loss, s=10, alpha=0.3, color=color)
            front = pareto_front(fit_time, loss)
            ax.step(
                fit_time[front],
                loss[front],
                where="post",
                color=color,
                marker="o",
                markersize=4,
                label=f"vector_leaf={vector_leaf}",
            )
        best = sub[loss_column].min()
        ax.set_ylim(best - 0.02 * abs(best), np.quantile(sub[loss_column], 0.6))
        ax.set_xscale("log")
        ax.set_xlabel("fit time (s)")
        ax.set_ylabel(f"{split_label} {task_loss[task]}")
        n_outputs = sub["n_outputs"].iloc[0]
        ax.set_title(f"{name} ({task}, {n_outputs} outputs)")
        ax.legend()
    for ax in axes.ravel()[len(datasets) :]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(figure, dpi=120)
    print(f"Saved {figure}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS))
    parser.add_argument("--output", default="hgb_vector_leaf.csv")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument(
        "--figure-prefix",
        default="hgb_vector_leaf_pareto",
        help="Plotting writes {prefix}_train.png and {prefix}_test.png.",
    )
    args = parser.parse_args()
    if args.plot:
        plot(
            args.output,
            f"{args.figure_prefix}_train.png",
            "train_loss",
            "training set",
        )
        plot(
            args.output,
            f"{args.figure_prefix}_test.png",
            "test_loss",
            f"held-out test (n={TEST_SIZE})",
        )
    else:
        run(args.datasets, args.output)
