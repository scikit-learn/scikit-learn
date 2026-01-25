"""
Manual grid search with warm_start for tree ensembles.

GridSearchCV fits each parameter setting from scratch and does not reuse a
previously fitted estimator. For estimators that support ``warm_start``, it is
possible to reuse the already fitted state when only increasing a "growable"
parameter such as ``n_estimators`` in tree ensembles.

This example shows a manual grid search over ``max_features`` and
``n_estimators`` for a RandomForestClassifier. For each value of
``max_features``, we progressively increase ``n_estimators`` with
``warm_start=True`` and record a validation score. We also compare the runtime
with a baseline GridSearchCV that evaluates the same parameter grid on the same
train/validation split.

Note: GridSearchCV is the recommended tool for most workflows because it
supports cross-validation, pipelines, and proper model selection practices.
This example is meant to illustrate a manual approach when warm-start can be
leveraged for a growing parameter.
"""

from __future__ import annotations

import time

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, PredefinedSplit, train_test_split


def manual_grid_search_warm_start(
    X_train,
    y_train,
    X_valid,
    y_valid,
    max_features_grid,
    n_estimators_grid,
    random_state=0
):
    """Manual grid search reusing fitted trees when increasing n_estimators."""
    n_estimators_grid = sorted(n_estimators_grid)
    scores = np.zeros((len(max_features_grid), len(n_estimators_grid)))

    start = time.perf_counter()

    for i, max_features in enumerate(max_features_grid):
        # New model per max_features value (cannot reuse across max_features).
        clf = RandomForestClassifier(
            warm_start=True,
            random_state=random_state,
            n_jobs=-1,
            max_features=max_features
        )

        # Reuse the fitted forest while increasing n_estimators.
        for j, n_estimators in enumerate(n_estimators_grid):
            clf.set_params(n_estimators=n_estimators)
            clf.fit(X_train, y_train)
            scores[i, j] = clf.score(X_valid, y_valid)

    elapsed = time.perf_counter() - start
    return scores, n_estimators_grid, elapsed


def gridsearchcv_baseline(
    X_train,
    y_train,
    X_valid,
    y_valid,
    max_features_grid,
    n_estimators_grid,
    random_state=0
):
    """Baseline GridSearchCV on the same train/validation split (fits from scratch)."""
    X_all = np.vstack([X_train, X_valid])
    y_all = np.hstack([y_train, y_valid])

    # PredefinedSplit: training rows are -1, validation rows are 0.
    test_fold = np.concatenate(
        [-np.ones(X_train.shape[0], dtype=int), np.zeros(X_valid.shape[0], dtype=int)]
    )
    cv = PredefinedSplit(test_fold=test_fold)

    param_grid = {
        "max_features": list(max_features_grid),
        "n_estimators": list(n_estimators_grid),
    }

    clf = RandomForestClassifier(
        random_state=random_state,
        n_jobs=-1,
    )

    gs = GridSearchCV(
        estimator=clf,
        param_grid=param_grid,
        cv=cv,
        scoring="accuracy",
        refit=False,
    )

    start = time.perf_counter()
    gs.fit(X_all, y_all)
    elapsed = time.perf_counter() - start

    best_index = int(np.argmax(gs.cv_results_["mean_test_score"]))
    best_score = float(gs.cv_results_["mean_test_score"][best_index])

    best_params = {
        "max_features": gs.cv_results_["param_max_features"][best_index],
        "n_estimators": int(gs.cv_results_["param_n_estimators"][best_index]),
    }

    return best_params, best_score, elapsed


def best_from_score_surface(scores, max_features_grid, n_estimators_grid):
    """Return best (max_features, n_estimators) from the score matrix."""
    best_flat = int(np.argmax(scores))
    i, j = np.unravel_index(best_flat, scores.shape)
    return (
        max_features_grid[i],
        n_estimators_grid[j],
        float(scores[i, j]),
    )


def plot_heatmap(scores, max_features_grid, n_estimators_grid):
    fig, ax = plt.subplots()
    im = ax.imshow(scores, aspect="auto", origin="lower")

    ax.set_xticks(np.arange(len(n_estimators_grid)))
    ax.set_xticklabels(n_estimators_grid)
    ax.set_xlabel("n_estimators")

    ax.set_yticks(np.arange(len(max_features_grid)))
    ax.set_yticklabels([str(v) for v in max_features_grid])
    ax.set_ylabel("max_features")

    ax.set_title("Validation score for each (max_features, n_estimators) pair")
    fig.colorbar(im, ax=ax, label="accuracy")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    X, y = make_classification(
        n_samples=2500,
        n_features=30,
        n_informative=10,
        n_redundant=10,
        random_state=0
    )

    # Split into train+valid and test, then split train+valid into train and valid.
    X_train_valid, X_test, y_train_valid, y_test = train_test_split(
        X, y, test_size=0.25, random_state=0, stratify=y
    )
    X_train, X_valid, y_train, y_valid = train_test_split(
        X_train_valid,
        y_train_valid,
        test_size=0.25,
        random_state=0,
        stratify=y_train_valid,
    )

    max_features_grid = ["sqrt", 0.5, 1.0]
    n_estimators_grid = [25, 50, 100, 150]

    scores, n_estimators_sorted, t_warm = manual_grid_search_warm_start(
        X_train,
        y_train,
        X_valid,
        y_valid,
        max_features_grid,
        n_estimators_grid,
        random_state=0,
    )

    best_mf, best_ne, best_score = best_from_score_surface(
        scores, max_features_grid, n_estimators_sorted
    )

    gs_best_params, gs_best_score, t_gs = gridsearchcv_baseline(
        X_train,
        y_train,
        X_valid,
        y_valid,
        max_features_grid,
        n_estimators_grid,
        random_state=0,
    )

    print("Manual warm_start grid search")
    print(f"  best params: max_features={best_mf}, n_estimators={best_ne}")
    print(f"  best valid accuracy: {best_score:.4f}")
    print(f"  runtime: {t_warm:.3f} s")
    print()

    print("GridSearchCV baseline (same split, fits from scratch)")
    print(
        "  best params: max_features={mf}, n_estimators={ne}".format(
            mf=gs_best_params["max_features"], ne=gs_best_params["n_estimators"]
        )
    )
    print(f"  best valid accuracy: {gs_best_score:.4f}")
    print(f"  runtime: {t_gs:.3f} s")
    print()

    # Optional: fit the best manual config on train+valid and evaluate on test.
    final_clf = RandomForestClassifier(
        random_state=0,
        n_jobs=-1,
        max_features=best_mf,
        n_estimators=best_ne,
    )
    final_clf.fit(X_train_valid, y_train_valid)
    test_score = final_clf.score(X_test, y_test)
    print(f"Test accuracy using manual best params: {test_score:.4f}")

    plot_heatmap(scores, max_features_grid, n_estimators_sorted)
