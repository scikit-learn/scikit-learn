from functools import partial

import numpy as np
import pytest

from sklearn.metrics import (
    log_loss,
    mean_absolute_error,
    mean_poisson_deviance,
    mean_squared_error,
)
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

CLF_CRITERIONS = ("gini", "log_loss")
REG_CRITERIONS = ("squared_error", "absolute_error", "friedman_mse", "poisson")


# TODO: use GRS;
@pytest.mark.parametrize("duplicates", ["x"])  # , "duplicates"])
@pytest.mark.parametrize("missing_values", ["x"])  # , "missing_values"])
@pytest.mark.parametrize(
    "criterion",
    ["gini", "log_loss", "squared_error", "absolute_error", "friedman_mse", "poisson"],
)
def test_best_split_optimality(duplicates, missing_values, criterion):
    is_clf = criterion in CLF_CRITERIONS
    with_duplicates = duplicates != "x"
    with_nans = missing_values != "x"
    n_classes = np.random.randint(2, 5)

    def gini_loss(y, y_pred, sample_weight):
        p = y_pred[0]
        return (p * (1 - p)).sum()

    losses = {
        "poisson": mean_poisson_deviance,
        "squared_error": mean_squared_error,
        "absolute_error": mean_absolute_error,
        "log_loss": partial(log_loss, labels=list(range(n_classes))),
        "gini": gini_loss,
    }

    def weighted_mean(y, w):
        return (y * w).sum() / w.sum()

    def weighted_median(y, w):
        sorter = np.argsort(y)
        wc = np.cumsum(w[sorter])
        idx = np.searchsorted(wc, wc[-1] / 2)
        return y[sorter[idx]]

    def class_ratios(y, w):
        return np.clip(np.bincount(y, w, minlength=n_classes) / w.sum(), None, 1)

    predictors = {
        "poisson": weighted_mean,
        "squared_error": weighted_mean,
        "absolute_error": weighted_median,
        "log_loss": class_ratios,
        "gini": class_ratios,
    }

    def compute_child_loss(y: np.ndarray, w: np.ndarray):
        if y.size == 0:
            return np.inf
        pred_dim = (y.size, n_classes) if is_clf else (y.size,)
        y_pred = np.empty(pred_dim)
        y_pred[:] = predictors[criterion](y, w)
        return w.sum() * losses[criterion](y, y_pred, sample_weight=w)

    def compute_split_loss(x, y, w, threshold, missing_left=False):
        mask = x < threshold
        if missing_left:
            mask |= np.isnan(x)
        if criterion == "friedman_mse":
            diff = weighted_mean(y[mask], w[mask]) - weighted_mean(y[~mask], w[~mask])
            return (-(diff**2) * w[mask].sum() * w[~mask].sum() / w.sum(),)
        return (
            compute_child_loss(y[mask], w[mask]),
            compute_child_loss(y[~mask], w[~mask]),
        )

    def compute_all_losses(x, y, w, missing_left=False):
        nan_mask = np.isnan(x)
        xu = np.unique(x[~nan_mask], sorted=True)
        thresholds = (xu[1:] + xu[:-1]) / 2
        if nan_mask.any() and not missing_left:
            thresholds = np.append(thresholds, xu.max() * 2)
        return thresholds, [
            sum(compute_split_loss(x, y, w, threshold, missing_left))
            for threshold in thresholds
        ]

    def best_split_naive(X, y, w):
        splits = []
        for f in range(X.shape[1]):
            thresholds, losses = compute_all_losses(X[:, f], y, w)
            if with_nans:
                thresholds_, losses_ = compute_all_losses(
                    X[:, f], y, w, missing_left=True
                )
                thresholds = np.concat((thresholds, thresholds_))
                losses = np.concat((losses, losses_))
            if len(losses) == 0:
                continue
            idx = np.argmin(losses)
            splits.append(
                (
                    losses[idx],
                    thresholds[idx],
                    with_nans and idx >= thresholds.size // 2,
                    f,
                )
            )
        return min(splits)

    d = 2
    for it, n in enumerate([5] * 10 + [10] * 10 + [30] * 5 + [100, 100, 200]):
        X = np.random.rand(n, d)
        y = np.random.rand(n) + X.sum(axis=1)
        w = np.random.rand(n)
        if is_clf:
            q = np.linspace(0, 1, num=n_classes + 1)[1:-1]
            y = np.searchsorted(np.quantile(y, q), y)
        if with_duplicates:
            X = X.round(1 if n < 50 else 2)
        if with_nans:
            for i in range(d):
                step = np.random.randint(2, 10)
                X[i::step, i] = np.nan

        Tree = DecisionTreeClassifier if is_clf else DecisionTreeRegressor
        tree = Tree(criterion=criterion, max_depth=1, max_features=d)
        tree.fit(X, y, sample_weight=w)
        best_split = best_split_naive(X, y, w)
        tree_loss = compute_split_loss(
            X[:, tree.tree_.feature[0]],
            y,
            w,
            tree.tree_.threshold[0],
            bool(tree.tree_.missing_go_to_left[0]),
        )
        tree_split = (
            sum(tree_loss),
            tree.tree_.threshold[0],
            bool(tree.tree_.missing_go_to_left[0]),
            tree.tree_.feature[0],
        )
        assert np.isclose(best_split[0], tree_split[0]), (it, best_split, tree_split)

        vals = tree.tree_.impurity * tree.tree_.weighted_n_node_samples
        if criterion == "log_loss":
            vals *= np.log(2)
        if criterion == "poisson":
            vals *= 2
        if criterion != "friedman_mse":
            assert np.allclose(vals[1:], tree_loss), it
