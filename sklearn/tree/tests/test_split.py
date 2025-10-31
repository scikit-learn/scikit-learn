from dataclasses import dataclass
from functools import partial

import numpy as np
import pytest
from scipy.sparse import csc_array

from sklearn.metrics import (
    log_loss,
    mean_absolute_error,
    mean_poisson_deviance,
    mean_squared_error,
)
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

CLF_CRITERIONS = ("gini", "log_loss")
REG_CRITERIONS = ("squared_error", "absolute_error", "friedman_mse", "poisson")


def gini_loss(y, y_pred, sample_weight):
    p = y_pred[0]
    return (p * (1 - p)).sum()


@dataclass
class NaiveSplitter:
    is_clf: bool
    criterion: str
    with_nans: bool
    n_classes: int

    @staticmethod
    def weighted_mean(y, w):
        return (y * w).sum() / w.sum()

    @staticmethod
    def weighted_median(y, w):
        sorter = np.argsort(y)
        wc = np.cumsum(w[sorter])
        idx = np.searchsorted(wc, wc[-1] / 2)
        return y[sorter[idx]]

    def class_ratios(self, y, w):
        return np.clip(np.bincount(y, w, minlength=self.n_classes) / w.sum(), None, 1)

    @property
    def loss(self):
        losses = {
            "poisson": mean_poisson_deviance,
            "squared_error": mean_squared_error,
            "absolute_error": mean_absolute_error,
            "log_loss": partial(log_loss, labels=list(range(self.n_classes))),
            "gini": gini_loss,
        }
        return losses[self.criterion]

    @property
    def predictor(self):
        if self.is_clf:
            return self.class_ratios
        elif self.criterion == "absolute_error":
            return self.weighted_median
        else:
            return self.weighted_mean

    def compute_child_loss(self, y: np.ndarray, w: np.ndarray):
        if y.size == 0:
            return np.inf
        pred_dim = (y.size, self.n_classes) if self.is_clf else (y.size,)
        y_pred = np.empty(pred_dim)
        y_pred[:] = self.predictor(y, w)
        return w.sum() * self.loss(y, y_pred, sample_weight=w)

    def compute_split_loss(self, x, y, w, threshold, missing_left=False):
        mask = x < threshold
        if missing_left:
            mask |= np.isnan(x)
        if self.criterion == "friedman_mse":
            diff = self.weighted_mean(y[mask], w[mask]) - self.weighted_mean(
                y[~mask], w[~mask]
            )
            return (-(diff**2) * w[mask].sum() * w[~mask].sum() / w.sum(),)
        return (
            self.compute_child_loss(y[mask], w[mask]),
            self.compute_child_loss(y[~mask], w[~mask]),
        )

    def compute_all_losses(self, x, y, w, missing_left=False):
        nan_mask = np.isnan(x)
        xu = np.unique(x[~nan_mask], sorted=True)
        thresholds = (xu[1:] + xu[:-1]) / 2
        if nan_mask.any() and not missing_left:
            thresholds = np.append(thresholds, xu.max() * 2)
        return thresholds, [
            sum(self.compute_split_loss(x, y, w, threshold, missing_left))
            for threshold in thresholds
        ]

    def best_split_naive(self, X, y, w):
        splits = []
        for f in range(X.shape[1]):
            thresholds, losses = self.compute_all_losses(X[:, f], y, w)
            if self.with_nans:
                thresholds_, losses_ = self.compute_all_losses(
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
                    self.with_nans and idx >= thresholds.size // 2,
                    f,
                )
            )
        return min(splits)


def sparsify(X, density):
    X -= 0.5
    th_low = np.quantile(X.ravel(), q=density / 2)
    th_up = np.quantile(X.ravel(), q=1 - density / 2)
    X[(th_low < X) & (X < th_up)] = 0
    return csc_array(X)


def make_simple_dataset(
    n, d, with_nans, is_sparse, is_clf, n_classes, rng: np.random.Generator
):
    X_dense = np.random.rand(n, d)
    y = np.random.rand(n) + X_dense.sum(axis=1)
    w = np.random.rand(n)

    with_duplicates = rng.integers(2) == 0
    if with_duplicates:
        X_dense = X_dense.round(1 if n < 50 else 2)
    if with_nans:
        for i in range(d):
            step = rng.integers(2, 10)
            X_dense[i::step, i] = np.nan
    if is_sparse:
        density = rng.uniform(0.05, 0.99)
        X = sparsify(X_dense, density)
    else:
        X = X_dense

    if is_clf:
        q = np.linspace(0, 1, num=n_classes + 1)[1:-1]
        y = np.searchsorted(np.quantile(y, q), y)

    return X_dense, X, y, w


@pytest.mark.parametrize("sparse", ["x", "sparse"])
@pytest.mark.parametrize("missing_values", ["x", "missing_values"])
@pytest.mark.parametrize(
    "criterion",
    ["gini", "log_loss", "squared_error", "absolute_error", "friedman_mse", "poisson"],
)
def test_best_split_optimality(sparse, missing_values, criterion, global_random_seed):
    is_clf = criterion in CLF_CRITERIONS
    with_nans = missing_values != "x"
    is_sparse = sparse != "x"
    if is_sparse and with_nans:
        pytest.skip("Sparse + missing values not supported yet")

    rng = np.random.default_rng(global_random_seed)

    d = 2
    for it, n in enumerate([5] * 10 + [10] * 10 + [30] * 5 + [100, 100, 200]):
        n_classes = rng.integers(2, 5)
        X_dense, X, y, w = make_simple_dataset(
            n, d, with_nans, is_sparse, is_clf, n_classes, rng
        )

        naive_splitter = NaiveSplitter(is_clf, criterion, with_nans, n_classes)

        Tree = DecisionTreeClassifier if is_clf else DecisionTreeRegressor
        tree = Tree(criterion=criterion, max_depth=1, max_features=d)
        tree.fit(X, y, sample_weight=w)
        best_split = naive_splitter.best_split_naive(X_dense, y, w)
        tree_loss = naive_splitter.compute_split_loss(
            X_dense[:, tree.tree_.feature[0]],
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
