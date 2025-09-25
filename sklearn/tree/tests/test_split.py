from dataclasses import dataclass
from functools import cached_property
from itertools import chain, combinations
from operator import itemgetter

import numpy as np
import pytest
from scipy.sparse import csc_array
from scipy.special import xlogy

from sklearn.metrics import (
    mean_absolute_error,
    mean_poisson_deviance,
    mean_squared_error,
)
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor

CLF_CRITERIONS = ("gini", "log_loss")
REG_CRITERIONS = ("squared_error", "absolute_error", "friedman_mse", "poisson")


def powerset(iterable):
    s = list(iterable)  # allows handling sets too
    return chain.from_iterable(
        (list(c) for c in combinations(s, r)) for r in range(1, (len(s) + 1) // 2 + 1)
    )


@dataclass
class NaiveSplitter:
    is_clf: bool
    criterion: str
    with_nans: bool
    n_classes: int
    is_categorical: np.ndarray

    @staticmethod
    def weighted_median(y, w):
        sorter = np.argsort(y)
        wc = np.cumsum(w[sorter])
        idx = np.searchsorted(wc, wc[-1] / 2)
        return y[sorter[idx]]

    @staticmethod
    def log_loss(y, y_pred, sample_weight):
        """the one from sklearn.metrics is too slow, due to input-validation"""
        eps = np.finfo(y_pred.dtype).eps
        y_pred = np.clip(y_pred, eps, 1 - eps)
        y_aligned = np.zeros_like(y_pred)
        y_aligned[np.arange(len(y)), y] = 1
        loss = -xlogy(y_aligned, y_pred).sum(axis=1)
        return np.average(loss, weights=sample_weight)

    @staticmethod
    def gini_loss(y, y_pred, sample_weight):
        p = y_pred[0]
        return (p * (1 - p)).sum()

    @cached_property
    def loss(self):
        losses = {
            "poisson": mean_poisson_deviance,
            "squared_error": mean_squared_error,
            "absolute_error": mean_absolute_error,
            "log_loss": self.log_loss,
            "gini": self.gini_loss,
        }
        return losses[self.criterion]

    def class_ratios(self, y, w):
        return np.clip(np.bincount(y, w, minlength=self.n_classes) / w.sum(), None, 1)

    @cached_property
    def predictor(self):
        if self.is_clf:
            return self.class_ratios
        elif self.criterion == "absolute_error":
            return self.weighted_median
        else:
            return lambda y, w: np.average(y, weights=w)

    def compute_child_loss(self, y: np.ndarray, w: np.ndarray):
        if y.size == 0:
            return np.inf
        pred_dim = (y.size, self.n_classes) if self.is_clf else (y.size,)
        y_pred = np.empty(pred_dim)
        y_pred[:] = self.predictor(y, w)
        return w.sum() * self.loss(y, y_pred, sample_weight=w)

    def compute_split_loss(
        self, x, y, w, threshold=None, categories=None, missing_left=False
    ):
        if categories is not None:
            mask_c = np.zeros(int(x.max() + 1), dtype=bool)
            mask_c[categories] = True
            mask = mask_c[x.astype(int)]
        else:
            mask = x < threshold
        if missing_left:
            mask |= np.isnan(x)
        if self.criterion == "friedman_mse":
            diff = np.average(y[mask], weights=w[mask]) - np.average(
                y[~mask], weights=w[~mask]
            )
            return (-(diff**2) * w[mask].sum() * w[~mask].sum() / w.sum(),)
        return (
            self.compute_child_loss(y[mask], w[mask]),
            self.compute_child_loss(y[~mask], w[~mask]),
        )

    def compute_all_losses(self, x, y, w, is_categorical=False, missing_left=False):
        if is_categorical:
            return self.compute_all_losses_categorical(x, y, w)
        nan_mask = np.isnan(x)
        xu = np.unique(x[~nan_mask], sorted=True)
        thresholds = (xu[1:] + xu[:-1]) / 2
        if nan_mask.any() and not missing_left:
            thresholds = np.append(thresholds, xu.max() * 2)
        return thresholds, [
            sum(self.compute_split_loss(x, y, w, threshold, missing_left=missing_left))
            for threshold in thresholds
        ]

    def compute_all_losses_categorical(self, x, y, w):
        cat_splits = list(powerset(np.unique(x).astype(int)))
        return cat_splits, [
            sum(self.compute_split_loss(x, y, w, categories=left_cat))
            for left_cat in cat_splits
        ]

    def best_split_naive(self, X, y, w):
        splits = []
        for f in range(X.shape[1]):
            thresholds, losses = self.compute_all_losses(
                X[:, f], y, w, is_categorical=self.is_categorical[f]
            )
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
        return min(splits, key=itemgetter(0))


def sparsify(X, density):
    X -= 0.5
    th_low = np.quantile(X.ravel(), q=density / 2)
    th_up = np.quantile(X.ravel(), q=1 - density / 2)
    X[(th_low < X) & (X < th_up)] = 0
    return csc_array(X)


def to_categorical(x, nc, rng: np.random.Generator):
    q = np.linspace(0, 1, num=nc + 1)[1:-1]
    quantiles = np.quantile(x, q)
    cats = np.searchsorted(quantiles, x)
    return rng.permutation(nc)[cats]


def make_simple_dataset(
    n,
    d,
    with_nans,
    is_sparse,
    is_categorical,
    is_clf,
    n_classes,
    rng: np.random.Generator,
):
    X_dense = rng.random((n, d))
    y = rng.random(n) + X_dense.sum(axis=1)
    w = rng.random(n)

    for idx in np.where(is_categorical)[0]:
        nc = rng.integers(2, 6)  # cant go to high or test will be too slow
        X_dense[:, idx] = to_categorical(X_dense[:, idx], nc, rng)
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


def bitset_to_set(v: np.uint64):
    return [c for c in range(64) if v & (1 << c)]


@pytest.mark.parametrize("sparse", ["x", "sparse"])
@pytest.mark.parametrize("categorical", ["x"])
@pytest.mark.parametrize("missing_values", ["x", "missing_values"])
@pytest.mark.parametrize(
    "criterion",
    ["gini", "log_loss", "squared_error", "absolute_error", "friedman_mse", "poisson"],
)
def test_best_split_optimality(
    sparse, categorical, missing_values, criterion, global_random_seed
):
    is_clf = criterion in CLF_CRITERIONS
    with_nans = missing_values != "x"
    is_sparse = sparse != "x"
    with_categoricals = categorical != "x"
    if is_sparse and with_nans:
        pytest.skip("Sparse + missing values not supported yet")
    if with_categoricals and (is_sparse or criterion == "absolute_error" or with_nans):
        pytest.skip("Categorical features not supported in this case")

    rng = np.random.default_rng(global_random_seed)

    ns = [5] * 5 + [10] * 5 + [30, 30]
    if not with_categoricals:  # and criterion != "log_loss":
        ns.extend([30, 30, 30, 100, 100, 200])

    for it, n in enumerate(ns):
        d = rng.integers(1, 4)
        n_classes = 2 if with_categoricals else rng.integers(2, 5)
        if with_categoricals:
            is_categorical = rng.random(d) < 0.5
        else:
            is_categorical = np.zeros(d, dtype=bool)
        X_dense, X, y, w = make_simple_dataset(
            n, d, with_nans, is_sparse, is_categorical, is_clf, n_classes, rng
        )

        naive_splitter = NaiveSplitter(
            is_clf, criterion, with_nans, n_classes, is_categorical
        )
        best_split = naive_splitter.best_split_naive(X_dense, y, w)

        is_categorical = is_categorical if is_categorical.any() else None
        Tree = DecisionTreeClassifier if is_clf else DecisionTreeRegressor
        tree = Tree(
            criterion=criterion,
            max_depth=1,
            max_features=d,
            random_state=global_random_seed,
        )
        tree.fit(X, y, sample_weight=w)

        split_feature = tree.tree_.feature[0]
        split = (
            {"threshold": tree.tree_.threshold[0]}
            if is_categorical is None or not is_categorical[split_feature]
            else {"categories": bitset_to_set(tree.tree_.categorical_bitset[0])}
        )
        tree_loss = naive_splitter.compute_split_loss(
            X_dense[:, tree.tree_.feature[0]],
            y,
            w,
            **split,
            missing_left=bool(tree.tree_.missing_go_to_left[0]),
        )
        tree_split = (
            sum(tree_loss),
            split,
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
