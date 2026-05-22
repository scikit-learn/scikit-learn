from dataclasses import dataclass
from itertools import product
from operator import itemgetter

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.sparse import csc_array
from scipy.special import xlogy

from sklearn.metrics import mean_poisson_deviance
from sklearn.tree import (
    DecisionTreeClassifier,
    DecisionTreeRegressor,
    ExtraTreeClassifier,
    ExtraTreeRegressor,
)
from sklearn.utils.stats import _weighted_percentile

CLF_CRITERIONS = ("gini", "log_loss")

REG_CRITERIONS = ("squared_error", "absolute_error", "poisson")

CLF_TREES = {
    "DecisionTreeClassifier": DecisionTreeClassifier,
    "ExtraTreeClassifier": ExtraTreeClassifier,
}

REG_TREES = {
    "DecisionTreeRegressor": DecisionTreeRegressor,
    "ExtraTreeRegressor": ExtraTreeRegressor,
}


@dataclass
class NaiveSplitter:
    criterion: str
    n_classes: int = 0

    def compute_node_value_and_impurity(self, y, w):
        sum_weights = np.sum(w)
        if sum_weights < 1e-7:
            return np.nan, np.inf  # invalid split
        if self.criterion in ["gini", "entropy", "log_loss"]:
            pred = np.bincount(y, weights=w, minlength=self.n_classes) / sum_weights
            if self.criterion == "gini":
                # 1 - sum(pk^2)
                loss = 1.0 - np.sum(pred**2)
            else:
                # -sum(pk * log2(pk))
                loss = -np.sum(xlogy(pred, pred)) / np.log(2)
        elif self.criterion == "squared_error":
            pred = np.average(y, weights=w)
            loss = np.average((y - pred) ** 2, weights=w)
        elif self.criterion == "absolute_error":
            pred = _weighted_percentile(y, w, percentile_rank=50, average=True)
            loss = np.average(np.abs(y - pred), weights=w)
        elif self.criterion == "poisson":
            pred = np.average(y, weights=w)
            loss = mean_poisson_deviance(y, np.repeat(pred, y.size), sample_weight=w)
            loss *= 1 / 2
        else:
            raise ValueError(f"Unknown criterion: {self.criterion}")
        return pred, loss * sum_weights

    def compute_split_nodes(self, X, y, w, feature, threshold=None, missing_left=False):
        x = X[:, feature]
        go_left = x <= threshold
        if missing_left:
            go_left |= np.isnan(x)
        return (
            self.compute_node_value_and_impurity(y[go_left], w[go_left]),
            self.compute_node_value_and_impurity(y[~go_left], w[~go_left]),
        )

    def compute_split_impurity(
        self, X, y, w, feature, threshold=None, missing_left=False
    ):
        nodes = self.compute_split_nodes(X, y, w, feature, threshold, missing_left)
        (_, left_impurity), (_, right_impurity) = nodes
        return left_impurity + right_impurity

    def _generate_all_splits(self, X):
        for f in range(X.shape[1]):
            x = X[:, f]
            nan_mask = np.isnan(x)
            thresholds = np.unique(x[~nan_mask])
            for th in thresholds:
                yield {
                    "feature": f,
                    "threshold": th,
                    "missing_left": False,
                }
            if not nan_mask.any():
                continue
            for th in [*thresholds, -np.inf]:
                # include -inf to test the split with only NaNs on the left node
                yield {
                    "feature": f,
                    "threshold": th,
                    "missing_left": True,
                }

    def best_split_naive(self, X, y, w):
        splits = list(self._generate_all_splits(X))
        if len(splits) == 0:
            return (np.inf, None)

        split_impurities = [
            self.compute_split_impurity(X, y, w, **split) for split in splits
        ]

        return min(zip(split_impurities, splits), key=itemgetter(0))


def make_simple_dataset(
    n,
    d,
    with_nans,
    is_sparse,
    is_clf,
    n_classes,
    rng,
):
    X_dense = rng.random((n, d))
    y = rng.random(n) + X_dense.sum(axis=1)
    w = rng.integers(0, 5, size=n) if rng.uniform() < 0.5 else rng.random(n)

    with_duplicates = rng.integers(2) == 0
    if with_duplicates:
        X_dense = X_dense.round(1 if n < 50 else 2)
    if with_nans:
        nan_density = rng.uniform(0.05, 0.8)
        mask = rng.random(X_dense.shape) < nan_density
        X_dense[mask] = np.nan
    if is_sparse:
        density = rng.uniform(0.05, 0.99)
        X_dense -= 0.5
        mask = rng.random(X_dense.shape) > density
        X_dense[mask] = 0
        X = csc_array(X_dense)
    else:
        X = X_dense

    if is_clf:
        q = np.linspace(0, 1, num=n_classes + 1)[1:-1]
        y = np.searchsorted(np.quantile(y, q), y)

    # Trees cast X to float32 internally; match that dtype here to avoid
    # routing/impurity mismatches from rounding with `<=`.
    return X_dense.astype("float32"), X, y, w


@pytest.mark.filterwarnings("ignore:.*friedman_mse.*:FutureWarning")
@pytest.mark.parametrize(
    "Tree, criterion",
    [
        *product(REG_TREES.values(), REG_CRITERIONS),
        *product(CLF_TREES.values(), CLF_CRITERIONS),
    ],
)
@pytest.mark.parametrize(
    "sparse, missing_values",
    [(False, False), (True, False), (False, True)],
    ids=["dense-without_missing", "sparse-without_missing", "dense-with_missing"],
)
def test_split_impurity(Tree, criterion, sparse, missing_values, global_random_seed):
    is_clf = criterion in CLF_CRITERIONS

    rng = np.random.default_rng(global_random_seed)

    ns = [5] * 5 + [10] * 5 + [20, 30, 50, 100]

    for it, n in enumerate(ns):
        d = rng.integers(1, 4)
        n_classes = rng.integers(2, 5)  # only used for classification
        X_dense, X, y, w = make_simple_dataset(
            n, d, missing_values, sparse, is_clf, n_classes, rng
        )

        naive_splitter = NaiveSplitter(criterion, n_classes)

        tree = Tree(
            criterion=criterion,
            max_depth=1,
            random_state=global_random_seed,
        )
        tree.fit(X, y, sample_weight=w)
        actual_impurity = tree.tree_.impurity * tree.tree_.weighted_n_node_samples
        actual_value = tree.tree_.value[:, 0]

        # Check root's impurity:
        # The root is 0, left child is 1 and right child is 2.
        root_val, root_impurity = naive_splitter.compute_node_value_and_impurity(y, w)
        assert_allclose(root_impurity, actual_impurity[0], atol=1e-12)
        assert_allclose(root_val, actual_value[0], atol=1e-12)

        if tree.tree_.node_count == 1:
            # if no splits was made assert that either:
            assert (
                "Extra" in Tree.__name__
                or root_impurity < 1e-12  # root impurity is 0
                # or no valid split can be made:
                or naive_splitter.best_split_naive(X_dense, y, w)[0] == np.inf
            )
            continue

        # Check children impurity:
        actual_split = {
            "feature": int(tree.tree_.feature[0]),
            "threshold": tree.tree_.threshold[0],
            "missing_left": bool(tree.tree_.missing_go_to_left[0]),
        }
        nodes = naive_splitter.compute_split_nodes(X_dense, y, w, **actual_split)
        (left_val, left_impurity), (right_val, right_impurity) = nodes
        assert_allclose(left_impurity, actual_impurity[1], atol=1e-12)
        assert_allclose(right_impurity, actual_impurity[2], atol=1e-12)
        assert_allclose(left_val, actual_value[1], atol=1e-12)
        assert_allclose(right_val, actual_value[2], atol=1e-12)

        if "Extra" in Tree.__name__:
            # The remainder of the test checks for optimality of the found split.
            # However, randomized trees are not guaranteed to find an optimal split
            # but only a "better-than-nothing" split.
            # Therefore, end the test here for these models.
            continue

        # Check that the selected split has the same impurity as the best split
        # found by the naive splitter. Note that there could exist multiple splits
        # with the same optimal impurity, so the assertion is made on the impurity
        # value: the split value is only displayed to help debugging in case
        # of assertion failure.
        best_impurity, best_split = naive_splitter.best_split_naive(X_dense, y, w)
        actual_split_impurity = actual_impurity[1:].sum()
        assert np.isclose(best_impurity, actual_split_impurity), (
            best_split,
            actual_split,
        )
