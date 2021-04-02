# Authors: Jasper Roebroek <roebroek.jasper@gmail.com>
# License: BSD 3 clause

"""
This module is inspired on the skgarden implementation of Forest Quantile Regression,
based on the following paper:

Nicolai Meinshausen, Quantile Regression Forests
http://www.jmlr.org/papers/volume7/meinshausen06a/meinshausen06a.pdf

Two implementations are available:
- based on the original paper (_DefaultForestQuantileRegressor)
- based on the adapted implementation in quantregForest,
  which provided substantial speed improvements
  (_RandomSampleForestQuantileRegressor)

Two algorithms for fitting are implemented (which are broadcasted)
- Random forest (RandomForestQuantileRegressor)
- Extra Trees (ExtraTreesQuantileRegressor)

RandomForestQuantileRegressor and ExtraTreesQuantileRegressor are therefore only
placeholders that link to the two implementations, passing on a parameter base_estimator
to pick the right training algorithm.
"""
from types import MethodType
from abc import ABCMeta, abstractmethod

import numpy as np
from numba import jit, float32, float64, int64, prange
from numpy.lib.function_base import _quantile_is_valid

from ..tree import DecisionTreeRegressor, ExtraTreeRegressor
from ..utils import check_array, check_X_y, check_random_state
from ._forest import ForestRegressor
from ._forest import _generate_sample_indices

__all__ = ["RandomForestQuantileRegressor", "ExtraTreesQuantileRegressor"]


@jit(float32[:](float32[:, :], float32, float32[:]), nopython=True)
def _weighted_quantile(a, q, weights):
    """
    Weighted quantile calculation.

    Parameters
    ----------
    a : array, shape = (n_sample, n_features)
        Data from which the quantiles are calculated. One quantile value
        per feature (n_features) is given. Should be float32.
    q : float
        Quantile in range [0, 1]. Should be a float32 value.
    weights : array, shape = (n_sample)
        Weights of each sample. Should be float32

    Returns
    -------
    quantiles : array, shape = (n_features)
        Quantile values

    References
    ----------
    1. https://en.wikipedia.org/wiki/Percentile#The_Weighted_Percentile_method

    Notes
    -----
    Note that weighted_quantile(a, q) is not equivalent to np.quantile(a, q).
    This is because in np.quantile sorted(a)[i] is assumed to be at quantile 0.0,
    while here we assume sorted(a)[i] is given a weight of 1.0 / len(a), hence
    it is at the 1.0 / len(a)th quantile.
    """
    nz = weights != 0
    a = a[nz]
    weights = weights[nz]

    n_features = a.shape[1]
    quantiles = np.full(n_features, np.nan, dtype=np.float32)
    if a.shape[0] == 1 or a.size == 0:
        return a[0]

    for i in range(n_features):
        sorted_indices = np.argsort(a[:, i])
        sorted_a = a[sorted_indices, i]
        sorted_weights = weights[sorted_indices]

        # Step 1
        sorted_cum_weights = np.cumsum(sorted_weights)
        total = sorted_cum_weights[-1]

        # Step 2
        partial_sum = 1 / total * (sorted_cum_weights - sorted_weights / 2.0)
        start = np.searchsorted(partial_sum, q) - 1
        if start == len(sorted_cum_weights) - 1:
            quantiles[i] = sorted_a[-1]
            continue
        if start == -1:
            quantiles[i] = sorted_a[0]
            continue

        # Step 3.
        fraction = (q - partial_sum[start]) / (partial_sum[start + 1] - partial_sum[start])
        quantiles[i] = sorted_a[start] + fraction * (sorted_a[start + 1] - sorted_a[start])
    return quantiles


def weighted_quantile(a, q, weights=None):
    """
    Returns the weighted quantile of a at q given weights.

    Parameters
    ----------
    a: array-like, shape=(n_samples, n_features)
        Samples from which the quantile is calculated

    q: float
        Quantile (in the range from 0-1)

    weights: array-like, shape=(n_samples,)
        Weights[i] is the weight given to point a[i] while computing the
        quantile. If weights[i] is zero, a[i] is simply ignored during the
        quantile computation.

    Returns
    -------
    quantile: array, shape = (n_features)
        Weighted quantile of a at q.

    References
    ----------
    1. https://en.wikipedia.org/wiki/Percentile#The_Weighted_Percentile_method

    Notes
    -----
    Note that weighted_quantile(a, q) is not equivalent to np.quantile(a, q).
    This is because in np.quantile sorted(a)[i] is assumed to be at quantile 0.0,
    while here we assume sorted(a)[i] is given a weight of 1.0 / len(a), hence
    it is at the 1.0 / len(a)th quantile.
    """
    if q > 1 or q < 0:
        raise ValueError("q should be in-between 0 and 1, "
                         "got %d" % q)

    a = np.asarray(a, dtype=np.float32)
    if a.ndim == 1:
        a = a.reshape((-1, 1))
    elif a.ndim > 2:
        raise ValueError("a should be in the format (n_samples, n_feature)")

    if weights is None:
        weights = np.ones(a.shape[0], dtype=np.float32)
    else:
        weights = np.asarray(weights, dtype=np.float32)
        if weights.ndim > 1:
            raise ValueError("weights need to be 1 dimensional")

    if a.shape[0] != weights.shape[0]:
        raise ValueError("a and weights should have the same length.")

    q = np.float32(q)

    quantiles = _weighted_quantile(a, q, weights)

    if quantiles.size == 1:
        return quantiles[0]
    else:
        return quantiles


@jit(float32[:, :](int64[:, :], float32[:, :], int64[:, :], float32[:, :], float32), parallel=True, nopython=True)
def _quantile_forest_predict(X_leaves, y_train, y_train_leaves, y_weights, q):
    quantiles = np.zeros((X_leaves.shape[0], y_train.shape[1]), dtype=np.float32)
    for i in prange(len(X_leaves)):
        x_leaf = X_leaves[i]
        x_weights = np.zeros(y_weights.shape[1], dtype=np.float32)
        for j in range(y_weights.shape[1]):
            x_weights[j] = (y_weights[:, j] * (y_train_leaves[:, j] == x_leaf)).sum()
        quantiles[i] = _weighted_quantile(y_train, q, x_weights)
    return quantiles


@jit(int64[:](int64[:], int64[:], float64[:], int64[:], float64[:]), parallel=True, nopython=True)
def _weighted_random_sample(leaves, unique_leaves, weights, idx, random_numbers):
    """
    Random sample for each unique leaf

    Parameters
    ----------
    leaves : array, shape = (n_samples)
        Leaves of a Regression tree, corresponding to weights and indices (idx)
    weights : array, shape = (n_samples)
        Weights for each observation. They need to sum up to 1 per unique leaf.
    idx : array, shape = (n_samples)
        Indices of original observations. The output will drawn from this.

    Returns
    -------
    unique_leaves, sampled_idx, shape = (n_unique_samples)
        Unique leaves (from 'leaves') and a randomly (and weighted) sample
        from 'idx' corresponding to the leaf.
    """
    sampled_idx = np.empty_like(unique_leaves, dtype=np.int64)

    for i in prange(len(unique_leaves)):
        mask = unique_leaves[i] == leaves
        c_weights = weights[mask]
        c_idx = idx[mask]

        if c_idx.size == 1:
            sampled_idx[i] = c_idx[0]
            continue

        p = 0
        r = random_numbers[i]
        for j in range(len(c_idx)):
            p += c_weights[j]
            if p > r:
                sampled_idx[i] = c_idx[j]
                break

    return sampled_idx


class _ForestQuantileRegressor(ForestRegressor, metaclass=ABCMeta):
    """
    A forest regressor providing quantile estimates.

    The generation of the forest can be either based on Random Forest or
    Extra Trees algorithms. The fitting and prediction of the forest can
    be based on the methods layed out in the original paper of Meinshausen,
    or on the adapted implementation of the R quantregForest package.

    Parameters
    ----------
    n_estimators : integer, optional (default=10)
        The number of trees in the forest.

    criterion : string, optional (default="mse")
        The function to measure the quality of a split. Supported criteria
        are "mse" for the mean squared error, which is equal to variance
        reduction as feature selection criterion, and "mae" for the mean
        absolute error.
        .. versionadded:: 0.18
           Mean Absolute Error (MAE) criterion.

    max_features : int, float, string or None, optional (default="auto")
        The number of features to consider when looking for the best split:
        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a percentage and
          `int(max_features * n_features)` features are considered at each
          split.
        - If "auto", then `max_features=n_features`.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None, then `max_features=n_features`.
        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int, float, optional (default=2)
        The minimum number of samples required to split an internal node:
        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a percentage and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.
        .. versionchanged:: 0.18
           Added float values for percentages.

    min_samples_leaf : int, float, optional (default=1)
        The minimum number of samples required to be at a leaf node:
        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a percentage and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.
        .. versionchanged:: 0.18
           Added float values for percentages.

    min_weight_fraction_leaf : float, optional (default=0.)
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_leaf_nodes : int or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    oob_score : bool, optional (default=False)
        whether to use out-of-bag samples to estimate
        the R^2 on unseen data.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel for both `fit` and `predict`.
        If -1, then the number of jobs is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controls the verbosity of the tree building process.

    warm_start : bool, optional (default=False)
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit a whole
        new forest.

    base_estimator : ``DecisionTreeRegressor``, optional
        Subclass of ``DecisionTreeRegressor`` as the base_estimator for the
        generation of the forest. Either DecisionTreeRegressor() or ExtraTreeRegressor().


    Attributes
    ----------
    estimators_ : list of DecisionTreeRegressor
        The collection of fitted sub-estimators.

    feature_importances_ : array of shape = [n_features]
        The feature importances (the higher, the more important the feature).

    n_features_ : int
        The number of features when ``fit`` is performed.

    n_outputs_ : int
        The number of outputs when ``fit`` is performed.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.

    oob_prediction_ : array of shape = [n_samples]
        Prediction computed with out-of-bag estimate on the training set.

    References
    ----------
    .. [1] Nicolai Meinshausen, Quantile Regression Forests
        http://www.jmlr.org/papers/volume7/meinshausen06a/meinshausen06a.pdf
    """
    # allowed options
    methods = ['default', 'sample']
    base_estimators = ['random_forest', 'extra_trees']

    def __init__(self,
                 n_estimators=10,
                 criterion='mse',
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_weight_fraction_leaf=0.0,
                 max_features='auto',
                 max_leaf_nodes=None,
                 bootstrap=True,
                 oob_score=False,
                 n_jobs=1,
                 random_state=None,
                 verbose=0,
                 warm_start=False,
                 base_estimator=DecisionTreeRegressor()):
        super(_ForestQuantileRegressor, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_samples_split",
                              "min_samples_leaf", "min_weight_fraction_leaf",
                              "max_features", "max_leaf_nodes",
                              "random_state"),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes

    @abstractmethod
    def fit(self, X, y, sample_weight):
        """
        Build a forest from the training set (X, y).

        Parameters
        ----------
        X : array-like or sparse matrix, shape = (n_samples, n_features)
            The training input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csc_matrix``.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            The target values

        sample_weight : array-like, shape = (n_samples) or None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. Splits are also
            ignored if they would result in any single class carrying a
            negative weight in either child node.

        Returns
        -------
        self : object
            Returns self.
        """
        raise NotImplementedError("This class is not meant of direct construction, the fitting method should be "
                                  "obtained from either _DefaultForestQuantileRegressor or "
                                  "_RandomSampleForestQuantileRegressor")

    @abstractmethod
    def predict(self, X, q=0.5):
        """
        Predict quantile regression values for X.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        q : float, optional
            Value ranging from 0 to 1. By default, the median is predicted

        Returns
        -------
        y : array of shape = (n_samples) or (n_samples, n_outputs)
            return y such that F(Y=y | x) = quantile.
        """
        raise NotImplementedError("This class is not meant of direct construction, the prediction method should be "
                                  "obtained from either _DefaultForestQuantileRegressor or "
                                  "_RandomSampleForestQuantileRegressor")

    def repr(self, method):
        s = super(_ForestQuantileRegressor, self).__repr__()

        if type(self.base_estimator) is DecisionTreeRegressor:
            c = "RandomForestQuantileRegressor"
        elif type(self.base_estimator) is ExtraTreeRegressor:
            c = "ExtraTreesQuantileRegressor"

        params = s[s.find("(") + 1:s.rfind(")")].split(", ")
        params.append(f"method='{method}'")
        params = [x for x in params if x[:14] != "base_estimator"]

        return f"{c}({', '.join(params)})"


class _DefaultForestQuantileRegressor(_ForestQuantileRegressor):
    """
    fit and predict functions for forest quantile regressors based on:
    Nicolai Meinshausen, Quantile Regression Forests
        http://www.jmlr.org/papers/volume7/meinshausen06a/meinshausen06a.pdf

    todo; if y is 1D, the attributes might be presorted. This could speed up the processes a lot
    """
    def fit(self, X, y, sample_weight=None):
        # apply method requires X to be of dtype np.float32
        X, y = check_X_y(
            X, y, accept_sparse="csc", dtype=np.float32, multi_output=True)
        super(_ForestQuantileRegressor, self).fit(X, y, sample_weight=sample_weight)

        self.n_samples_ = len(y)
        self.y_train_ = y.reshape((-1, self.n_outputs_)).astype(np.float32)
        self.y_train_leaves_ = self.apply(X).T
        self.y_weights_ = np.zeros_like(self.y_train_leaves_, dtype=np.float32)

        if sample_weight is None:
            sample_weight = np.ones(self.n_samples_)

        # todo; parallelization
        for i, est in enumerate(self.estimators_):
            if self.bootstrap:
                bootstrap_indices = _generate_sample_indices(
                    est.random_state, self.n_samples_, self.n_samples_)
            else:
                bootstrap_indices = np.arange(self.n_samples_)

            weights = sample_weight * np.bincount(bootstrap_indices, minlength=self.n_samples_)
            self.y_weights_[i] = weights / est.tree_.weighted_n_node_samples[self.y_train_leaves_[i]]

        self.y_train_leaves_[self.y_weights_ == 0] = -1
        return self

    def predict(self, X, q):
        if not 0 <= q <= 1:
            raise ValueError("q should be between 0 and 1")

        # apply method requires X to be of dtype np.float32
        X = check_array(X, dtype=np.float32, accept_sparse="csc")

        X_leaves = self.apply(X)
        return _quantile_forest_predict(X_leaves, self.y_train_, self.y_train_leaves_, self.y_weights_, q).squeeze()

    def __repr__(self):
        return super(_DefaultForestQuantileRegressor, self).repr(method='default')


class _RandomSampleForestQuantileRegressor(_ForestQuantileRegressor):
    """
    fit and predict functions for forest quantile regressors. Implementation based on quantregForest R packakge.
    """
    def fit(self, X, y, sample_weight=None):
        # apply method requires X to be of dtype np.float32
        X, y = check_X_y(
            X, y, accept_sparse="csc", dtype=np.float32, multi_output=True)
        super(_ForestQuantileRegressor, self).fit(X, y, sample_weight=sample_weight)

        self.n_samples_ = len(y)
        y = y.reshape((-1, self.n_outputs_))

        if sample_weight is None:
            sample_weight = np.ones(self.n_samples_)

        # todo; parallelisation
        for i, est in enumerate(self.estimators_):
            if self.verbose:
                print(f"Sampling tree {i}")

            if self.bootstrap:
                bootstrap_indices = _generate_sample_indices(
                    est.random_state, self.n_samples_, self.n_samples_)
            else:
                bootstrap_indices = np.arange(self.n_samples_)

            y_weights = np.bincount(bootstrap_indices, minlength=self.n_samples_) * sample_weight
            mask = y_weights > 0

            leaves = est.apply(X[mask])
            idx = np.arange(len(y), dtype=np.int64)[mask]

            weights = y_weights[mask] / est.tree_.weighted_n_node_samples[leaves]
            unique_leaves = np.unique(leaves)

            random_instance = check_random_state(est.random_state)
            random_numbers = random_instance.rand(len(unique_leaves))

            sampled_idx = _weighted_random_sample(leaves, unique_leaves, weights, idx, random_numbers)

            est.tree_.value[unique_leaves, :, 0] = y[sampled_idx]

        return self

    def predict(self, X, q):
        q = np.atleast_1d(q)
        if not _quantile_is_valid(q):
            raise ValueError("Quantiles must be in the range [0, 1]")

        if q.ndim > 2:
            raise ValueError("q must be a scalar or 1D")

        # apply method requires X to be of dtype np.float32
        X = check_array(X, dtype=np.float32, accept_sparse="csc")

        quantiles = np.empty((len(X), self.n_outputs_, self.n_estimators))

        # todo; parallelisation
        for i, est in enumerate(self.estimators_):
            if self.n_outputs_ == 1:
                quantiles[:, 0, i] = est.predict(X)
            else:
                quantiles[:, :, i] = est.predict(X)

        return np.quantile(quantiles, q=q, axis=-1).squeeze()

    def __repr__(self):
        return super(_RandomSampleForestQuantileRegressor, self).repr(method='sample')


class RandomForestQuantileRegressor:
    def __new__(cls, *, method='default', **kwargs):
        if method == 'default':
            return _DefaultForestQuantileRegressor(**kwargs)
        elif method == 'sample':
            return _RandomSampleForestQuantileRegressor(**kwargs)


class ExtraTreesQuantileRegressor:
    def __new__(cls, *, method='default', **kwargs):
        if method == 'default':
            return _DefaultForestQuantileRegressor(base_estimator=ExtraTreeRegressor(), **kwargs)
        elif method == 'sample':
            return _RandomSampleForestQuantileRegressor(base_estimator=ExtraTreeRegressor(), **kwargs)
