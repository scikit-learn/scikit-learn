# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Authors: Jasper Roebroek <roebroek.jasper@gmail.com>
# License: BSD 3 clause

"""
This module is inspired on the skgarden implementation of Forest Quantile Regression,
based on the following paper:

Nicolai Meinshausen, Quantile Regression Forests
http://www.jmlr.org/papers/volume7/meinshausen06a/meinshausen06a.pdf

Two implementations are available:
- based on the original paper (_DefaultForestQuantileRegressor). Suitable up to around 100.000 test samples.
- based on the adapted implementation in quantregForest (R), which provides substantial speed improvements
  (_RandomSampleForestQuantileRegressor). Is to be prefered above 100.000 test samples.

Two algorithms for fitting are implemented (which are broadcasted)
- Random forest (RandomForestQuantileRegressor)
- Extra Trees (ExtraTreesQuantileRegressor)

RandomForestQuantileRegressor and ExtraTreesQuantileRegressor are therefore only
placeholders that link to the two implementations, passing on a parameter base_estimator
to pick the right training algorithm.
"""
from abc import ABCMeta, abstractmethod

from cython.parallel cimport prange, parallel
cimport openmp
cimport numpy as np
from numpy cimport ndarray
from ..utils._weighted_quantile cimport _weighted_quantile_presorted_1D, _weighted_quantile_unchecked_1D, Interpolation

import numpy as np
from numpy.lib.function_base import _quantile_is_valid

import threading
from joblib import Parallel

from ._forest import ForestRegressor, _accumulate_prediction, _generate_sample_indices
from ._base import _partition_estimators
from ..utils.fixes import delayed
from ..utils.fixes import _joblib_parallel_args
from ..tree import DecisionTreeRegressor, ExtraTreeRegressor
from ..utils import check_array, check_X_y, check_random_state
from ..utils.validation import check_is_fitted


__all__ = ["RandomForestQuantileRegressor", "ExtraTreesQuantileRegressor"]


cdef void _quantile_forest_predict(long[:, ::1] X_leaves,
                                   float[:, ::1] y_train,
                                   long[:, ::1] y_train_leaves,
                                   float[:, ::1] y_weights,
                                   float[::1] q,
                                   float[:, :, ::1] quantiles):
    """
    X_leaves : (n_estimators, n_test_samples)
    y_train : (n_samples, n_outputs)
    y_train_leaves : (n_estimators, n_samples)
    y_weights : (n_estimators, n_samples)
    q : (n_q)
    quantiles : (n_q, n_test_samplse, n_outputs)
    
    Notes
    -----
    inspired by:
    https://stackoverflow.com/questions/42281886/cython-make-prange-parallelization-thread-safe
    """
    # todo: potential speedup (according to the article linked in the notes) by padding x_weights and x_a:
    #   "You get a little bit extra performance by avoiding padding the private parts of the array to 64 byte,
    #   which is a typical cache line size.". I am not sure how to deal with this

    cdef:
        int n_estimators = X_leaves.shape[0]
        int n_outputs = y_train.shape[1]
        int n_q = q.shape[0]
        int n_samples = y_train.shape[0]
        int n_test_samples = X_leaves.shape[1]

        int i, j, e, o, tid, count_samples
        float curr_weight
        bint sorted = y_train.shape[1] == 1

        int num_threads = openmp.omp_get_max_threads()
        float[::1] x_weights = np.empty(n_samples * num_threads, dtype=np.float32)
        float[:, ::1] x_a = np.empty((n_samples * num_threads, n_outputs), dtype=np.float32)

    with nogil, parallel():
        tid = openmp.omp_get_thread_num()
        for i in prange(n_test_samples):
            count_samples = 0
            for j in range(n_samples):
                curr_weight = 0
                for e in range(n_estimators):
                    if X_leaves[e, i] == y_train_leaves[e, j]:
                        curr_weight = curr_weight + y_weights[e, j]
                if curr_weight > 0:
                    x_weights[tid * n_samples + count_samples] = curr_weight
                    x_a[tid * n_samples + count_samples] = y_train[j]
                    count_samples = count_samples + 1
            if sorted:
                _weighted_quantile_presorted_1D(x_a[tid * n_samples: tid * n_samples + count_samples, 0],
                                                q, x_weights[tid * n_samples: tid * n_samples + count_samples],
                                                quantiles[:, i, 0], Interpolation.linear)
            else:
                for o in range(n_outputs):
                    with gil:
                        curr_x_weights = x_weights[tid * n_samples: tid * n_samples + count_samples].copy()
                        curr_x_a = x_a[tid * n_samples: tid * n_samples + count_samples, o].copy()
                        _weighted_quantile_unchecked_1D(curr_x_a, q, curr_x_weights, quantiles[:, i, o],
                                                        Interpolation.linear)



cdef void _weighted_random_sample(long[::1] leaves,
                                  long[::1] unique_leaves,
                                  float[::1] weights,
                                  long[::1] idx,
                                  double[::1] random_numbers,
                                  long[::1] sampled_idx):
    """
    Random sample for each unique leaf

    Parameters
    ----------
    leaves : array, shape = (n_samples)
        Leaves of a Regression tree, corresponding to weights and indices (idx)
    unique_leaves : array, shape = (n_unique_leaves)
    
    weights : array, shape = (n_samples)
        Weights for each observation. They need to sum up to 1 per unique leaf.
    idx : array, shape = (n_samples)
        Indices of original observations. The output will drawn from this.
    random numbers : array, shape = (n_unique_leaves)
    
    sampled_idx : shape = (n_unique_leaves)
    """
    cdef:
        long c_leaf
        float p, r
        int i, j

        int n_unique_leaves = unique_leaves.shape[0]
        int n_samples = weights.shape[0]

    for i in prange(n_unique_leaves, nogil=True):
        p = 0
        r = random_numbers[i]
        c_leaf = unique_leaves[i]

        for j in range(n_samples):
            if leaves[j] == c_leaf:
                p = p + weights[j]
                if p > r:
                    sampled_idx[i] = idx[j]
                    break


def _accumulate_prediction(predict, X, i, out, lock):
    """
    From sklearn.ensemble._forest
    """
    prediction = predict(X, check_input=False)
    with lock:
        if out.shape[1] == 1:
            out[:, 0, i] = prediction
        else:
            out[..., i] = prediction


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

    quantiles : array-like, optional
        Value ranging from 0 to 1

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
                 base_estimator,
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
                 quantiles=None):
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
        self.quantiles = quantiles

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
    def predict(self, X):
        """
        Predict quantile regression values for X.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : array of shape = (n_quantiles, n_samples, n_outputs)
            return y such that F(Y=y | x) = quantile. If n_quantiles is 1, the array is reduced to
            (n_samples, n_outputs) and if n_outputs is 1, the array is reduced to (n_samples)
        """
        raise NotImplementedError("This class is not meant of direct construction, the prediction method should be "
                                  "obtained from either _DefaultForestQuantileRegressor or "
                                  "_RandomSampleForestQuantileRegressor")

    def validate_quantiles(self):
        if self.quantiles is None:
            raise AttributeError("Quantiles are not set. Please provide them with `model.quantiles = quantiles`")
        q = np.asarray(self.quantiles, dtype=np.float32)
        q = np.atleast_1d(q)
        if q.ndim > 2:
            raise ValueError("q must be a scalar or 1D")

        if not _quantile_is_valid(q):
            raise ValueError("Quantiles must be in the range [0, 1]")

        return q

    def repr(self, method):
        # not terribly pretty, but it works...
        s = super(_ForestQuantileRegressor, self).__repr__()

        if type(self.base_estimator) is DecisionTreeRegressor:
            c = "RandomForestQuantileRegressor"
        elif type(self.base_estimator) is ExtraTreeRegressor:
            c = "ExtraTreesQuantileRegressor"
        else:
            raise TypeError("base_estimator needs to be either DecisionTreeRegressor or ExtraTreeRegressor")

        params = s[s.find("(") + 1:s.rfind(")")].split(", ")
        params.append(f"method='{method}'")
        params = [x for x in params if x[:14] != "base_estimator"]

        return f"{c}({', '.join(params)})"


class _DefaultForestQuantileRegressor(_ForestQuantileRegressor):
    """
    fit and predict functions for forest quantile regressors based on:
    Nicolai Meinshausen, Quantile Regression Forests
        http://www.jmlr.org/papers/volume7/meinshausen06a/meinshausen06a.pdf
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

        for i, est in enumerate(self.estimators_):
            est.y_train_ = self.y_train_
            est.y_train_leaves_ = self.y_train_leaves_[i]
            est.y_weights_ = self.y_weights_[i]
            est.verbose = self.verbose
            est.n_samples_ = self.n_samples_
            est.bootstrap = self.bootstrap
            est._i = i

        for i, est in enumerate(self.estimators_):
            if self.bootstrap:
                bootstrap_indices = _generate_sample_indices(
                    est.random_state, self.n_samples_, self.n_samples_)
            else:
                bootstrap_indices = np.arange(self.n_samples_)

            weights = sample_weight * np.bincount(bootstrap_indices, minlength=self.n_samples_)
            self.y_weights_[i] = weights / est.tree_.weighted_n_node_samples[self.y_train_leaves_[i]]

        self.y_train_leaves_[self.y_weights_ == 0] = -1

        if self.n_outputs_ == 1:
            sort_ind = np.argsort(y)
            self.y_train_[:] = self.y_train_[sort_ind]
            self.y_weights_[:] = self.y_weights_[:, sort_ind]
            self.y_train_leaves_[:] = self.y_train_leaves_[:, sort_ind]
            self.y_sorted_ = True
        else:
            self.y_sorted_ = False

        return self

    def predict(self, X):
        check_is_fitted(self)
        q = self.validate_quantiles()

        # apply method requires X to be of dtype np.float32
        X = check_array(X, dtype=np.float32, accept_sparse="csc")
        X_leaves = self.apply(X).T
        n_test_samples = X.shape[0]

        quantiles = np.empty((q.size, n_test_samples, self.n_outputs_), dtype=np.float32)

        _quantile_forest_predict(X_leaves, self.y_train_, self.y_train_leaves_, self.y_weights_,
                                 q, quantiles)
        if q.size == 1:
            quantiles = quantiles[0]
        if self.n_outputs_ == 1:
            quantiles = quantiles[..., 0]

        return quantiles

    def __repr__(self):
        return super(_DefaultForestQuantileRegressor, self).repr(method='default')


class _RandomSampleForestQuantileRegressor(_DefaultForestQuantileRegressor):
    """
    fit and predict functions for forest quantile regressors. Implementation based on quantregForest R packakge.
    """
    def fit(self, X, y, sample_weight=None):
        super(_RandomSampleForestQuantileRegressor, self).fit(X, y, sample_weight=sample_weight)

        for i, est in enumerate(self.estimators_):
            if self.verbose:
                print(f"Sampling tree {i} of {self.n_estimators}")

            mask = est.y_weights_ > 0

            leaves = est.y_train_leaves_[mask]
            idx = np.arange(self.n_samples_)[mask]

            unique_leaves = np.unique(leaves)

            random_instance = check_random_state(est.random_state)
            random_numbers = random_instance.rand(len(unique_leaves))

            sampled_idx = np.empty(len(unique_leaves), dtype=np.int64)
            _weighted_random_sample(leaves, unique_leaves, est.y_weights_[mask], idx, random_numbers, sampled_idx)

            est.tree_.value[unique_leaves, :, 0] = self.y_train_[sampled_idx]

        return self

    def predict(self, X):
        check_is_fitted(self)
        q = self.validate_quantiles()

        # Assign chunk of trees to jobs
        n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # apply method requires X to be of dtype np.float32
        X = check_array(X, dtype=np.float32, accept_sparse="csc")

        predictions = np.empty((len(X), self.n_outputs_, self.n_estimators))

        lock = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose,
                 **_joblib_parallel_args(require="sharedmem"))(
            delayed(_accumulate_prediction)(est.predict, X, i, predictions, lock)
            for i, est in enumerate(self.estimators_))

        quantiles = np.quantile(predictions, q=q, axis=-1)
        if q.size == 1:
            quantiles = quantiles[0]
        if self.n_outputs_ == 1:
            quantiles = quantiles[..., 0]

        return quantiles

    def __repr__(self):
        return super(_RandomSampleForestQuantileRegressor, self).repr(method='sample')


class RandomForestQuantileRegressor:
    def __new__(cls, *, method='default', **kwargs):
        if method == 'default':
            return _DefaultForestQuantileRegressor(base_estimator=DecisionTreeRegressor(), **kwargs)
        elif method == 'sample':
            return _RandomSampleForestQuantileRegressor(base_estimator=DecisionTreeRegressor(), **kwargs)


class ExtraTreesQuantileRegressor:
    def __new__(cls, *, method='default', **kwargs):
        if method == 'default':
            return _DefaultForestQuantileRegressor(base_estimator=ExtraTreeRegressor(), **kwargs)
        elif method == 'sample':
            return _RandomSampleForestQuantileRegressor(base_estimator=ExtraTreeRegressor(), **kwargs)
