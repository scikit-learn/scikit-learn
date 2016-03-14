# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

from __future__ import division

import numbers
import numpy as np
from warnings import warn

from scipy.sparse import issparse

from ..externals import six
from ..externals.joblib import Parallel, delayed
from ..tree import ExtraTreeRegressor
from ..utils import check_random_state, check_array

from .bagging import BaseBagging
from .base import _partition_estimators

__all__ = ["IsolationForest"]


class IsolationForest(BaseBagging):
    """Isolation Forest Algorithm

    Return the anomaly score of each sample using the IsolationForest algorithm

    The IsolationForest 'isolates' observations by randomly selecting a feature
    and then randomly selecting a split value between the maximum and minimum
    values of the selected feature.

    Since recursive partitioning can be represented by a tree structure, the
    number of splittings required to isolate a sample is equivalent to the path
    length from the root node to the terminating node.

    This path length, averaged over a forest of such random trees, is a
    measure of abnormality and our decision function.

    Random partitioning produces noticeably shorter paths for anomalies.
    Hence, when a forest of random trees collectively produce shorter path
    lengths for particular samples, they are highly likely to be anomalies.

    Read more in the :ref:`User Guide <isolation_forest>`.

    Parameters
    ----------
    n_estimators : int, optional (default=100)
        The number of base estimators in the ensemble.

    max_samples : int or float, optional (default="auto")
        The number of samples to draw from X to train each base estimator.
            - If int, then draw `max_samples` samples.
            - If float, then draw `max_samples * X.shape[0]` samples.
            - If "auto", then `max_samples=min(256, n_samples)`.
        If max_samples is larger than the number of samples provided,
        all samples will be used for all trees (no sampling).

    max_features : int or float, optional (default=1.0)
        The number of features to draw from X to train each base estimator.
            - If int, then draw `max_features` features.
            - If float, then draw `max_features * X.shape[1]` features.

    bootstrap : boolean, optional (default=False)
        Whether samples are drawn with replacement.

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


    Attributes
    ----------
    estimators_ : list of DecisionTreeClassifier
        The collection of fitted sub-estimators.

    estimators_samples_ : list of arrays
        The subset of drawn samples (i.e., the in-bag samples) for each base
        estimator.

    max_samples_ : integer
        The actual number of samples

    References
    ----------
    .. [1] Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. "Isolation forest."
           Data Mining, 2008. ICDM'08. Eighth IEEE International Conference on.
    .. [2] Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. "Isolation-based
           anomaly detection." ACM Transactions on Knowledge Discovery from
           Data (TKDD) 6.1 (2012): 3.
    """

    def __init__(self,
                 n_estimators=100,
                 max_samples="auto",
                 max_features=1.,
                 bootstrap=False,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):
        super(IsolationForest, self).__init__(
            base_estimator=ExtraTreeRegressor(
                max_features=1,
                splitter='random',
                random_state=random_state),
            # here above max_features has no links with self.max_features
            bootstrap=bootstrap,
            bootstrap_features=False,
            n_estimators=n_estimators,
            max_samples=max_samples,
            max_features=max_features,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def _set_oob_score(self, X, y):
        raise NotImplementedError("OOB score not supported by iforest")

    def fit(self, X, y=None, sample_weight=None):
        """Fit estimator.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Use ``dtype=np.float32`` for maximum
            efficiency. Sparse matrices are also supported, use sparse
            ``csc_matrix`` for maximum efficiency.

        Returns
        -------
        self : object
            Returns self.
        """
        # ensure_2d=False because there are actually unit test checking we fail
        # for 1d.
        X = check_array(X, accept_sparse=['csc'], ensure_2d=False)
        if issparse(X):
            # Pre-sort indices to avoid that each individual tree of the
            # ensemble sorts the indices.
            X.sort_indices()

        rnd = check_random_state(self.random_state)
        y = rnd.uniform(size=X.shape[0])

        # ensure that max_sample is in [1, n_samples]:
        n_samples = X.shape[0]

        if isinstance(self.max_samples, six.string_types):
            if self.max_samples == 'auto':
                max_samples = min(256, n_samples)
            else:
                raise ValueError('max_samples (%s) is not supported.'
                                 'Valid choices are: "auto", int or'
                                 'float' % self.max_samples)

        elif isinstance(self.max_samples, six.integer_types):
            if self.max_samples > n_samples:
                warn("max_samples (%s) is greater than the "
                     "total number of samples (%s). max_samples "
                     "will be set to n_samples for estimation."
                     % (self.max_samples, n_samples))
                max_samples = n_samples
            else:
                max_samples = self.max_samples
        else: # float
            if not (0. < self.max_samples <= 1.):
                raise ValueError("max_samples must be in (0, 1]")
            max_samples = int(self.max_samples * X.shape[0])

        self.max_samples_ = max_samples
        max_depth = int(np.ceil(np.log2(max(max_samples, 2))))
        super(IsolationForest, self)._fit(X, y, max_samples,
                                          max_depth=max_depth,
                                          sample_weight=sample_weight)
        return self

    def predict(self, X):
        """Predict anomaly score of X with the IsolationForest algorithm.

        The anomaly score of an input sample is computed as
        the mean anomaly score of the trees in the forest.

        The measure of normality of an observation given a tree is the depth
        of the leaf containing this observation, which is equivalent to
        the number of splittings required to isolate this point. In case of
        several observations n_left in the leaf, the average path length of
        a n_left samples isolation tree is added.

        Parameters
        ----------
        X : array-like or sparse matrix of shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        scores : array of shape (n_samples,)
            The anomaly score of the input samples.
            The lower, the more normal.
        """
        # code structure from ForestClassifier/predict_proba
        # Check data
        X = self.estimators_[0]._validate_X_predict(X, check_input=True)
        n_samples = X.shape[0]


        n_samples_leaf = np.zeros((n_samples, self.n_estimators), order="f")
        depths = np.zeros((n_samples, self.n_estimators), order="f")

        for i, tree in enumerate(self.estimators_):
            leaves_index = tree.apply(X)
            node_indicator = tree.decision_path(X)
            n_samples_leaf[:, i] = tree.tree_.n_node_samples[leaves_index]
            depths[:, i] = np.asarray(node_indicator.sum(axis=1)).reshape(-1) - 1

        depths += _average_path_length(n_samples_leaf)

        scores = 2 ** (-depths.mean(axis=1) / _average_path_length(self.max_samples_))

        return scores

    def decision_function(self, X):
        """Average of the decision functions of the base classifiers.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        Returns
        -------
        score : array, shape (n_samples,)
            The decision function of the input samples.

        """
        # minus as bigger is better (here less abnormal):
        return - self.predict(X)


def _average_path_length(n_samples_leaf):
    """ The average path length in a n_samples iTree, which is equal to
    the average path length of an unsuccessful BST search since the
    latter has the same structure as an isolation tree.
    Parameters
    ----------
    n_samples_leaf : array-like of shape (n_samples, n_estimators), or int.
        The number of training samples in each test sample leaf, for
        each estimators.
    
    Returns
    -------
    average_path_length : array, same shape as n_samples_leaf

    """
    if isinstance(n_samples_leaf, six.integer_types):
        if n_samples_leaf <= 1:
            return 1.
        else:
            return 2. * (np.log(n_samples_leaf) + 0.5772156649) - 2. * (
                n_samples_leaf - 1.) / n_samples_leaf

    else:

        n_samples_leaf_shape = n_samples_leaf.shape
        n_samples_leaf = n_samples_leaf.reshape((1, -1))
        average_path_length = np.zeros(n_samples_leaf.shape)

        mask = (n_samples_leaf <= 1)
        not_mask = np.logical_not(mask)

        average_path_length[mask] = 1.
        average_path_length[not_mask] = 2. * (
            np.log(n_samples_leaf[not_mask]) + 0.5772156649) - 2. * (
                n_samples_leaf[not_mask] - 1.) / n_samples_leaf[not_mask]

        return average_path_length.reshape(n_samples_leaf_shape)
