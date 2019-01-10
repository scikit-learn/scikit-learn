# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

from __future__ import division

import numpy as np
from warnings import warn
from sklearn.utils.fixes import euler_gamma

from scipy.sparse import issparse

import numbers
from ..tree import ExtraTreeRegressor
from ..utils import check_random_state, check_array
from ..utils.fixes import _joblib_parallel_args
from ..utils.validation import check_is_fitted
from ..base import OutlierMixin

from .bagging import BaseBagging

__all__ = ["IsolationForest"]

INTEGER_TYPES = (numbers.Integral, np.integer)


class IsolationForest(BaseBagging, OutlierMixin):
    """Isolation Forest Algorithm

    Return the anomaly score of each sample using the IsolationForest algorithm

    The IsolationForest 'isolates' observations by randomly selecting a feature
    and then randomly selecting a split value between the maximum and minimum
    values of the selected feature.

    Since recursive partitioning can be represented by a tree structure, the
    number of splittings required to isolate a sample is equivalent to the path
    length from the root node to the terminating node.

    This path length, averaged over a forest of such random trees, is a
    measure of normality and our decision function.

    Random partitioning produces noticeably shorter paths for anomalies.
    Hence, when a forest of random trees collectively produce shorter path
    lengths for particular samples, they are highly likely to be anomalies.

    Read more in the :ref:`User Guide <isolation_forest>`.

    .. versionadded:: 0.18

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

    contamination : float in (0., 0.5), optional (default=0.1)
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set. Used when fitting to define the threshold
        on the decision function. If 'auto', the decision function threshold is
        determined as in the original paper.

        .. versionchanged:: 0.20
           The default value of ``contamination`` will change from 0.1 in 0.20
           to ``'auto'`` in 0.22.

    max_features : int or float, optional (default=1.0)
        The number of features to draw from X to train each base estimator.

            - If int, then draw `max_features` features.
            - If float, then draw `max_features * X.shape[1]` features.

    bootstrap : boolean, optional (default=False)
        If True, individual trees are fit on random subsets of the training
        data sampled with replacement. If False, sampling without replacement
        is performed.

    n_jobs : int or None, optional (default=None)
        The number of jobs to run in parallel for both `fit` and `predict`.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    behaviour : str, default='old'
        Behaviour of the ``decision_function`` which can be either 'old' or
        'new'. Passing ``behaviour='new'`` makes the ``decision_function``
        change to match other anomaly detection algorithm API which will be
        the default behaviour in the future. As explained in details in the
        ``offset_`` attribute documentation, the ``decision_function`` becomes
        dependent on the contamination parameter, in such a way that 0 becomes
        its natural threshold to detect outliers.

        .. versionadded:: 0.20
           ``behaviour`` is added in 0.20 for back-compatibility purpose.

        .. deprecated:: 0.20
           ``behaviour='old'`` is deprecated in 0.20 and will not be possible
           in 0.22.

        .. deprecated:: 0.22
           ``behaviour`` parameter will be deprecated in 0.22 and removed in
           0.24.

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

    offset_ : float
        Offset used to define the decision function from the raw scores.
        We have the relation: ``decision_function = score_samples - offset_``.
        Assuming behaviour == 'new', ``offset_`` is defined as follows.
        When the contamination parameter is set to "auto", the offset is equal
        to -0.5 as the scores of inliers are close to 0 and the scores of
        outliers are close to -1. When a contamination parameter different
        than "auto" is provided, the offset is defined in such a way we obtain
        the expected number of outliers (samples with decision function < 0)
        in training.
        Assuming the behaviour parameter is set to 'old', we always have
        ``offset_ = -0.5``, making the decision function independent from the
        contamination parameter.

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
                 contamination="legacy",
                 max_features=1.,
                 bootstrap=False,
                 n_jobs=None,
                 behaviour='old',
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

        self.behaviour = behaviour
        self.contamination = contamination

    def _set_oob_score(self, X, y):
        raise NotImplementedError("OOB score not supported by iforest")

    def _parallel_args(self):
        # ExtraTreeRegressor releases the GIL, so it's more efficient to use
        # a thread-based backend rather than a process-based backend so as
        # to avoid suffering from communication overhead and extra memory
        # copies.
        return _joblib_parallel_args(prefer='threads')

    def fit(self, X, y=None, sample_weight=None):
        """Fit estimator.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Use ``dtype=np.float32`` for maximum
            efficiency. Sparse matrices are also supported, use sparse
            ``csc_matrix`` for maximum efficiency.

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted.

        y : Ignored
            not used, present for API consistency by convention.

        Returns
        -------
        self : object
        """
        if self.contamination == "legacy":
            warn('default contamination parameter 0.1 will change '
                 'in version 0.22 to "auto". This will change the '
                 'predict method behavior.',
                 FutureWarning)
            self._contamination = 0.1
        else:
            self._contamination = self.contamination

        if self.behaviour == 'old':
            warn('behaviour="old" is deprecated and will be removed '
                 'in version 0.22. Please use behaviour="new", which '
                 'makes the decision_function change to match '
                 'other anomaly detection algorithm API.',
                 FutureWarning)

        X = check_array(X, accept_sparse=['csc'])
        if issparse(X):
            # Pre-sort indices to avoid that each individual tree of the
            # ensemble sorts the indices.
            X.sort_indices()

        rnd = check_random_state(self.random_state)
        y = rnd.uniform(size=X.shape[0])

        # ensure that max_sample is in [1, n_samples]:
        n_samples = X.shape[0]

        if isinstance(self.max_samples, str):
            if self.max_samples == 'auto':
                max_samples = min(256, n_samples)
            else:
                raise ValueError('max_samples (%s) is not supported.'
                                 'Valid choices are: "auto", int or'
                                 'float' % self.max_samples)

        elif isinstance(self.max_samples, INTEGER_TYPES):
            if self.max_samples > n_samples:
                warn("max_samples (%s) is greater than the "
                     "total number of samples (%s). max_samples "
                     "will be set to n_samples for estimation."
                     % (self.max_samples, n_samples))
                max_samples = n_samples
            else:
                max_samples = self.max_samples
        else:  # float
            if not (0. < self.max_samples <= 1.):
                raise ValueError("max_samples must be in (0, 1], got %r"
                                 % self.max_samples)
            max_samples = int(self.max_samples * X.shape[0])

        self.max_samples_ = max_samples
        max_depth = int(np.ceil(np.log2(max(max_samples, 2))))
        super(IsolationForest, self)._fit(X, y, max_samples,
                                          max_depth=max_depth,
                                          sample_weight=sample_weight)

        if self.behaviour == 'old':
            # in this case, decision_function = 0.5 + self.score_samples(X):
            if self._contamination == "auto":
                raise ValueError("contamination parameter cannot be set to "
                                 "'auto' when behaviour == 'old'.")

            self.offset_ = -0.5
            self._threshold_ = np.percentile(self.decision_function(X),
                                             100. * self._contamination)

            return self

        # else, self.behaviour == 'new':
        if self._contamination == "auto":
            # 0.5 plays a special role as described in the original paper.
            # we take the opposite as we consider the opposite of their score.
            self.offset_ = -0.5
            return self

        # else, define offset_ wrt contamination parameter, so that the
        # threshold_ attribute is implicitly 0 and is not needed anymore:
        self.offset_ = np.percentile(self.score_samples(X),
                                     100. * self._contamination)

        return self

    def predict(self, X):
        """Predict if a particular sample is an outlier or not.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        is_inlier : array, shape (n_samples,)
            For each observation, tells whether or not (+1 or -1) it should
            be considered as an inlier according to the fitted model.
        """
        check_is_fitted(self, ["offset_"])
        X = check_array(X, accept_sparse='csr')
        is_inlier = np.ones(X.shape[0], dtype=int)
        threshold = self.threshold_ if self.behaviour == 'old' else 0
        is_inlier[self.decision_function(X) < threshold] = -1
        return is_inlier

    def decision_function(self, X):
        """Average anomaly score of X of the base classifiers.

        The anomaly score of an input sample is computed as
        the mean anomaly score of the trees in the forest.

        The measure of normality of an observation given a tree is the depth
        of the leaf containing this observation, which is equivalent to
        the number of splittings required to isolate this point. In case of
        several observations n_left in the leaf, the average path length of
        a n_left samples isolation tree is added.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        Returns
        -------
        scores : array, shape (n_samples,)
            The anomaly score of the input samples.
            The lower, the more abnormal. Negative scores represent outliers,
            positive scores represent inliers.

        """
        # We subtract self.offset_ to make 0 be the threshold value for being
        # an outlier:

        return self.score_samples(X) - self.offset_

    def score_samples(self, X):
        """Opposite of the anomaly score defined in the original paper.

        The anomaly score of an input sample is computed as
        the mean anomaly score of the trees in the forest.

        The measure of normality of an observation given a tree is the depth
        of the leaf containing this observation, which is equivalent to
        the number of splittings required to isolate this point. In case of
        several observations n_left in the leaf, the average path length of
        a n_left samples isolation tree is added.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        Returns
        -------
        scores : array, shape (n_samples,)
            The anomaly score of the input samples.
            The lower, the more abnormal.
        """
        # code structure from ForestClassifier/predict_proba
        check_is_fitted(self, ["estimators_"])

        # Check data
        X = check_array(X, accept_sparse='csr')
        if self.n_features_ != X.shape[1]:
            raise ValueError("Number of features of the model must "
                             "match the input. Model n_features is {0} and "
                             "input n_features is {1}."
                             "".format(self.n_features_, X.shape[1]))
        n_samples = X.shape[0]

        n_samples_leaf = np.zeros((n_samples, self.n_estimators), order="f")
        depths = np.zeros((n_samples, self.n_estimators), order="f")

        if self._max_features == X.shape[1]:
            subsample_features = False
        else:
            subsample_features = True

        for i, (tree, features) in enumerate(zip(self.estimators_,
                                                 self.estimators_features_)):
            if subsample_features:
                X_subset = X[:, features]
            else:
                X_subset = X
            leaves_index = tree.apply(X_subset)
            node_indicator = tree.decision_path(X_subset)
            n_samples_leaf[:, i] = tree.tree_.n_node_samples[leaves_index]
            depths[:, i] = np.ravel(node_indicator.sum(axis=1))
            depths[:, i] -= 1

        depths += _average_path_length(n_samples_leaf)

        scores = 2 ** (-depths.mean(axis=1) / _average_path_length(
            self.max_samples_))

        # Take the opposite of the scores as bigger is better (here less
        # abnormal)
        return -scores

    @property
    def threshold_(self):
        if self.behaviour != 'old':
            raise AttributeError("threshold_ attribute does not exist when "
                                 "behaviour != 'old'")
        warn("threshold_ attribute is deprecated in 0.20 and will"
             " be removed in 0.22.", DeprecationWarning)
        return self._threshold_


def _average_path_length(n_samples_leaf):
    """ The average path length in a n_samples iTree, which is equal to
    the average path length of an unsuccessful BST search since the
    latter has the same structure as an isolation tree.
    Parameters
    ----------
    n_samples_leaf : array-like, shape (n_samples, n_estimators), or int.
        The number of training samples in each test sample leaf, for
        each estimators.

    Returns
    -------
    average_path_length : array, same shape as n_samples_leaf

    """
    if isinstance(n_samples_leaf, INTEGER_TYPES):
        if n_samples_leaf <= 1:
            return 1.
        else:
            return 2. * (np.log(n_samples_leaf - 1.) + euler_gamma) - 2. * (
                n_samples_leaf - 1.) / n_samples_leaf

    else:

        n_samples_leaf_shape = n_samples_leaf.shape
        n_samples_leaf = n_samples_leaf.reshape((1, -1))
        average_path_length = np.zeros(n_samples_leaf.shape)

        mask = (n_samples_leaf <= 1)
        not_mask = np.logical_not(mask)

        average_path_length[mask] = 1.
        average_path_length[not_mask] = 2. * (
            np.log(n_samples_leaf[not_mask] - 1.) + euler_gamma) - 2. * (
                n_samples_leaf[not_mask] - 1.) / n_samples_leaf[not_mask]

        return average_path_length.reshape(n_samples_leaf_shape)
