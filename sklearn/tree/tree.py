"""
This module gathers tree-based methods, including decision, regression and
randomized trees. Single and multi-output problems are both handled.
"""

# Code is originally adapted from MILK: Machine Learning Toolkit
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# License: MIT. See COPYING.MIT file in the milk distribution

# Authors: Brian Holt,
#          Peter Prettenhofer,
#          Satrajit Ghosh,
#          Gilles Louppe,
#          Noel Dawe
# License: BSD Style.

from __future__ import division

import numbers
import numpy as np
from abc import ABCMeta, abstractmethod
from warnings import warn

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..externals import six
from ..externals.six.moves import xrange
from ..feature_selection.selector_mixin import SelectorMixin
from ..utils import array2d, check_random_state
from ..utils.validation import check_arrays

from . import _tree


__all__ = ["DecisionTreeClassifier",
           "DecisionTreeRegressor",
           "ExtraTreeClassifier",
           "ExtraTreeRegressor"]

DTYPE = _tree.DTYPE
DOUBLE = _tree.DOUBLE

CLASSIFICATION = {
    "gini": _tree.Gini,
    "entropy": _tree.Entropy,
}

REGRESSION = {
    "mse": _tree.MSE,
}


class BaseDecisionTree(six.with_metaclass(ABCMeta, BaseEstimator, SelectorMixin)):
    """Base class for decision trees.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """

    @abstractmethod
    def __init__(self,
                 criterion,
                 max_depth,
                 min_samples_split,
                 min_samples_leaf,
                 min_density,
                 max_features,
                 compute_importances,
                 random_state):
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features

        if compute_importances:
            warn("Setting compute_importances=True is no longer "
                 "required. Variable importances are now computed on the fly "
                 "when accessing the feature_importances_ attribute. This "
                 "parameter will be removed in 0.15.", DeprecationWarning)

        self.compute_importances = compute_importances
        self.random_state = random_state

        self.n_features_ = None
        self.n_outputs_ = None
        self.classes_ = None
        self.n_classes_ = None
        self.find_split_ = _tree.TREE_SPLIT_BEST

        self.tree_ = None

    def fit(self, X, y,
            sample_mask=None, X_argsorted=None,
            check_input=True, sample_weight=None):
        """Build a decision tree from the training set (X, y).

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples. Use ``dtype=np.float32``
            and ``order='F'`` for maximum efficiency.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).
            Use ``dtype=np.float64`` and ``order='C'`` for maximum
            efficiency.

        sample_mask : array-like, shape = [n_samples], dtype = bool or None
            A bit mask that encodes the rows of ``X`` that should be
            used to build the decision tree. It can be used for bagging
            without the need to create of copy of ``X``.
            If None a mask will be created that includes all samples.

        X_argsorted : array-like, shape = [n_samples, n_features] or None
            Each column of ``X_argsorted`` holds the row indices of ``X``
            sorted according to the value of the corresponding feature
            in ascending order.
            I.e. ``X[X_argsorted[i, k], k] <= X[X_argsorted[j, k], k]``
            for each j > i.
            If None, ``X_argsorted`` is computed internally.
            The argument is supported to enable multiple decision trees
            to share the data structure and to avoid re-computation in
            tree ensembles. For maximum efficiency use dtype np.int32.

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. In the case of
            classification, splits are also ignored if they would result in any
            single class carrying a negative weight in either child node.

        check_input : boolean, (default=True)
            Allow to bypass several input checking.
            Don't use this parameter unless you know what you do.

        Returns
        -------
        self : object
            Returns self.
        """
        if check_input:
            X, y = check_arrays(X, y)
        random_state = check_random_state(self.random_state)

        # Convert data
        if (getattr(X, "dtype", None) != DTYPE or
                X.ndim != 2 or
                not X.flags.fortran):
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, self.n_features_ = X.shape
        is_classification = isinstance(self, ClassifierMixin)

        y = np.atleast_1d(y)
        if y.ndim == 1:
            # reshape is necessary to preserve the data contiguity against vs
            # [:, np.newaxis] that does not.
            y = np.reshape(y, (-1, 1))

        self.n_outputs_ = y.shape[1]

        if is_classification:
            y = np.copy(y)

            self.classes_ = []
            self.n_classes_ = []

            for k in xrange(self.n_outputs_):
                unique = np.unique(y[:, k])
                self.classes_.append(unique)
                self.n_classes_.append(unique.shape[0])
                y[:, k] = np.searchsorted(unique, y[:, k])

        else:
            self.classes_ = [None] * self.n_outputs_
            self.n_classes_ = [1] * self.n_outputs_

        if getattr(y, "dtype", None) != DOUBLE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        if is_classification:
            criterion = CLASSIFICATION[self.criterion](self.n_outputs_,
                                                       self.n_classes_)
        else:
            criterion = REGRESSION[self.criterion](self.n_outputs_)

        # Check parameters
        max_depth = np.inf if self.max_depth is None else self.max_depth

        if isinstance(self.max_features, six.string_types):
            if self.max_features == "auto":
                if is_classification:
                    max_features = max(1, int(np.sqrt(self.n_features_)))
                else:
                    max_features = self.n_features_
            elif self.max_features == "sqrt":
                max_features = max(1, int(np.sqrt(self.n_features_)))
            elif self.max_features == "log2":
                max_features = max(1, int(np.log2(self.n_features_)))
            else:
                raise ValueError(
                    'Invalid value for max_features. Allowed string '
                    'values are "auto", "sqrt" or "log2".')
        elif self.max_features is None:
            max_features = self.n_features_
        elif isinstance(self.max_features, (numbers.Integral, np.integer)):
            max_features = self.max_features
        else:  # float
            max_features = int(self.max_features * self.n_features_)

        if len(y) != n_samples:
            raise ValueError("Number of labels=%d does not match "
                             "number of samples=%d" % (len(y), n_samples))
        if self.min_samples_split <= 0:
            raise ValueError("min_samples_split must be greater than zero.")
        if self.min_samples_leaf <= 0:
            raise ValueError("min_samples_leaf must be greater than zero.")
        if max_depth <= 0:
            raise ValueError("max_depth must be greater than zero. ")
        if self.min_density < 0.0 or self.min_density > 1.0:
            raise ValueError("min_density must be in [0, 1]")
        if not (0 < max_features <= self.n_features_):
            raise ValueError("max_features must be in (0, n_features]")

        if sample_mask is not None:
            sample_mask = np.asarray(sample_mask, dtype=np.bool)

            if sample_mask.shape[0] != n_samples:
                raise ValueError("Length of sample_mask=%d does not match "
                                 "number of samples=%d"
                                 % (sample_mask.shape[0], n_samples))

        if sample_weight is not None:
            if (getattr(sample_weight, "dtype", None) != DOUBLE or
                    not sample_weight.flags.contiguous):
                sample_weight = np.ascontiguousarray(
                    sample_weight, dtype=DOUBLE)
            if len(sample_weight.shape) > 1:
                raise ValueError("Sample weights array has more "
                                 "than one dimension: %d" %
                                 len(sample_weight.shape))
            if len(sample_weight) != n_samples:
                raise ValueError("Number of weights=%d does not match "
                                 "number of samples=%d" %
                                 (len(sample_weight), n_samples))

        if X_argsorted is not None:
            X_argsorted = np.asarray(X_argsorted, dtype=np.int32,
                                     order='F')
            if X_argsorted.shape != X.shape:
                raise ValueError("Shape of X_argsorted does not match "
                                 "the shape of X")

        # Set min_samples_split sensibly
        min_samples_split = max(self.min_samples_split,
                                2 * self.min_samples_leaf)

        # Build tree
        self.tree_ = _tree.Tree(self.n_features_, self.n_classes_,
                                self.n_outputs_, criterion, max_depth,
                                min_samples_split, self.min_samples_leaf,
                                self.min_density, max_features,
                                self.find_split_, random_state)

        self.tree_.build(X, y,
                         sample_weight=sample_weight,
                         sample_mask=sample_mask,
                         X_argsorted=X_argsorted)

        if self.n_outputs_ == 1:
            self.n_classes_ = self.n_classes_[0]
            self.classes_ = self.classes_[0]

        return self

    def predict(self, X):
        """Predict class or regression value for X.

        For a classification model, the predicted class for each sample in X is
        returned. For a regression model, the predicted value based on X is
        returned.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted classes, or the predict values.
        """
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, n_features = X.shape

        if self.tree_ is None:
            raise Exception("Tree not initialized. Perform a fit first")

        if self.n_features_ != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features_, n_features))

        proba = self.tree_.predict(X)

        # Classification
        if isinstance(self, ClassifierMixin):
            if self.n_outputs_ == 1:
                return self.classes_.take(np.argmax(proba[:, 0], axis=1),
                                          axis=0)

            else:
                predictions = np.zeros((n_samples, self.n_outputs_))

                for k in xrange(self.n_outputs_):
                    predictions[:, k] = self.classes_[k].take(
                        np.argmax(proba[:, k], axis=1),
                        axis=0)

                return predictions

        # Regression
        else:
            if self.n_outputs_ == 1:
                return proba[:, 0, 0]

            else:
                return proba[:, :, 0]

    @property
    def feature_importances_(self):
        """Return the feature importances.

        The importance of a feature is computed as the (normalized) total
        reduction of the criterion brought by that feature.
        It is also known as the Gini importance [4]_.

        Returns
        -------
        feature_importances_ : array, shape = [n_features]
        """
        if self.tree_ is None:
            raise ValueError("Estimator not fitted, "
                             "call `fit` before `feature_importances_`.")

        return self.tree_.compute_feature_importances()


class DecisionTreeClassifier(BaseDecisionTree, ClassifierMixin):
    """A decision tree classifier.

    Parameters
    ----------
    criterion : string, optional (default="gini")
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "entropy" for the information gain.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:
          - If int, then consider `max_features` features at each split.
          - If float, then `max_features` is a percentage and
            `int(max_features * n_features)` features are considered at each
            split.
          - If "auto", then `max_features=sqrt(n_features)`.
          - If "sqrt", then `max_features=sqrt(n_features)`.
          - If "log2", then `max_features=log2(n_features)`.
          - If None, then `max_features=n_features`.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : integer, optional (default=2)
        The minimum number of samples required to split an internal node.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples required to be at a leaf node.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `tree_` : Tree object
        The underlying Tree object.

    `classes_` : array of shape = [n_classes] or a list of such arrays
        The classes labels (single output problem),
        or a list of arrays of class labels (multi-output problem).

    `n_classes_` : int or list
        The number of classes (for single output problems),
        or a list containing the number of classes for each
        output (for multi-output problems).

    `feature_importances_` : array of shape = [n_features]
        The feature importances. The higher, the more important the
        feature. The importance of a feature is computed as the (normalized)
        total reduction of the criterion brought by that feature.  It is also
        known as the Gini importance [4]_.

    See also
    --------
    DecisionTreeRegressor

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    .. [4] L. Breiman, and A. Cutler, "Random Forests",
           http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.cross_validation import cross_val_score
    >>> from sklearn.tree import DecisionTreeClassifier

    >>> clf = DecisionTreeClassifier(random_state=0)
    >>> iris = load_iris()

    >>> cross_val_score(clf, iris.data, iris.target, cv=10)
    ...                             # doctest: +SKIP
    ...
    array([ 1.     ,  0.93...,  0.86...,  0.93...,  0.93...,
            0.93...,  0.93...,  1.     ,  0.93...,  1.      ])
    """
    def __init__(self,
                 criterion="gini",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_density=0.1,
                 max_features=None,
                 compute_importances=False,
                 random_state=None):
        super(DecisionTreeClassifier, self).__init__(criterion,
                                                     max_depth,
                                                     min_samples_split,
                                                     min_samples_leaf,
                                                     min_density,
                                                     max_features,
                                                     compute_importances,
                                                     random_state)

    def predict_proba(self, X):
        """Predict class probabilities of the input samples X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are ordered
            by arithmetical order.
        """
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, n_features = X.shape

        if self.tree_ is None:
            raise Exception("Tree not initialized. Perform a fit first.")

        if self.n_features_ != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features_, n_features))

        proba = self.tree_.predict(X)

        if self.n_outputs_ == 1:
            proba = proba[:, 0, :self.n_classes_]
            normalizer = proba.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            proba /= normalizer

            return proba

        else:
            all_proba = []

            for k in xrange(self.n_outputs_):
                proba_k = proba[:, k, :self.n_classes_[k]]
                normalizer = proba_k.sum(axis=1)[:, np.newaxis]
                normalizer[normalizer == 0.0] = 1.0
                proba_k /= normalizer
                all_proba.append(proba_k)

            return all_proba

    def predict_log_proba(self, X):
        """Predict class log-probabilities of the input samples X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            such arrays if n_outputs > 1.
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        proba = self.predict_proba(X)

        if self.n_outputs_ == 1:
            return np.log(proba)

        else:
            for k in xrange(self.n_outputs_):
                proba[k] = np.log(proba[k])

            return proba


class DecisionTreeRegressor(BaseDecisionTree, RegressorMixin):
    """A tree regressor.

    Parameters
    ----------
    criterion : string, optional (default="mse")
        The function to measure the quality of a split. The only supported
        criterion is "mse" for the mean squared error.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:
          - If int, then consider `max_features` features at each split.
          - If float, then `max_features` is a percentage and
            `int(max_features * n_features)` features are considered at each
            split.
          - If "auto", then `max_features=n_features`.
          - If "sqrt", then `max_features=sqrt(n_features)`.
          - If "log2", then `max_features=log2(n_features)`.
          - If None, then `max_features=n_features`.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : integer, optional (default=2)
        The minimum number of samples required to split an internal node.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples required to be at a leaf node.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `tree_` : Tree object
        The underlying Tree object.

    `feature_importances_` : array of shape = [n_features]
        The feature importances.
        The higher, the more important the feature.
        The importance of a feature is computed as the
        (normalized)total reduction of the criterion brought
        by that feature. It is also known as the Gini importance [4]_.

    See also
    --------
    DecisionTreeClassifier

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    .. [4] L. Breiman, and A. Cutler, "Random Forests",
           http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm

    Examples
    --------
    >>> from sklearn.datasets import load_boston
    >>> from sklearn.cross_validation import cross_val_score
    >>> from sklearn.tree import DecisionTreeRegressor

    >>> boston = load_boston()
    >>> regressor = DecisionTreeRegressor(random_state=0)

    R2 scores (a.k.a. coefficient of determination) over 10-folds CV:

    >>> cross_val_score(regressor, boston.data, boston.target, cv=10)
    ...                    # doctest: +SKIP
    ...
    array([ 0.61..., 0.57..., -0.34..., 0.41..., 0.75...,
            0.07..., 0.29..., 0.33..., -1.42..., -1.77...])
    """
    def __init__(self,
                 criterion="mse",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_density=0.1,
                 max_features=None,
                 compute_importances=False,
                 random_state=None):
        super(DecisionTreeRegressor, self).__init__(criterion,
                                                    max_depth,
                                                    min_samples_split,
                                                    min_samples_leaf,
                                                    min_density,
                                                    max_features,
                                                    compute_importances,
                                                    random_state)


class ExtraTreeClassifier(DecisionTreeClassifier):
    """An extremely randomized tree classifier.

    Extra-trees differ from classic decision trees in the way they are built.
    When looking for the best split to separate the samples of a node into two
    groups, random splits are drawn for each of the `max_features` randomly
    selected features and the best split among those is chosen. When
    `max_features` is set 1, this amounts to building a totally random
    decision tree.

    Warning: Extra-trees should only be used within ensemble methods.

    See also
    --------
    ExtraTreeRegressor, ExtraTreesClassifier, ExtraTreesRegressor

    References
    ----------

    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self,
                 criterion="gini",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_density=0.1,
                 max_features="auto",
                 compute_importances=False,
                 random_state=None):
        super(ExtraTreeClassifier, self).__init__(criterion,
                                                  max_depth,
                                                  min_samples_split,
                                                  min_samples_leaf,
                                                  min_density,
                                                  max_features,
                                                  compute_importances,
                                                  random_state)

        self.find_split_ = _tree.TREE_SPLIT_RANDOM


class ExtraTreeRegressor(DecisionTreeRegressor):
    """An extremely randomized tree regressor.

    Extra-trees differ from classic decision trees in the way they are built.
    When looking for the best split to separate the samples of a node into two
    groups, random splits are drawn for each of the `max_features` randomly
    selected features and the best split among those is chosen. When
    `max_features` is set 1, this amounts to building a totally random
    decision tree.

    Warning: Extra-trees should only be used within ensemble methods.

    See also
    --------
    ExtraTreeClassifier : A classifier base on extremely randomized trees
    sklearn.ensemble.ExtraTreesClassifier : An ensemble of extra-trees for
        classification
    sklearn.ensemble.ExtraTreesRegressor : An ensemble of extra-trees for
        regression

    References
    ----------

    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self,
                 criterion="mse",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_density=0.1,
                 max_features="auto",
                 compute_importances=False,
                 random_state=None):
        super(ExtraTreeRegressor, self).__init__(criterion,
                                                 max_depth,
                                                 min_samples_split,
                                                 min_samples_leaf,
                                                 min_density,
                                                 max_features,
                                                 compute_importances,
                                                 random_state)

        self.find_split_ = _tree.TREE_SPLIT_RANDOM
