# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Vincent Michel <vincent.michel@inria.fr>
#          Gilles Louppe <g.louppe@gmail.com>
#
# License: BSD 3 clause

"""Recursive feature elimination for feature ranking"""

import numpy as np
from ..utils import check_X_y, safe_sqr
from ..base import BaseEstimator
from ..base import MetaEstimatorMixin
from ..base import clone
from ..base import is_classifier
from ..cross_validation import _check_cv as check_cv
from ..cross_validation import _safe_split, _score
from .base import SelectorMixin
from ..metrics.scorer import check_scoring


class RFE(BaseEstimator, MetaEstimatorMixin, SelectorMixin):
    """Feature ranking with recursive feature elimination.

    Given an external estimator that assigns weights to features (e.g., the
    coefficients of a linear model), the goal of recursive feature elimination
    (RFE) is to select features by recursively considering smaller and smaller
    sets of features. First, the estimator is trained on the initial set of
    features and weights are assigned to each one of them. Then, features whose
    absolute weights are the smallest are pruned from the current set features.
    That procedure is recursively repeated on the pruned set until the desired
    number of features to select is eventually reached.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a `fit` method that updates a
        `coef_` attribute that holds the fitted parameters. Important features
        must correspond to high absolute values in the `coef_` array.

        For instance, this is the case for most supervised learning
        algorithms such as Support Vector Classifiers and Generalized
        Linear Models from the `svm` and `linear_model` modules.

    n_features_to_select : int or None (default=None)
        The number of features to select. If `None`, half of the features
        are selected.

    step : int or float, optional (default=1)
        If greater than or equal to 1, then `step` corresponds to the (integer)
        number of features to remove at each iteration.
        If within (0.0, 1.0), then `step` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    estimator_params : dict
        Parameters for the external estimator.
        Useful for doing grid searches when an `RFE` object is passed as an
        argument to, e.g., a `sklearn.grid_search.GridSearchCV` object.

    Attributes
    ----------
    n_features_ : int
        The number of selected features.

    support_ : array of shape [n_features]
        The mask of selected features.

    ranking_ : array of shape [n_features]
        The feature ranking, such that `ranking_[i]` corresponds to the \
        ranking position of the i-th feature. Selected (i.e., estimated \
        best) features are assigned rank 1.

    estimator_ : object
        The external estimator fit on the reduced dataset.

    Examples
    --------
    The following example shows how to retrieve the 5 right informative
    features in the Friedman #1 dataset.

    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.feature_selection import RFE
    >>> from sklearn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFE(estimator, 5, step=1)
    >>> selector = selector.fit(X, y)
    >>> selector.support_ # doctest: +NORMALIZE_WHITESPACE
    array([ True,  True,  True,  True,  True,
            False, False, False, False, False], dtype=bool)
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

    References
    ----------

    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    def __init__(self, estimator, n_features_to_select=None, step=1,
                 estimator_params={}, verbose=0):
        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.step = step
        self.estimator_params = estimator_params
        self.verbose = verbose

    def fit(self, X, y):
        """Fit the RFE model and then the underlying estimator on the selected
           features.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values.
        """
        X, y = check_X_y(X, y, "csc")
        # Initialization
        n_features = X.shape[1]
        if self.n_features_to_select is None:
            n_features_to_select = n_features / 2
        else:
            n_features_to_select = self.n_features_to_select

        if 0.0 < self.step < 1.0:
            step = int(self.step * n_features)
        else:
            step = int(self.step)
        if step <= 0:
            raise ValueError("Step must be >0")

        support_ = np.ones(n_features, dtype=np.bool)
        ranking_ = np.ones(n_features, dtype=np.int)
        # Elimination
        while np.sum(support_) > n_features_to_select:
            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator = clone(self.estimator)
            estimator.set_params(**self.estimator_params)
            if self.verbose > 0:
                print("Fitting estimator with %d features." % np.sum(support_))

            estimator.fit(X[:, features], y)

            if estimator.coef_.ndim > 1:
                ranks = np.argsort(safe_sqr(estimator.coef_).sum(axis=0))
            else:
                ranks = np.argsort(safe_sqr(estimator.coef_))

            # for sparse case ranks is matrix
            ranks = np.ravel(ranks)

            # Eliminate the worse features
            threshold = min(step, np.sum(support_) - n_features_to_select)
            support_[features[ranks][:threshold]] = False
            ranking_[np.logical_not(support_)] += 1

        # Set final attributes
        self.estimator_ = clone(self.estimator)
        self.estimator_.set_params(**self.estimator_params)
        self.estimator_.fit(X[:, support_], y)
        self.n_features_ = support_.sum()
        self.support_ = support_
        self.ranking_ = ranking_

        return self

    def predict(self, X):
        """Reduce X to the selected features and then predict using the
           underlying estimator.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape [n_samples]
            The predicted target values.
        """
        return self.estimator_.predict(self.transform(X))

    def score(self, X, y):
        """Reduce X to the selected features and then return the score of the
           underlying estimator.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        y : array of shape [n_samples]
            The target values.
        """
        return self.estimator_.score(self.transform(X), y)

    def _get_support_mask(self):
        return self.support_

    def decision_function(self, X):
        return self.estimator_.decision_function(self.transform(X))

    def predict_proba(self, X):
        return self.estimator_.predict_proba(self.transform(X))


class RFECV(RFE, MetaEstimatorMixin):
    """Feature ranking with recursive feature elimination and cross-validated
       selection of the best number of features.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a `fit` method that updates a
        `coef_` attribute that holds the fitted parameters. Important features
        must correspond to high absolute values in the `coef_` array.

        For instance, this is the case for most supervised learning
        algorithms such as Support Vector Classifiers and Generalized
        Linear Models from the `svm` and `linear_model` modules.

    step : int or float, optional (default=1)
        If greater than or equal to 1, then `step` corresponds to the (integer)
        number of features to remove at each iteration.
        If within (0.0, 1.0), then `step` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    cv : int or cross-validation generator, optional (default=None)
        If int, it is the number of folds.
        If None, 3-fold cross-validation is performed by default.
        Specific cross-validation objects can also be passed, see
        `sklearn.cross_validation module` for details.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    estimator_params : dict
        Parameters for the external estimator.
        Useful for doing grid searches when an `RFE` object is passed as an
        argument to, e.g., a `sklearn.grid_search.GridSearchCV` object.

    verbose : int, default=0
        Controls verbosity of output.

    Attributes
    ----------
    n_features_ : int
        The number of selected features with cross-validation.
    support_ : array of shape [n_features]
        The mask of selected features.

    ranking_ : array of shape [n_features]
        The feature ranking, such that `ranking_[i]`
        corresponds to the ranking
        position of the i-th feature.
        Selected (i.e., estimated best)
        features are assigned rank 1.

    grid_scores_ : array of shape [n_subsets_of_features]
        The cross-validation scores such that
        `grid_scores_[i]` corresponds to
        the CV score of the i-th subset of features.

    estimator_ : object
        The external estimator fit on the reduced dataset.

    Examples
    --------
    The following example shows how to retrieve the a-priori not known 5
    informative features in the Friedman #1 dataset.

    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.feature_selection import RFECV
    >>> from sklearn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFECV(estimator, step=1, cv=5)
    >>> selector = selector.fit(X, y)
    >>> selector.support_ # doctest: +NORMALIZE_WHITESPACE
    array([ True,  True,  True,  True,  True,
            False, False, False, False, False], dtype=bool)
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

    References
    ----------

    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    def __init__(self, estimator, step=1, cv=None, scoring=None,
                 estimator_params={}, verbose=0):
        self.estimator = estimator
        self.step = step
        self.cv = cv
        self.scoring = scoring
        self.estimator_params = estimator_params
        self.verbose = verbose

    def fit(self, X, y):
        """Fit the RFE model and automatically tune the number of selected
           features.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the total number of features.

        y : array-like, shape = [n_samples]
            Target values (integers for classification, real numbers for
            regression).
        """
        X, y = check_X_y(X, y, "csr")
        # Initialization
        rfe = RFE(estimator=self.estimator, n_features_to_select=1,
                  step=self.step, estimator_params=self.estimator_params,
                  verbose=self.verbose - 1)

        cv = check_cv(self.cv, X, y, is_classifier(self.estimator))
        scorer = check_scoring(self.estimator, scoring=self.scoring)
        scores = np.zeros(X.shape[1])

        # Cross-validation
        for n, (train, test) in enumerate(cv):
            X_train, y_train = _safe_split(self.estimator, X, y, train)
            X_test, y_test = _safe_split(self.estimator, X, y, test, train)

            # Compute a full ranking of the features
            ranking_ = rfe.fit(X_train, y_train).ranking_
            # Score each subset of features
            for k in range(0, max(ranking_)):
                mask = np.where(ranking_ <= k + 1)[0]
                estimator = clone(self.estimator)
                estimator.fit(X_train[:, mask], y_train)
                score = _score(estimator, X_test[:, mask], y_test, scorer)

                if self.verbose > 0:
                    print("Finished fold with %d / %d feature ranks, score=%f"
                          % (k + 1, max(ranking_), score))
                scores[k] += score

        # Pick the best number of features on average
        k = np.argmax(scores)
        best_score = scores[k]

        # Re-execute an elimination with best_k over the whole set
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=k+1,
                  step=self.step, estimator_params=self.estimator_params)

        rfe.fit(X, y)

        # Set final attributes
        self.support_ = rfe.support_
        self.n_features_ = rfe.n_features_
        self.ranking_ = rfe.ranking_
        self.estimator_ = clone(self.estimator)
        self.estimator_.set_params(**self.estimator_params)
        self.estimator_.fit(self.transform(X), y)

        # Fixing a normalization error, n is equal to len(cv) - 1
        # here, the scores are normalized by len(cv)
        self.grid_scores_ = scores / len(cv)
        return self
