# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Vincent Michel <vincent.michel@inria.fr>
#          Gilles Louppe <g.louppe@gmail.com>
#
# License: BSD Style.

"""Recursive feature elimination for feature ranking"""

import numpy as np
from ..base import BaseEstimator
from ..base import clone
from ..base import is_classifier
from ..cross_validation import check_cv


class RFE(BaseEstimator):
    """Feature ranking with recursive feature elimination.

    Given an external estimator that assigns weights to features (e.g., the
    coefficients of a linear model), the goal of recursive feature elimination
    (RFE) is to select features by recursively considering smaller and smaller
    sets of features.  First, the estimator is trained on the initial set of
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

    n_features_to_select : int
        The number of features to select.

    step : int or float, optional (default=1)
        If greater than or equal to 1, then `step` corresponds to the (integer)
        number of features to remove at each iteration.
        If within (0.0, 1.0), then `step` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    Attributes
    ----------
    `n_features_` : int
        The number of selected features.

    `support_` : array of shape [n_features]
        The mask of selected features.

    `ranking_` : array of shape [n_features]
        The feature ranking, such that `ranking_[i]` corresponds to the ranking \
        position of the i-th feature. Selected (i.e., estimated best) features \
        are assigned rank 1.

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

    Notes
    -----
    **References**:

    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    def __init__(self, estimator, n_features_to_select, step=1):
        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.step = step

    def fit(self, X, y):
        """Fit the RFE model and then the underlying estimator on the selected
           features.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The training input samples.

        y : array of shape [n_samples]
            The target values.
        """
        # Initialization
        n_features = X.shape[1]

        if 0.0 < self.step < 1.0:
            step = int(self.step * n_features)
        else:
            step = int(self.step)
        assert step > 0

        support_ = np.ones(n_features, dtype=np.bool)
        ranking_ = np.ones(n_features, dtype=np.int)

        # Elimination
        while np.sum(support_) > self.n_features_to_select:
            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator = clone(self.estimator)
            estimator.fit(X[:, features], y)

            if estimator.coef_.ndim > 1:
                ranks = np.argsort(np.sum(estimator.coef_ ** 2, axis=0))
            else:
                ranks = np.argsort(estimator.coef_ ** 2)

            # Eliminate the worse features
            threshold = min(step, np.sum(support_) - self.n_features_to_select)
            support_[features[ranks][:threshold]] = False
            ranking_[np.logical_not(support_)] += 1

        # Set final attributes
        self.estimator.fit(X[:, support_], y)
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
        return self.estimator.predict(X[:, self.support_])

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
        return self.estimator.score(X[:, self.support_], y)

    def transform(self, X):
        """Reduce X to the selected features during the elimination.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the features selected during the \
            elimination.
        """
        return X[:, self.support_]


class RFECV(RFE):
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
        `scikits.learn.cross_validation module` for details.

    loss_function : function, optional (default=None)
        The loss function to minimize by cross-validation. If None, then the
        score function of the estimator is maximized.

    Attributes
    ----------
    `n_features_` : int
        The number of selected features with cross-validation.
    `support_` : array of shape [n_features]
        The mask of selected features.

    `ranking_` : array of shape [n_features]
        The feature ranking, such that `ranking_[i]`
        corresponds to the ranking
        position of the i-th feature.
        Selected (i.e., estimated best)
        features are assigned rank 1.

    `cv_scores_` : array of shape [n_subsets_of_features]
        The cross-validation scores such that
        `cv_scores_[i]` corresponds to
        the CV score of the i-th subset of features.

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

    Notes
    -----
    **References**:

    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    def __init__(self, estimator, step=1, cv=None, loss_func=None):
        self.estimator = estimator
        self.step = 1
        self.cv = cv
        self.loss_func = loss_func

    def fit(self, X, y):
        """Fit the RFE model and automatically tune the number of selected
           features.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the total number of features.

        y : array of shape [n_samples]
            Target values (integers for classification, real numbers for
            regression).
        """
        # Initialization
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=1,
                  step=self.step)

        cv = check_cv(self.cv, X, y, is_classifier(self.estimator))
        scores = {}

        # Cross-validation
        n = 0

        for train, test in cv:
            # Compute a full ranking of the features
            ranking_ = rfe.fit(X[train], y[train]).ranking_

            # Score each subset of features
            for k in xrange(1, max(ranking_)):
                mask = ranking_ <= k
                estimator = clone(self.estimator)
                estimator.fit(X[train][:, mask], y[train])

                if self.loss_func is None:
                    score_k = 1.0 - estimator.score(
                                  X[test][:, mask],
                                  y[test])
                else:
                    score_k = self.loss_func(
                                  y[test],
                                  estimator.predict(X[test][:, mask]))

                if not k in scores:
                    scores[k] = 0.0

                scores[k] += score_k

            n += 1

        # Pick the best number of features on average
        best_score = np.inf
        best_k = None

        for k, score in sorted(scores.iteritems()):
            if score < best_score:
                best_score = score
                best_k = k

        # Re-execute an elimination with best_k over the whole set
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=best_k,
                  step=self.step)

        rfe.fit(X, y)

        # Set final attributes
        self.estimator.fit(X[:, rfe.support_], y)
        self.n_features_ = rfe.n_features_
        self.support_ = rfe.support_
        self.ranking_ = rfe.ranking_

        self.cv_scores_ = [0] * len(scores)
        for k, score in scores.iteritems():
            self.cv_scores_[k - 1] = score / n

        return self
