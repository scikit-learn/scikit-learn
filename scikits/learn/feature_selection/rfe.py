# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Vincent Michel <vincent.michel@inria.fr>
#          Gilles Louppe <g.louppe@gmail.com>
#
# License: BSD Style.

"""Recursive feature elimination for feature ranking"""

import numpy as np
from ..base import BaseEstimator
from ..base import is_classifier
from ..cross_val import check_cv


class RFE(BaseEstimator):
    """Feature ranking with recursive feature elimination.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a `fit` method that updates a
        `coef_` attribute that holds the fitted parameters. The first
        dimension of the `coef_` array must be equal to the number of features
        of the input dataset of the estimator. Important features must
        correspond to high absolute values in the `coef_` array.

        For instance, this is the case for most supervised learning
        algorithms such as Support Vector Classifiers and Generalized
        Linear Models from the `svm` and `linear_model` modules.

    n_features_to_select : int
        The number of features to select.

    step : int or float, optional (default=1)
        If int, then `step` corresponds to the number of features to remove
        at each iteration.
        If float, then `step` corresponds to the percentage (rounded down) of
        features to remove at each iteration. This value should be in (0, 1].

    Attributes
    ----------
    `support_` : array of shape [n_features]
        The mask of selected features.

    `ranking_` : array of shape [n_features]
        The feature ranking, such that `ranking_[i]` corresponds to the ranking
        position of the i-th feature. Selected (i.e., estimated best) features
        are assigned rank 1.

    Examples
    --------
    The following example shows how to retrieve the 5 right informative
    features in the Friedman #1 dataset.

    >>> from scikits.learn.datasets import make_friedman1
    >>> from scikits.learn.feature_selection import RFE
    >>> from scikits.learn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFE(estimator, 5, step=1)
    >>> selector = selector.fit(X, y)
    >>> selector.support_
    array([ True,  True,  True,  True,  True, False, False, False, False, False], dtype=bool)
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

    References
    ----------
    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    def __init__(self, estimator, n_features_to_select, step=1):
        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.step = step

    def fit(self, X, y):
        """Fit the RFE model.

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
        n_features = X.shape[1]

        if 0.0 < self.step < 1.0
            step = int(self.step * n_features)
        else:
            step = int(self.step)
        assert step > 0

        estimator = self.estimator
        support_ = np.ones(n_features, dtype=np.bool)
        ranking_ = np.ones(n_features, dtype=np.int)

        # Elimination
        while np.sum(support_) > self.n_features_to_select:
            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator.fit(X[:, features], y)
            ranks = np.argsort(np.sum(estimator.coef_ ** 2, axis=0))

            # Eliminate the worse features
            threshold = min(step, np.sum(support_) - self.n_features_to_select)
            support_[features[ranks][:threshold]] = False
            ranking_[~support_] += 1

        # Set final attributes
        self.support_ = support_
        self.ranking_ = ranking_

        return self

    def transform(self, X, copy=True):
        """Reduce X to the features selected during the fit.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            Vector, where n_samples in the number of samples and
            n_features is the number of features.

        copy : boolean, optional (default=True)
            If True, return a new array whose values are copied from X.
            If False, return a masked view of X.
        """
        X_r = X[:, self.support_]
        return X_r.copy() if copy else X_r


class RFECV(RFE):
    """Feature ranking with recursive feature elimination and cross-validated
       selection of the best number of features.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a `fit` method that updates a
        `coef_` attribute that holds the fitted parameters. The first
        dimension of the `coef_` array must be equal to the number of features
        of the input dataset of the estimator. Important features must
        correspond to high absolute values in the `coef_` array.

        For instance, this is the case for most supervised learning
        algorithms such as Support Vector Classifiers and Generalized
        Linear Models from the `svm` and `linear_model` modules.

    step : int or float, optional (default=1)
        If int, then `step` corresponds to the number of features to remove
        at each iteration.
        If float, then `step` corresponds to the percentage (rounded down) of
        features to remove at each iteration. This value should be in (0, 1].

    cv : int or cross-validation generator, optional (default=None)
        If int, it is the number of folds.
        If None, 3-fold cross-validation is performed by default.
        Specific cross-validation objects can also be passed, see
        `scikits.learn.cross_val module` for details.

    loss_function : function, optional (default=None)
        The loss function to minimize by cross-validation. If None, then the
        score function of the estimator is maximized.

    Attributes
    ----------
    `support_` : array of shape [n_features]
        The mask of selected features.

    `ranking_` : array of shape [n_features]
        The feature ranking, such that `ranking_[i]` corresponds to the ranking
        position of the i-th feature. Selected (i.e., estimated best) features
        are assigned rank 1.

    Examples
    --------
    The following example shows how to retrieve the a-priori not known 5
    informative features in the Friedman #1 dataset.

    >>> from scikits.learn.datasets import make_friedman1
    >>> from scikits.learn.feature_selection import RFECV
    >>> from scikits.learn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFECV(estimator, step=1, cv=5)
    >>> selector = selector.fit(X, y)
    >>> selector.support_
    array([ True,  True,  True,  True,  True, False, False, False, False, False], dtype=bool)
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

    References
    ----------
    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    def __init__(self, estimator, step=1, cv=None, loss_func=None):
        self.estimator = estimator
        self.step = 1
        self.cv = cv
        self.loss_func = None

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
        for train, test in cv:
            # Compute a full ranking of the features
            ranking_ = rfe.fit(X[train], y[train]).ranking_

            # Score each subset of features
            for k in xrange(1, max(ranking_)):
                mask = ranking_ <= k
                self.estimator.fit(X[train][:, mask], y[train])

                if self.loss_func is None:
                    score_k = 1.0 - self.estimator.score(
                                  X[test][:, mask],
                                  y[test])
                else:
                    score_k = self.loss_func(
                                  y[test],
                                  self.estimator.predict(X[test][:, mask]))

                if not k in scores:
                    scores[k] = 0.0

                scores[k] += score_k

        # Pick the best number of features on average
        best_score = np.inf
        best_k = None

        for k, score in scores.iteritems():
            if score < best_score:
                best_score = score
                best_k = k

        # Re-execute an elimination with best_k over the whole set
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=best_k,
                  step=self.step)

        rfe.fit(X, y)
        self.support_ = rfe.support_
        self.ranking_ = rfe.ranking_

        return self
