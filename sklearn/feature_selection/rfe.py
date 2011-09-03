# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Vincent Michel <vincent.michel@inria.fr>
#
# License: BSD Style.

"""Recursive feature elimination for feature ranking"""

import numpy as np
from ..base import BaseEstimator


class RFE(BaseEstimator):
    """Feature ranking with Recursive feature elimination

    Parameters
    ----------
    estimator : object
         A supervised learning estimator with a fit method that updates a
         coef_ attributes that holds the fitted parameters. The first
         dimension of the coef_ array must be equal n_features an important
         features must yield high absolute values in the coef_ array.

         For instance this is the case for most supervised learning
         algorithms such as Support Vector Classifiers and Generalized
         Linear Models from the svm and linear_model package.

    n_features : int
        Number of features to select

    percentage : float
        The percentage of features to remove at each iteration
        Should be between (0, 1].  By default 0.1 will be taken.

    Attributes
    ----------
    `support_` : array-like, shape = [n_features]
        Mask of estimated support

    `ranking_` : array-like, shape = [n_features]
        Mask of the ranking of features

    Methods
    -------
    fit(X, y) : self
        Fit the model

    transform(X) : array
        Reduce X to support

    Examples
    --------
    >>> # TODO!

    References
    ----------
    Guyon, I., Weston, J., Barnhill, S., & Vapnik, V. (2002). Gene
    selection for cancer classification using support vector
    machines. Mach. Learn., 46(1-3), 389--422.
    """

    def __init__(self, estimator=None, n_features=None, percentage=0.1):
        self.n_features = n_features
        self.percentage = percentage
        self.estimator = estimator

    def fit(self, X, y):
        """Fit the RFE model

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
        """
        n_features_total = X.shape[1]
        estimator = self.estimator
        support_ = np.ones(n_features_total, dtype=np.bool)
        ranking_ = np.ones(n_features_total, dtype=np.int)
        while np.sum(support_) > self.n_features:
            estimator.fit(X[:, support_], y)
            # rank features based on coef_ (handle multi class)
            abs_coef_ = np.sum(estimator.coef_ ** 2, axis=0)
            sorted_abs_coef_ = np.sort(abs_coef_)
            threshold = sorted_abs_coef_[np.int(np.sum(support_) *
                                         self.percentage)]
            support_[support_] = abs_coef_ > threshold
            ranking_[support_] += 1
        self.support_ = support_
        self.ranking_ = ranking_
        return self

    def transform(self, X, copy=True):
        """Reduce X to the features selected during the fit

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Vector, where n_samples in the number of samples and
            n_features is the number of features.
        """
        X_r = X[:, self.support_]
        return X_r.copy() if copy else X_r


class RFECV(RFE):
    """Feature ranking with Recursive feature elimination and cross validation

    Parameters
    ----------
    estimator : object
         A supervised learning estimator with a fit method that updates a
         coef_ attributes that holds the fitted parameters. The first
         dimension of the coef_ array must be equal n_features an important
         features must yield high absolute values in the coef_ array.

         For instance this is the case for most supervised learning
         algorithms such as Support Vector Classifiers and Generalized
         Linear Models from the svm and linear_model package.

    n_features : int
        Number of features to select

    percentage : float
        The percentage of features to remove at each iteration
        Should be between (0, 1].  By default 0.1 will be taken.

    Attributes
    ----------
    `support_` : array-like, shape = [n_features]
        Mask of estimated support

    `ranking_` : array-like, shape = [n_features]
        Mask of the ranking of features

    Methods
    -------
    fit(X, y) : self
        Fit the model

    transform(X) : array
        Reduce X to support

    Examples
    --------
    >>> # TODO!

    References
    ----------
    Guyon, I., Weston, J., Barnhill, S., & Vapnik, V. (2002). Gene
    selection for cancer classification using support vector
    machines. Mach. Learn., 46(1-3), 389--422.
    """

    def __init__(self, estimator=None, n_features=None, percentage=0.1,
                  loss_func=None):
        self.n_features = n_features
        self.percentage = percentage
        self.estimator = estimator
        self.loss_func = loss_func

    def fit(self, X, y, cv=None):
        """Fit the RFE model with cross-validation

        The final size of the support is tuned by cross validation.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
        cv : cross-validation instance
        """
        rfe = RFE(estimator=self.estimator, n_features=self.n_features,
                          percentage=self.percentage)
        self.ranking_ = rfe.fit(X, y).ranking_
        clf = self.estimator
        n_models = np.max(self.ranking_)
        self.cv_scores_ = np.zeros(n_models)
        self.n_features_ = np.bincount(self.ranking_)[::-1].cumsum()[-2::-1]

        for train, test in cv:
            ranking_ = rfe.fit(X[train], y[train]).ranking_

            # assert n_models == np.max(ranking_)
            for k in range(n_models):
                mask = ranking_ >= (k+1)
                clf.fit(X[train][:, mask], y[train])
                y_pred = clf.predict(X[test][:, mask])
                self.cv_scores_[k] += self.loss_func(y[test], y_pred)

        # Take the best model (if multiple models have the same accuracy
        # use the last one ie the one with minimum number of features)
        min_score = n_models - np.argmin(self.cv_scores_[::-1])
        self.support_ = self.ranking_ >= min_score
        return self
