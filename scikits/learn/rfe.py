"""Recursive feature elimination
for feature ranking
"""

import numpy as np
from .base import BaseEstimator
# from base import BaseEstimator

class RFE(BaseEstimator):
    """
    Feature ranking with Recursive feature elimination

    Parameters
    ----------
    estimator : object
         object 

    n_features : int
        Number of features to select

    percentage : float
        The percentage of features to remove at each iteration
        Should be between (0, 1].  By default 0.1 will be taken.

    Attributes
    ----------
    `support_` : array-like, shape = [nfeatures]
        Mask of estimated support

    Methods
    -------
    fit(X, y) : self
        Fit the model

    transform(X) : array
        Reduce X to support

    Examples
    --------
    >>> 

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
        """Fit the RFE model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [nsamples, nfeatures]
            Training vector, where nsamples in the number of samples and
            nfeatures is the number of features.
        y : array, shape = [nsamples]
            Target values (integers in classification, real numbers in
            regression)
        """
        n_features_total = X.shape[1]
        estimator = self.estimator
        support_ = np.ones(n_features_total, dtype=np.bool)
        while np.sum(support_) > self.n_features:
            estimator.fit(X[:,support_], y)
            # rank features based on coef_ (handle multi class)
            abs_coef_ = np.sum(estimator.coef_ ** 2, axis=0)
            sorted_abs_coef_ = np.sort(abs_coef_)
            thresh = sorted_abs_coef_[-self.n_features]
            support_[support_] = abs_coef_ > thresh
        self.support_ = support_
        return self

    def transform(self, X, copy=True):
        """Reduce X to the features selected during the fit

        Parameters
        ----------
        X : array-like, shape = [nsamples, nfeatures]
            Vector, where nsamples in the number of samples and
            nfeatures is the number of features.
        """
        X_r = X[:,self.support_]
        return X_r.copy() if copy else X_r

if __name__ == '__main__':
    from scikits.learn.svm import SVC
    from scikits.learn import datasets
    iris = datasets.load_iris()

    # Add the noisy data to the informative features
    X = iris.data
    y = iris.target

    svc = SVC(kernel='linear')
    rfe = RFE(estimator=svc, n_features=2, percentage=0.1)
    rfe.fit(X, y)
    print rfe.support_
