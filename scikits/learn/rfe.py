"""Recursive feature elimination
for feature ranking
"""

import numpy as np
from .base import BaseEstimator
# from base import BaseEstimator

class RFE(BaseEstimator):
    """Recursive feature elimination

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

    def transform(self, X):
        return X[:,self.support_]

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
