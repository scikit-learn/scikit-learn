# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$
"""Implementation of Ridge regression
"""

import numpy as np
import scipy.linalg as linalg

class Ridge(object):
    """Ridge"""

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def fit(self, X, y):
        """Fit Ridge regression model"""
        nsamples, nfeatures = X.shape

        if nsamples > nfeatures:
            # w = inv(X^t X + alpha*Id) * X.T y
            self.w = linalg.solve(np.dot(X.T,X) + self.alpha * np.eye(nfeatures),
                                  np.dot(X.T,y))
        else:
            # w = X.T * inv(X X^t + alpha*Id) y
            self.w = np.dot(X.T,
                    linalg.solve(np.dot(X, X.T) + self.alpha * np.eye(nsamples), y))

        return self

    def predict(self, X):
        """Predict with Linear Model
        """
        y = np.dot(X,self.w)
        return y

if __name__ == '__main__':
    """Tests Ridge regression
    """

    alpha = 1.0

    # With more samples than features
    nsamples, nfeatures = 10, 5
    np.random.seed(0)
    y = np.random.randn(nsamples)
    X = np.random.randn(nsamples, nfeatures)

    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)

    # With more features than samples
    nsamples, nfeatures = 5, 10
    np.random.seed(0)
    y = np.random.randn(nsamples)
    X = np.random.randn(nsamples, nfeatures)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
