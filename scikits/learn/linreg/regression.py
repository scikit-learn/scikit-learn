# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

# $Id$
"""
Regression using linear models: Linear, Ridge.
"""

import numpy as np
import scipy.linalg

class LinearRegression(object):
    """
    Linear Regression.

    Parameters
    ----------
    This class takes no parameters

    This is just plain linear regression wrapped is a Predictor object.
    """

    def fit(self,X,Y):
        """
        Fit linear model
        """
        self.w, self.residues, self.rank, self.singular = \
                scipy.linalg.lstsq(X, Y)
        return self

    def predict(self, T):
        """
        Predict using linear model
        """
        return np.dot(T, self.w)


class RidgeRegression(object):
    """
    Ridge regression.


    Parameters
    ----------
    alpha : ridge parameter. Small positive values of alpha improve
    the coditioning of the problem and reduce the variance of the
    estimates.
    
    Examples
    --------
    # With more samples than features
    >>> import numpy as np
    >>> nsamples, nfeatures = 10, 5
    >>> np.random.seed(0)
    >>> Y = np.random.randn(nsamples)
    >>> X = np.random.randn(nsamples, nfeatures)
    >>> clf = Ridge(alpha=alpha)
    >>> clf.fit(X, Y)
    ?

    See also
    --------
    http://scikit-learn.sourceforge.net/doc/modules/glm.html

    """

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def fit(self, X, Y):
        """Fit Ridge regression model"""
        nsamples, nfeatures = X.shape

        if nsamples > nfeatures:
            # w = inv(X^t X + alpha*Id) * X.T y
            self.w = scipy.linalg.solve(np.dot(X.T,X) + self.alpha * np.eye(nfeatures),
                                  np.dot(X.T,y))
        else:
            # w = X.T * inv(X X^t + alpha*Id) y
            self.w = np.dot(X.T,
                    scipy.linalg.solve(np.dot(X, X.T) + self.alpha * np.eye(nsamples), y))

        return self

    def predict(self, X):
        """Predict using Linear Model
        """
        return np.dot(X,self.w)
