import numpy as np
from scipy import linalg

from .base import LinearModel

class Ridge(LinearModel):
    """
    Ridge regression.

    Parameters
    ----------
    alpha : float
        Small positive values of alpha improve the coditioning of the
        problem and reduce the variance of the estimates.
    fit_intercept : boolean
        wether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    Examples
    --------
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> Y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = Ridge(alpha=1.0)
    >>> clf.fit(X, Y)
    Ridge(alpha=1.0, fit_intercept=True)
    """

    def __init__(self, alpha=1.0, fit_intercept=True):
        self.alpha = alpha
        self.fit_intercept = fit_intercept


    def fit(self, X, Y, **params):
        """
        Fit Ridge regression model

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data
        Y : numpy array of shape [n_samples]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)

        X = np.asanyarray(X, dtype=np.float)
        Y = np.asanyarray(Y, dtype=np.float)

        n_samples, n_features = X.shape

        X, Y, Xmean, Ymean = self._center_data (X, Y)

        if n_samples > n_features:
            # w = inv(X^t X + alpha*Id) * X.T y
            self.coef_ = linalg.solve(
                np.dot(X.T, X) + self.alpha * np.eye(n_features),
                np.dot(X.T, Y))
        else:
            # w = X.T * inv(X X^t + alpha*Id) y
            self.coef_ = np.dot(X.T, linalg.solve(
                np.dot(X, X.T) + self.alpha * np.eye(n_samples), Y))

        self._set_intercept(Xmean, Ymean)
        return self


