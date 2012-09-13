# Authors: Fabian Pedregosa <fabian@fseoane.net>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD Style.

import numpy as np
from scipy import interpolate
from ..base import BaseEstimator, TransformerMixin
from ..utils import as_float_array, check_arrays


def isotonic_regression(y, w=None, x_min=None, x_max=None):
    """
    Solve the isotonic regression model:

        min Sum w_i (y_i - x_i) ** 2

        subject to x_min = x_1 <= x_2 ... <= x_n = x_max

    where each w_i is strictly positive and each y_i is an arbitrary
    real number.

    Parameters
    ----------
    y: iterable of floating-point values

    w: iterable of floating-point values
        If None w is set to 1 (no weights)

    x_min: optional, default: None
        if not None, set the lowest value of the fit to x_min

    x_max: optional, default: None
        if not None, set the highest value of the fit to x_max

    Returns
    -------
    x: list of floating-point values
    """
    if w is None:
        w = np.ones(len(y), dtype=y.dtype)
    if x_min is not None or x_max is not None:
        y = np.copy(y)
        w = np.copy(w)
        C = np.dot(w, y * y) * 10  # upper bound on the cost function
        if x_min is not None:
            y[0] = x_min
            w[0] = C
        if x_max is not None:
            y[-1] = x_max
            w[-1] = C

    J = [(w[i] * y[i], w[i], [i, ]) for i in range(len(y))]
    cur = 0

    while cur < len(J) - 1:
        v0, v1, v2 = 0, 0, np.inf
        w0, w1, w2 = 1, 1, 1
        while v0 * w1 <= v1 * w0 and cur < len(J) - 1:
            v0, w0, idx0 = J[cur]
            v1, w1, idx1 = J[cur + 1]
            if v0 * w1 <= v1 * w0:
                cur += 1

        if cur == len(J) - 1:
            break

        # merge two groups
        v0, w0, idx0 = J.pop(cur)
        v1, w1, idx1 = J.pop(cur)
        J.insert(cur, (v0 + v1, w0 + w1, idx0 + idx1))
        while v2 * w0 > v0 * w2 and cur > 0:
            v0, w0, idx0 = J[cur]
            v2, w2, idx2 = J[cur - 1]
            if w0 * v2 >= w2 * v0:
                J.pop(cur)
                J[cur - 1] = (v0 + v2, w0 + w2, idx0 + idx2)
                cur -= 1

    sol = np.empty(len(y))
    for v, w, idx in J:
        sol[idx] = v / w
    return sol


class IsotonicRegression(BaseEstimator, TransformerMixin):
    """Solve the isotonic regression optimization problem

    The isotonic regression optimization problem is defined by:
        min Sum w_i (y[i] - y_[i]) ** 2

        subject to y_min = y_[1] <= y_[2] ... <= y_[n] = y_max

    where each w_i is strictly positive and each y_i is an arbitrary
    real number.

    Parameters
    ----------
    y_min: optional, default: None
        if not None, set the lowest value of the fit to y_min

    y_max: optional, default: None
        if not None, set the highest value of the fit to y_max

    Attributes
    ----------
    `X_`: ndarray (n, )
        A copy of the input X

    `y_`: ndarray (n, )
        Estimated y

    Notes
    -----
    Isotonic Median Regression: A Linear Programming Approach
    Nilotpal Chakravarti
    Mathematics of Operations Research
    Vol. 14, No. 2 (May, 1989), pp. 303-308
    """
    def __init__(self, x_min=None, x_max=None):
        self.x_min = x_min
        self.x_max = x_max

    def _check_fit_data(self, X, y, w=None):
        if w is not None:
            if len(X) != len(w):
                raise ValueError("Shapes of X and w do not match")
        if len(X.shape) != 1:
            raise ValueError("X should be a vector")
        if len(X) != len(y):
            raise ValueError("X and y should have the same length")

    def fit(self, X, y, w=None):
        """Fit the model using X as training data

        Parameters
        ----------
        X: array-like, shape=(n_samples,)
            training data

        y: array-like, shape=(n_samples,)
            training target

        w: array-like, shape=(n_samples,)
            weights

        Returns
        -------
        self; object
            returns an instance of self
        """
        X, y = check_arrays(X, y, sparse_format='dense')
        y = as_float_array(y)
        self.X_ = as_float_array(X, copy=True)
        self._check_fit_data(self.X_, y, w)
        self.y_ = isotonic_regression(y, w, self.x_min, self.x_max)
        return self

    def transform(self, T):
        """Transform new data by linear interpolation along

        Parameters
        ----------
        T: array-like, shape=(n_samples,)
            data to transform

        Returns
        -------
        T_: array, shape=(n_samples,)
            The transformed data
        """
        if len(T.shape) != 1:
            raise ValueError("X should be a vector")

        f = interpolate.interp1d(self.X_, self.y_, kind='linear',
                                 bounds_error=True)
        return f(T)

    def fit_transform(self, X, y, w=None):
        """Transform by linear interpolation

        Parameters
        ----------
        X: array-like, shape=(n_samples,)
            training data

        y: array-like, shape=(n_samples,)
            training target

        w: array-like, shape=(n_samples,)
            weights

        Returns
        -------
        y_: array, shape=(n_samples,)
            The transformed data
        """
        self.fit(X, y, w)
        return self.y_
