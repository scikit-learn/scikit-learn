# Authors: Fabian Pedregosa <fabian@fseoane.net>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD Style.

import numpy as np
from scipy import interpolate
from ..base import BaseEstimator, TransformerMixin, RegressorMixin
from ..utils import as_float_array, check_arrays


def isotonic_regression(y, weight=None, x_min=None, x_max=None):
    """solve the isotonic regression model:

        min Sum w_i (y_i - x_i) ** 2

        subject to x_min = x_1 <= x_2 ... <= x_n = x_max

    where each w_i is strictly positive and each y_i is an arbitrary
    real number.

    Parameters
    ----------
    y: iterable of floating-point values
        The data

    weight: iterable of floating-point values, optional, default: None
        Weights on each point of the regression.
        If None, weight is set to 1 (equal weights)

    x_min: optional, default: None
        if not None, set the lowest value of the fit to x_min

    x_max: optional, default: None
        if not None, set the highest value of the fit to x_max

    Returns
    -------
    x: list of floating-point values

    References
    ----------
    Isotonic Median Regression: A Linear Programming Approach
    Nilotpal Chakravarti
    Mathematics of Operations Research
    Vol. 14, No. 2 (May, 1989), pp. 303-308
    """
    if weight is None:
        weight = np.ones(len(y), dtype=y.dtype)
    if x_min is not None or x_max is not None:
        y = np.copy(y)
        weight = np.copy(weight)
        C = np.dot(weight, y * y) * 10  # upper bound on the cost function
        if x_min is not None:
            y[0] = x_min
            weight[0] = C
        if x_max is not None:
            y[-1] = x_max
            weight[-1] = C

    active_set = [(weight[i] * y[i], weight[i], [i, ])
                  for i in range(len(y))]
    current = 0

    while current < len(active_set) - 1:
        value0, value1, value2 = 0, 0, np.inf
        weight0, weight1, weight2 = 1, 1, 1
        while value0 * weight1 <= value1 * weight0 and \
              current < len(active_set) - 1:
            value0, weight0, idx0 = active_set[current]
            value1, weight1, idx1 = active_set[current + 1]
            if value0 * weight1 <= value1 * weight0:
                current += 1

        if current == len(active_set) - 1:
            break

        # merge two groups
        value0, weight0, idx0 = active_set.pop(current)
        value1, weight1, idx1 = active_set.pop(current)
        active_set.insert(current,
                          (value0 + value1,
                           weight0 + weight1, idx0 + idx1))
        while value2 * weight0 > value0 * weight2 and current > 0:
            value0, weight0, idx0 = active_set[current]
            value2, weight2, idx2 = active_set[current - 1]
            if weight0 * value2 >= weight2 * value0:
                active_set.pop(current)
                active_set[current - 1] = (value0 + value2,
                                  weight0 + weight2,
                                  idx0 + idx2)
                current -= 1

    solution = np.empty(len(y))
    for value, weight, idx in active_set:
        solution[idx] = value / weight
    return solution


class IsotonicRegression(BaseEstimator, TransformerMixin, RegressorMixin):
    """solve the isotonic regression optimization problem

    The isotonic regression optimization problem is defined by::
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

    References
    ----------
    Isotonic Median Regression: A Linear Programming Approach
    Nilotpal Chakravarti
    Mathematics of Operations Research
    Vol. 14, No. 2 (May, 1989), pp. 303-308
    """
    def __init__(self, x_min=None, x_max=None):
        self.x_min = x_min
        self.x_max = x_max

    def _check_fit_data(self, X, y, weight=None):
        if len(X.shape) != 1:
            raise ValueError("X should be a vector")

    def fit(self, X, y, weight=None):
        """Fit the model using X as training data

        Parameters
        ----------
        X: array-like, shape=(n_samples,)
            training data

        y: array-like, shape=(n_samples,)
            training target

        weight: array-like, shape=(n_samples,), optional, default: None
            weights. If set to None, all weights will be set to 1 (equal
            weights)

        Returns
        -------
        self; object
            returns an instance of self

        Note
        ----
        X doesn't influence the result of `fit`. It is however stored
        for future use, as `transform` needs X to interpolate new
        input data.
        """
        X, y, weight = check_arrays(X, y, weight, sparse_format='dense')
        y = as_float_array(y)
        self.X_ = as_float_array(X, copy=True)
        self._check_fit_data(self.X_, y, weight)
        self.y_ = isotonic_regression(y, weight, self.x_min, self.x_max)
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
        T = as_float_array(T)
        if len(T.shape) != 1:
            raise ValueError("X should be a vector")

        f = interpolate.interp1d(self.X_, self.y_, kind='linear',
                                 bounds_error=True)
        return f(T)

    def fit_transform(self, X, y, weight=None):
        """Transform by linear interpolation

        Parameters
        ----------
        X: array-like, shape=(n_samples,)
            training data

        y: array-like, shape=(n_samples,)
            training target

        weight: array-like, shape=(n_samples,), optional, default: None
            weights. If set to None, all weights will be equal to 1 (equal
            weights)

        Returns
        -------
        y_: array, shape=(n_samples,)
            The transformed data

        Note
        ----
        X doesn't influence the result of `fit_transform`. It is however stored
        for future use, as `transform` needs X to interpolate new input
        data.
        """
        self.fit(X, y, weight)
        return self.y_

    def predict(self, T):
        """Predict new data by linear interpolation along

        Parameters
        ----------
        T: array-like, shape=(n_samples,)
            data to transform

        Returns
        -------
        T_: array, shape=(n_samples,)
            The transformed data
        """
        return self.transform(T)
