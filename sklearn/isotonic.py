# Authors: Fabian Pedregosa <fabian@fseoane.net>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD Style.

import numpy as np
from scipy import interpolate
from .base import BaseEstimator, TransformerMixin, RegressorMixin
from .utils import as_float_array, check_arrays


def isotonic_regression(y, weight=None, y_min=None, y_max=None):
    """Solve the isotonic regression model::

        min sum w[i] (y[i] - y_[i]) ** 2

        subject to y_min = y_[1] <= y_[2] ... <= y_[n] = y_max

    where:
        - y[i] are inputs (real numbers)
        - y_[i] are fitted
        - w[i] are optional strictly positive weights (default to 1.0)

    Parameters
    ----------
    y : iterable of floating-point values
        The data.

    weight : iterable of floating-point values, optional, default: None
        Weights on each point of the regression.
        If None, weight is set to 1 (equal weights).

    y_min : optional, default: None
        If not None, set the lowest value of the fit to y_min.

    y_max : optional, default: None
        If not None, set the highest value of the fit to y_max.

    Returns
    -------
    `y_` : list of floating-point values
        Isotonic fit of y.

    References
    ----------
    "Active set algorithms for isotonic regression; A unifying framework"
    by Michael J. Best and Nilotpal Chakravarti, section 3.
    """
    if weight is None:
        weight = np.ones(len(y), dtype=y.dtype)
    if y_min is not None or y_max is not None:
        y = np.copy(y)
        weight = np.copy(weight)
        C = np.dot(weight, y * y) * 10  # upper bound on the cost function
        if y_min is not None:
            y[0] = y_min
            weight[0] = C
        if y_max is not None:
            y[-1] = y_max
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
                active_set[current - 1] = (value0 + value2, weight0 + weight2,
                                           idx0 + idx2)
                current -= 1

    solution = np.empty(len(y))
    for value, weight, idx in active_set:
        solution[idx] = value / weight
    return solution


class IsotonicRegression(BaseEstimator, TransformerMixin, RegressorMixin):
    """Isotonic regression model.

    The isotonic regression optimization problem is defined by::

        min sum w_i (y[i] - y_[i]) ** 2

        subject to y_[i] <= y_[j] whenever X[i] <= X[j]
        and min(y_) = y_min, max(y_) = y_max

    where:
        - ``y[i]`` are inputs (real numbers)
        - ``y_[i]`` are fitted
        - ``X`` specifies the order.
          If ``X`` is non-decreasing then ``y_`` is non-decreasing.
        - ``w[i]`` are optional strictly positive weights (default to 1.0)

    Parameters
    ----------
    y_min : optional, default: None
        If not None, set the lowest value of the fit to y_min.

    y_max : optional, default: None
        If not None, set the highest value of the fit to y_max.

    Attributes
    ----------
    `X_` : ndarray (n_samples, )
        A copy of the input X.

    `y_` : ndarray (n_samples, )
        Isotonic fit of y.

    References
    ----------
    Isotonic Median Regression: A Linear Programming Approach
    Nilotpal Chakravarti
    Mathematics of Operations Research
    Vol. 14, No. 2 (May, 1989), pp. 303-308
    """
    def __init__(self, y_min=None, y_max=None):
        self.y_min = y_min
        self.y_max = y_max

    def _check_fit_data(self, X, y, weight=None):
        if len(X.shape) != 1:
            raise ValueError("X should be a vector")

    def fit(self, X, y, weight=None):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like, shape=(n_samples,)
            Training data.

        y : array-like, shape=(n_samples,)
            Training target.

        weight : array-like, shape=(n_samples,), optional, default: None
            Weights. If set to None, all weights will be set to 1 (equal
            weights).

        Returns
        -------
        self : object
            Returns an instance of self.

        Notes
        -----
        X is stored for future use, as `transform` needs X to interpolate
        new input data.
        """
        X, y, weight = check_arrays(X, y, weight, sparse_format='dense')
        y = as_float_array(y)
        self._check_fit_data(X, y, weight)
        order = np.argsort(X)
        self.X_ = as_float_array(X[order], copy=False)
        self.y_ = isotonic_regression(y[order], weight, self.y_min, self.y_max)
        return self

    def transform(self, T):
        """Transform new data by linear interpolation

        Parameters
        ----------
        T : array-like, shape=(n_samples,)
            Data to transform.

        Returns
        -------
        `T_` : array, shape=(n_samples,)
            The transformed data
        """
        T = as_float_array(T)
        if len(T.shape) != 1:
            raise ValueError("X should be a vector")

        f = interpolate.interp1d(self.X_, self.y_, kind='linear',
                                 bounds_error=True)
        return f(T)

    def fit_transform(self, X, y, weight=None):
        """Fit model and transform y by linear interpolation.

        Parameters
        ----------
        X : array-like, shape=(n_samples,)
            Training data.

        y : array-like, shape=(n_samples,)
            Training target.

        weight : array-like, shape=(n_samples,), optional, default: None
            Weights. If set to None, all weights will be equal to 1 (equal
            weights).

        Returns
        -------
        `y_` : array, shape=(n_samples,)
            The transformed data.

        Notes
        -----
        X doesn't influence the result of `fit_transform`. It is however stored
        for future use, as `transform` needs X to interpolate new input
        data.
        """
        X, y, weight = check_arrays(X, y, weight, sparse_format='dense')
        y = as_float_array(y)
        self._check_fit_data(X, y, weight)
        order = np.lexsort((y, X))
        order_inv = np.argsort(order)
        self.X_ = as_float_array(X[order], copy=False)
        self.y_ = isotonic_regression(y[order], weight, self.y_min, self.y_max)
        return self.y_[order_inv]

    def predict(self, T):
        """Predict new data by linear interpolation.

        Parameters
        ----------
        T : array-like, shape=(n_samples,)
            Data to transform.

        Returns
        -------
        `T_` : array, shape=(n_samples,)
            Transformed data.
        """
        return self.transform(T)
