# Authors: Fabian Pedregosa <fabian@fseoane.net>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD 3 clause

import numpy as np
from scipy import interpolate
from scipy.stats import spearmanr
from .base import BaseEstimator, TransformerMixin, RegressorMixin
from .utils import as_float_array, check_array, check_consistent_length
from ._isotonic import _isotonic_regression
import warnings
import math


def check_increasing(x, y):
    """Determine whether y is monotonically correlated with x.

    y is found increasing or decreasing with respect to x based on a Spearman
    correlation test.

    Parameters
    ----------
    x : array-like, shape=(n_samples,)
            Training data.

    y : array-like, shape=(n_samples,)
        Training target.

    Returns
    -------
    `increasing_bool` : boolean
        Whether the relationship is increasing or decreasing.

    Notes
    -----
    The Spearman correlation coefficient is estimated from the data, and the
    sign of the resulting estimate is used as the result.

    In the event that the 95% confidence interval based on Fisher transform
    spans zero, a warning is raised.

    References
    ----------
    Fisher transformation. Wikipedia.
    http://en.wikipedia.org/w/index.php?title=Fisher_transformation
    """

    # Calculate Spearman rho estimate and set return accordingly.
    rho, _ = spearmanr(x, y)
    increasing_bool = rho >= 0

    # Run Fisher transform to get the rho CI, but handle rho=+/-1
    if rho not in [-1.0, 1.0]:
        F = 0.5 * math.log((1. + rho) / (1. - rho))
        F_se = 1 / math.sqrt(len(x) - 3)

        # Use a 95% CI, i.e., +/-1.96 S.E.
        # http://en.wikipedia.org/wiki/Fisher_transformation
        rho_0 = math.tanh(F - 1.96 * F_se)
        rho_1 = math.tanh(F + 1.96 * F_se)

        # Warn if the CI spans zero.
        if np.sign(rho_0) != np.sign(rho_1):
            warnings.warn("Confidence interval of the Spearman "
                          "correlation coefficient spans zero. "
                          "Determination of ``increasing`` may be "
                          "suspect.")

    return increasing_bool


def isotonic_regression(y, sample_weight=None, y_min=None, y_max=None,
                        increasing=True):
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

    sample_weight : iterable of floating-point values, optional, default: None
        Weights on each point of the regression.
        If None, weight is set to 1 (equal weights).

    y_min : optional, default: None
        If not None, set the lowest value of the fit to y_min.

    y_max : optional, default: None
        If not None, set the highest value of the fit to y_max.

    increasing : boolean, optional, default: True
        Whether to compute ``y_`` is increasing (if set to True) or decreasing
        (if set to False)

    Returns
    -------
    y_ : list of floating-point values
        Isotonic fit of y.

    References
    ----------
    "Active set algorithms for isotonic regression; A unifying framework"
    by Michael J. Best and Nilotpal Chakravarti, section 3.
    """
    y = np.asarray(y, dtype=np.float)
    if sample_weight is None:
        sample_weight = np.ones(len(y), dtype=y.dtype)
    else:
        sample_weight = np.asarray(sample_weight, dtype=np.float)
    if not increasing:
        y = y[::-1]
        sample_weight = sample_weight[::-1]

    if y_min is not None or y_max is not None:
        y = np.copy(y)
        sample_weight = np.copy(sample_weight)
        # upper bound on the cost function
        C = np.dot(sample_weight, y * y) * 10
        if y_min is not None:
            y[0] = y_min
            sample_weight[0] = C
        if y_max is not None:
            y[-1] = y_max
            sample_weight[-1] = C

    solution = np.empty(len(y))
    y_ = _isotonic_regression(y, sample_weight, solution)
    if increasing:
        return y_
    else:
        return y_[::-1]


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

    increasing : boolean or string, optional, default: True
        If boolean, whether or not to fit the isotonic regression with y
        increasing or decreasing.

        The string value "auto" determines whether y should
        increase or decrease based on the Spearman correlation estimate's
        sign.

    out_of_bounds : string, optional, default: "nan"
        The ``out_of_bounds`` parameter handles how x-values outside of the
        training domain are handled.  When set to "nan", predicted y-values
        will be NaN.  When set to "clip", predicted y-values will be
        set to the value corresponding to the nearest train interval endpoint.
        When set to "raise", allow ``interp1d`` to throw ValueError.


    Attributes
    ----------
    X_ : ndarray (n_samples, )
        A copy of the input X.

    y_ : ndarray (n_samples, )
        Isotonic fit of y.

    X_min_ : float
        Minimum value of input array `X_` for left bound.

    X_max_ : float
        Maximum value of input array `X_` for right bound.

    f_ : function
        The stepwise interpolating function that covers the domain
        X_.

    References
    ----------
    Isotonic Median Regression: A Linear Programming Approach
    Nilotpal Chakravarti
    Mathematics of Operations Research
    Vol. 14, No. 2 (May, 1989), pp. 303-308
    """
    def __init__(self, y_min=None, y_max=None, increasing=True,
                 out_of_bounds='nan'):
        self.y_min = y_min
        self.y_max = y_max
        self.increasing = increasing
        self.out_of_bounds = out_of_bounds

    def _check_fit_data(self, X, y, sample_weight=None):
        if len(X.shape) != 1:
            raise ValueError("X should be a 1d array")

    def _build_f(self, X, y):
        """Build the f_ interp1d function."""

        # Handle the out_of_bounds argument by setting bounds_error
        if self.out_of_bounds not in ["raise", "nan", "clip"]:
            raise ValueError("The argument ``out_of_bounds`` must be in "
                             "'nan', 'clip', 'raise'; got {0}"
                             .format(self.out_of_bounds))

        bounds_error = self.out_of_bounds == "raise"
        self.f_ = interpolate.interp1d(X, y, kind='linear',
                                       bounds_error=bounds_error)

    def _build_y(self, X, y, sample_weight):
        """Build the y_ IsotonicRegression."""
        check_consistent_length(X, y, sample_weight)
        X, y = [check_array(x, ensure_2d=False) for x in [X, y]]
        if sample_weight is not None:
            sample_weight = check_array(sample_weight, ensure_2d=False)

        y = as_float_array(y)
        self._check_fit_data(X, y, sample_weight)

        # Determine increasing if auto-determination requested
        if self.increasing == 'auto':
            self.increasing_ = check_increasing(X, y)
        else:
            self.increasing_ = self.increasing

        order = np.lexsort((y, X))
        order_inv = np.argsort(order)
        self.X_ = as_float_array(X[order], copy=False)
        self.y_ = isotonic_regression(y[order], sample_weight, self.y_min,
                                      self.y_max, increasing=self.increasing_)

        return order_inv

    def fit(self, X, y, sample_weight=None):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like, shape=(n_samples,)
            Training data.

        y : array-like, shape=(n_samples,)
            Training target.

        sample_weight : array-like, shape=(n_samples,), optional, default: None
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
        # Build y_
        self._build_y(X, y, sample_weight)

        # Handle the left and right bounds on X
        self.X_min_ = np.min(self.X_)
        self.X_max_ = np.max(self.X_)

        # Build f_
        self._build_f(self.X_, self.y_)

        return self

    def transform(self, T):
        """Transform new data by linear interpolation

        Parameters
        ----------
        T : array-like, shape=(n_samples,)
            Data to transform.

        Returns
        -------
        T_ : array, shape=(n_samples,)
            The transformed data
        """
        T = as_float_array(T)
        if len(T.shape) != 1:
            raise ValueError("Isotonic regression input should be a 1d array")

        # Handle the out_of_bounds argument by clipping if needed
        if self.out_of_bounds not in ["raise", "nan", "clip"]:
            raise ValueError("The argument ``out_of_bounds`` must be in "
                             "'nan', 'clip', 'raise'; got {0}"
                             .format(self.out_of_bounds))

        if self.out_of_bounds == "clip":
            T = np.clip(T, self.X_min_, self.X_max_)
        return self.f_(T)

    def fit_transform(self, X, y, sample_weight=None):
        """Fit model and transform y by linear interpolation.

        Parameters
        ----------
        X : array-like, shape=(n_samples,)
            Training data.

        y : array-like, shape=(n_samples,)
            Training target.

        sample_weight : array-like, shape=(n_samples,), optional, default: None
            Weights. If set to None, all weights will be equal to 1 (equal
            weights).

        Returns
        -------
        y_ : array, shape=(n_samples,)
            The transformed data.

        Notes
        -----
        X doesn't influence the result of `fit_transform`. It is however stored
        for future use, as `transform` needs X to interpolate new input
        data.
        """
        # Build y_
        order_inv = self._build_y(X, y, sample_weight)

        # Handle the left and right bounds on X
        self.X_min_ = np.min(self.X_)
        self.X_max_ = np.max(self.X_)

        # Build f_
        self._build_f(self.X_, self.y_)

        return self.y_[order_inv]

    def predict(self, T):
        """Predict new data by linear interpolation.

        Parameters
        ----------
        T : array-like, shape=(n_samples,)
            Data to transform.

        Returns
        -------
        T_ : array, shape=(n_samples,)
            Transformed data.
        """
        return self.transform(T)
