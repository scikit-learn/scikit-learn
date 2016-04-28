# Authors: Fabian Pedregosa <fabian@fseoane.net>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Nelle Varoquaux <nelle.varoquaux@gmail.com>
# License: BSD 3 clause

import numpy as np
from scipy import interpolate
from scipy.stats import spearmanr
from .base import BaseEstimator, TransformerMixin, RegressorMixin
from .utils import as_float_array, check_array, check_consistent_length
from .utils import deprecated
from .utils.fixes import astype
from ._isotonic import _isotonic_regression, _make_unique
import warnings
import math


__all__ = ['check_increasing', 'isotonic_regression',
           'IsotonicRegression']


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
    https://en.wikipedia.org/wiki/Fisher_transformation
    """

    # Calculate Spearman rho estimate and set return accordingly.
    rho, _ = spearmanr(x, y)
    increasing_bool = rho >= 0

    # Run Fisher transform to get the rho CI, but handle rho=+/-1
    if rho not in [-1.0, 1.0]:
        F = 0.5 * math.log((1. + rho) / (1. - rho))
        F_se = 1 / math.sqrt(len(x) - 3)

        # Use a 95% CI, i.e., +/-1.96 S.E.
        # https://en.wikipedia.org/wiki/Fisher_transformation
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

    Read more in the :ref:`User Guide <isotonic>`.

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
    y = np.asarray(y, dtype=np.float64)
    if sample_weight is None:
        sample_weight = np.ones(len(y), dtype=y.dtype)
    else:
        sample_weight = np.asarray(sample_weight, dtype=np.float64)
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

    Read more in the :ref:`User Guide <isotonic>`.

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
    X_min_ : float
        Minimum value of input array `X_` for left bound.

    X_max_ : float
        Maximum value of input array `X_` for right bound.

    f_ : function
        The stepwise interpolating function that covers the domain `X_`.

    Notes
    -----
    Ties are broken using the secondary method from Leeuw, 1977.

    References
    ----------
    Isotonic Median Regression: A Linear Programming Approach
    Nilotpal Chakravarti
    Mathematics of Operations Research
    Vol. 14, No. 2 (May, 1989), pp. 303-308

    Isotone Optimization in R : Pool-Adjacent-Violators
    Algorithm (PAVA) and Active Set Methods
    Leeuw, Hornik, Mair
    Journal of Statistical Software 2009

    Correctness of Kruskal's algorithms for monotone regression with ties
    Leeuw, Psychometrica, 1977
    """
    def __init__(self, y_min=None, y_max=None, increasing=True,
                 out_of_bounds='nan'):
        self.y_min = y_min
        self.y_max = y_max
        self.increasing = increasing
        self.out_of_bounds = out_of_bounds

    @property
    @deprecated("Attribute ``X_`` is deprecated in version 0.18 and will be"
                " removed in version 0.20.")
    def X_(self):
        return self._X_

    @X_.setter
    def X_(self, value):
        self._X_ = value

    @X_.deleter
    def X_(self):
        del self._X_

    @property
    @deprecated("Attribute ``y_`` is deprecated in version 0.18 and will"
                " be removed in version 0.20.")
    def y_(self):
        return self._y_

    @y_.setter
    def y_(self, value):
        self._y_ = value

    @y_.deleter
    def y_(self):
        del self._y_

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
        if len(y) == 1:
            # single y, constant prediction
            self.f_ = lambda x: y.repeat(x.shape)
        else:
            self.f_ = interpolate.interp1d(X, y, kind='linear',
                                           bounds_error=bounds_error)

    def _build_y(self, X, y, sample_weight, trim_duplicates=True):
        """Build the y_ IsotonicRegression."""
        check_consistent_length(X, y, sample_weight)
        X, y = [check_array(x, ensure_2d=False) for x in [X, y]]

        y = as_float_array(y)
        self._check_fit_data(X, y, sample_weight)

        # Determine increasing if auto-determination requested
        if self.increasing == 'auto':
            self.increasing_ = check_increasing(X, y)
        else:
            self.increasing_ = self.increasing

        # If sample_weights is passed, removed zero-weight values and clean
        # order
        if sample_weight is not None:
            sample_weight = check_array(sample_weight, ensure_2d=False)
            mask = sample_weight > 0
            X, y, sample_weight = X[mask], y[mask], sample_weight[mask]
        else:
            sample_weight = np.ones(len(y))

        order = np.lexsort((y, X))
        X, y, sample_weight = [astype(array[order], np.float64, copy=False)
                               for array in [X, y, sample_weight]]
        unique_X, unique_y, unique_sample_weight = _make_unique(
            X, y, sample_weight)

        # Store _X_ and _y_ to maintain backward compat during the deprecation
        # period of X_ and y_
        self._X_ = X = unique_X
        self._y_ = y = isotonic_regression(unique_y, unique_sample_weight,
                                           self.y_min, self.y_max,
                                           increasing=self.increasing_)

        # Handle the left and right bounds on X
        self.X_min_, self.X_max_ = np.min(X), np.max(X)

        if trim_duplicates:
            # Remove unnecessary points for faster prediction
            keep_data = np.ones((len(y),), dtype=bool)
            # Aside from the 1st and last point, remove points whose y values
            # are equal to both the point before and the point after it.
            keep_data[1:-1] = np.logical_or(
                np.not_equal(y[1:-1], y[:-2]),
                np.not_equal(y[1:-1], y[2:])
            )
            return X[keep_data], y[keep_data]
        else:
            # The ability to turn off trim_duplicates is only used to it make
            # easier to unit test that removing duplicates in y does not have
            # any impact the resulting interpolation function (besides
            # prediction speed).
            return X, y

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
        # Transform y by running the isotonic regression algorithm and
        # transform X accordingly.
        X, y = self._build_y(X, y, sample_weight)

        # It is necessary to store the non-redundant part of the training set
        # on the model to make it possible to support model persistence via
        # the pickle module as the object built by scipy.interp1d is not
        # picklable directly.
        self._necessary_X_, self._necessary_y_ = X, y

        # Build the interpolation function
        self._build_f(X, y)
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

    def __getstate__(self):
        """Pickle-protocol - return state of the estimator. """
        # copy __dict__
        state = dict(self.__dict__)
        # remove interpolation method
        state.pop('f_', None)
        return state

    def __setstate__(self, state):
        """Pickle-protocol - set state of the estimator.

        We need to rebuild the interpolation function.
        """
        self.__dict__.update(state)
        if hasattr(self, '_necessary_X_') and hasattr(self, '_necessary_y_'):
            self._build_f(self._necessary_X_, self._necessary_y_)
