"""
This file contains preprocessing tools based on polynomials.
"""
import numbers

import numpy as np
from scipy.interpolate import BSpline

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.fixes import linspace
from ..utils.validation import check_is_fitted, FLOAT_DTYPES


__all__ = [
    "SplineTransformer",
]


# TODO:
# - sparse support (either scipy or own cython solution)?
# - extrapolation (cyclic)
class SplineTransformer(TransformerMixin, BaseEstimator):
    """Generate univariate B-spline bases for features.

    Generate a new feature matrix consisting of
    `n_splines=n_knots + degree - 1` spline basis functions (B-splines) of
    polynomial order=`degree` for each feature.

    Read more in the :ref:`User Guide <spline_transformer>`.

    .. versionadded:: 1.0

    Parameters
    ----------
    n_knots : int, default=5
        Number of knots of the splines if `knots` equals one of
        {'uniform', 'quantile'}. Must be larger or equal 2.

    degree : int, default=3
        The polynomial degree of the spline basis. Must be a non-negative
        integer.

    knots : {'uniform', 'quantile'} or array-like of shape \
        (n_knots, n_features), default='uniform'
        Set knot positions such that first knot <= features <= last knot.

        - If 'uniform', `n_knots` number of knots are distributed uniformly
          from min to max values of the features.
        - If 'quantile', they are distributed uniformly along the quantiles of
          the features.
        - If an array-like is given, it directly specifies the sorted knot
          positions including the boundary knots. Note that, internally,
          `degree` number of knots are added before the first knot, the same
          after the last knot.

    extrapolation : {'error', 'constant', 'linear', 'continue'}, \
        default='constant'
        If 'error', values outside the min and max values of the training
        features raises a `ValueError`. If 'constant', the value of the
        splines at minimum and maximum value of the features is used as
        constant extrapolation. If 'linear', a linear extrapolation is used.
        If 'continue', the splines are extrapolated as is, i.e. option
        `extrapolate=True` in :class:`scipy.interpolate.BSpline`.

    include_bias : bool, default=True
        If True (default), then the last spline element inside the data range
        of a feature is dropped. As B-splines sum to one over the spline basis
        functions for each data point, they implicitly include a bias term,
        i.e. a column of ones. It acts as an intercept term in a linear models.

    order : {'C', 'F'}, default='C'
        Order of output array. 'F' order is faster to compute, but may slow
        down subsequent estimators.

    Attributes
    ----------
    bsplines_ : list of shape (n_features,)
        List of BSplines objects, one for each feature.

    n_features_in_ : int
        The total number of input features.

    n_features_out_ : int
        The total number of output features, which is computed as
        `n_features * n_splines`, where `n_splines` is
        the number of bases elements of the B-splines, `n_knots + degree - 1`.
        If `include_bias=False`, then it is only
        `n_features * (n_splines - 1)`.

    See Also
    --------
    KBinsDiscretizer : Transformer that bins continuous data into intervals.

    PolynomialFeatures : Transformer that generates polynomial and interaction
        features.

    Notes
    -----
    High degrees and a high number of knots can cause overfitting.

    See :ref:`examples/linear_model/plot_polynomial_interpolation.py
    <sphx_glr_auto_examples_linear_model_plot_polynomial_interpolation.py>`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.preprocessing import SplineTransformer
    >>> X = np.arange(6).reshape(6, 1)
    >>> spline = SplineTransformer(degree=2, n_knots=3)
    >>> spline.fit_transform(X)
    array([[0.5 , 0.5 , 0.  , 0.  ],
           [0.18, 0.74, 0.08, 0.  ],
           [0.02, 0.66, 0.32, 0.  ],
           [0.  , 0.32, 0.66, 0.02],
           [0.  , 0.08, 0.74, 0.18],
           [0.  , 0.  , 0.5 , 0.5 ]])
    """

    def __init__(
        self,
        n_knots=5,
        degree=3,
        *,
        knots="uniform",
        extrapolation="constant",
        include_bias=True,
        order="C",
    ):
        self.n_knots = n_knots
        self.degree = degree
        self.knots = knots
        self.extrapolation = extrapolation
        self.include_bias = include_bias
        self.order = order

    @staticmethod
    def _get_base_knot_positions(X, n_knots=10, knots="uniform"):
        """Calculate base knot positions.

        Base knots such that first knot <= feature <= last knot. For the
        B-spline construction with scipy.interpolate.BSpline, 2*degree knots
        beyond the base interval are added.

        Returns
        -------
        knots : ndarray of shape (n_knots, n_features), dtype=np.float64
            Knot positions (points) of base interval.
        """
        if knots == "quantile":
            knots = np.percentile(
                X,
                100
                * np.linspace(start=0, stop=1, num=n_knots, dtype=np.float64),
                axis=0,
            )
        else:
            # knots == 'uniform':
            # Note that the variable `knots` has already been validated and
            # `else` is therefore safe.
            x_min = np.amin(X, axis=0)
            x_max = np.amax(X, axis=0)
            knots = linspace(
                start=x_min,
                stop=x_max,
                num=n_knots,
                endpoint=True,
                dtype=np.float64,
            )

        return knots

    def get_feature_names(self, input_features=None):
        """Return feature names for output features.

        Parameters
        ----------
        input_features : list of str of shape (n_features,), default=None
            String names for input features if available. By default,
            "x0", "x1", ... "xn_features" is used.

        Returns
        -------
        output_feature_names : list of str of shape (n_output_features,)
        """
        n_splines = self.bsplines_[0].c.shape[0]
        if input_features is None:
            input_features = ["x%d" % i for i in range(self.n_features_in_)]
        feature_names = []
        for i in range(self.n_features_in_):
            for j in range(n_splines - 1 + self.include_bias):
                feature_names.append(f"{input_features[i]}_sp_{j}")
        return feature_names

    def fit(self, X, y=None):
        """Compute knot positions of splines.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data.

        y : None
            Ignored.

        Returns
        -------
        self : object
            Fitted transformer.
        """
        X = self._validate_data(
            X,
            reset=True,
            accept_sparse=False,
            ensure_min_samples=2,
            ensure_2d=True,
        )
        n_samples, n_features = X.shape

        if not (
            isinstance(self.degree, numbers.Integral) and self.degree >= 0
        ):
            raise ValueError("degree must be a non-negative integer.")

        if not (
            isinstance(self.n_knots, numbers.Integral) and self.n_knots >= 2
        ):
            raise ValueError("n_knots must be a positive integer >= 2.")

        if isinstance(self.knots, str) and self.knots in [
            "uniform",
            "quantile",
        ]:
            base_knots = self._get_base_knot_positions(
                X, n_knots=self.n_knots, knots=self.knots
            )
        else:
            base_knots = check_array(self.knots)
            if base_knots.shape[0] < 2:
                raise ValueError(
                    "Number of knots, knots.shape[0], must be >= " "2."
                )
            elif base_knots.shape[1] != n_features:
                raise ValueError("knots.shape[1] == n_features is violated.")
            elif not np.all(np.diff(base_knots, axis=0) > 0):
                raise ValueError("knots must be sorted without duplicates.")

        if self.extrapolation not in (
            "error",
            "constant",
            "linear",
            "continue",
        ):
            raise ValueError(
                "extrapolation must be one of 'error', "
                "'constant', 'linear' or 'continue'."
            )

        if not isinstance(self.include_bias, (bool, np.bool_)):
            raise ValueError("include_bias must be bool.")

        # number of knots for base interval
        n_knots = base_knots.shape[0]
        # number of splines basis functions
        n_splines = n_knots + self.degree - 1
        degree = self.degree
        n_out = n_features * n_splines
        # We have to add degree number of knots below, and degree number knots
        # above the base knots in order to make the spline basis complete.
        # Eilers & Marx in "Flexible smoothing with B-splines and  penalties"
        # https://doi.org/10.1214/ss/1038425655 advice against repeating first
        # and last knot several times, which would have inferior behaviour at
        # boundaries if combined with a penalty (hence P-Spline). We follow
        # this advice even if our splines are unpenalized.
        # Meaning we do not:
        # knots = np.r_[np.tile(base_knots.min(axis=0), reps=[degree, 1]),
        #              base_knots,
        #              np.tile(base_knots.max(axis=0), reps=[degree, 1])
        #              ]
        # Instead, we reuse the distance of the 2 fist/last knots.
        dist_min = base_knots[1] - base_knots[0]
        dist_max = base_knots[-1] - base_knots[-2]
        knots = np.r_[
            linspace(
                base_knots[0] - degree * dist_min,
                base_knots[0] - dist_min,
                num=degree,
            ),
            base_knots,
            linspace(
                base_knots[-1] + dist_max,
                base_knots[-1] + degree * dist_max,
                num=degree,
            ),
        ]

        # With a diagonal coefficient matrix, we get back the spline basis
        # elements, i.e. the design matrix of the spline.
        # Note, BSpline appreciates C-contiguous float64 arrays as c=coef.
        coef = np.eye(n_knots + self.degree - 1, dtype=np.float64)
        extrapolate = self.extrapolation == "continue"
        bsplines = [
            BSpline.construct_fast(
                knots[:, i], coef, self.degree, extrapolate=extrapolate
            )
            for i in range(n_features)
        ]
        self.bsplines_ = bsplines

        self.n_features_out_ = n_out - n_features * (1 - self.include_bias)
        return self

    def transform(self, X):
        """Transform each feature data to B-splines.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to transform.

        Returns
        -------
        XBS : ndarray of shape (n_samples, n_features * n_splines)
            The matrix of features, where n_splines is the number of bases
            elements of the B-splines, n_knots + degree - 1.
        """
        check_is_fitted(self)

        X = self._validate_data(
            X, reset=False, accept_sparse=False, ensure_2d=True
        )

        n_samples, n_features = X.shape
        n_splines = self.bsplines_[0].c.shape[0]
        degree = self.degree

        # Note that scipy BSpline returns float64 arrays and converts input
        # x=X[:, i] to c-contiguous float64.
        n_out = self.n_features_out_ + n_features * (1 - self.include_bias)
        if X.dtype in FLOAT_DTYPES:
            dtype = X.dtype
        else:
            dtype = np.float64
        XBS = np.zeros((n_samples, n_out), dtype=dtype, order=self.order)

        for i in range(n_features):
            spl = self.bsplines_[i]

            if self.extrapolation in ("continue", "error"):
                XBS[:, (i * n_splines):((i + 1) * n_splines)] = spl(X[:, i])
            else:
                xmin = spl.t[degree]
                xmax = spl.t[-degree - 1]
                mask = (xmin <= X[:, i]) & (X[:, i] <= xmax)
                XBS[mask, (i * n_splines):((i + 1) * n_splines)] = spl(
                    X[mask, i]
                )

            # Note for extrapolation:
            # 'continue' is already returned as is by scipy BSplines
            if self.extrapolation == "error":
                # BSpline with extrapolate=False does not raise an error, but
                # output np.nan.
                if np.any(
                    np.isnan(XBS[:, (i * n_splines):((i + 1) * n_splines)])
                ):
                    raise ValueError(
                        "X contains values beyond the limits of the knots."
                    )
            elif self.extrapolation == "constant":
                # Set all values beyond xmin and xmax to the value of the
                # spline basis functions at those two positions.
                # Only the first degree and last degree number of splines
                # have non-zero values at the boundaries.

                # spline values at boundaries
                f_min = spl(xmin)
                f_max = spl(xmax)
                mask = X[:, i] < xmin
                if np.any(mask):
                    XBS[
                        mask, (i * n_splines):(i * n_splines + degree)
                    ] = f_min[:degree]

                mask = X[:, i] > xmax
                if np.any(mask):
                    XBS[
                        mask,
                        ((i + 1) * n_splines - degree):((i + 1) * n_splines),
                    ] = f_max[-degree:]
            elif self.extrapolation == "linear":
                # Continue the degree first and degree last spline bases
                # linearly beyond the boundaries, with slope = derivative at
                # the boundary.
                # Note that all others have derivative = value = 0 at the
                # boundaries.

                # spline values at boundaries
                f_min, f_max = spl(xmin), spl(xmax)
                # spline derivatives = slopes at boundaries
                fp_min, fp_max = spl(xmin, nu=1), spl(xmax, nu=1)
                # Compute the linear continuation.
                if degree <= 1:
                    # For degree=1, the derivative of 2nd spline is not zero at
                    # boundary. For degree=0 it is the same as 'constant'.
                    degree += 1
                for j in range(degree):
                    mask = X[:, i] < xmin
                    if np.any(mask):
                        XBS[mask, i * n_splines + j] = (
                            f_min[j] + (X[mask, i] - xmin) * fp_min[j]
                        )

                    mask = X[:, i] > xmax
                    if np.any(mask):
                        k = n_splines - 1 - j
                        XBS[mask, i * n_splines + k] = (
                            f_max[k] + (X[mask, i] - xmax) * fp_max[k]
                        )

        if self.include_bias:
            return XBS
        else:
            # We throw away one spline basis per feature.
            # We chose the last one.
            indices = [
                j for j in range(XBS.shape[1]) if (j + 1) % n_splines != 0
            ]
            return XBS[:, indices]
