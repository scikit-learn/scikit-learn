"""Module for RBF interpolation."""
import warnings
from types import GenericAlias

import numpy as np
from scipy.spatial import KDTree

from . import _rbfinterp_np
from . import _rbfinterp_xp

from scipy._lib._array_api import (
    _asarray, array_namespace, xp_size, is_numpy, xp_capabilities
)
import scipy._lib.array_api_extra as xpx


__all__ = ["RBFInterpolator"]


# These RBFs are implemented.
_AVAILABLE = {
    "linear",
    "thin_plate_spline",
    "cubic",
    "quintic",
    "multiquadric",
    "inverse_multiquadric",
    "inverse_quadratic",
    "gaussian"
    }


# The shape parameter does not need to be specified when using these RBFs.
_SCALE_INVARIANT = {"linear", "thin_plate_spline", "cubic", "quintic"}


# For RBFs that are conditionally positive definite of order m, the interpolant
# should include polynomial terms with degree >= m - 1. Define the minimum
# degrees here. These values are from Chapter 8 of Fasshauer's "Meshfree
# Approximation Methods with MATLAB". The RBFs that are not in this dictionary
# are positive definite and do not need polynomial terms.
_NAME_TO_MIN_DEGREE = {
    "multiquadric": 0,
    "linear": 0,
    "thin_plate_spline": 1,
    "cubic": 1,
    "quintic": 2
    }


def _get_backend(xp):
    if is_numpy(xp):
        return _rbfinterp_np
    return _rbfinterp_xp


extra_note="""Only the default ``neighbors=None`` is Array API compatible.
    If a non-default value of ``neighbors`` is given, the behavior is NumPy -only.

"""

@xp_capabilities(
    skip_backends=[
        ("dask.array", "linalg.lu is broken; array_api_extra#488"),
        ("array_api_strict", "array-api#977, diag, view")
    ],
    extra_note=extra_note
)
class RBFInterpolator:
    """Radial basis function interpolator in N ≥ 1 dimensions.

    Parameters
    ----------
    y : (npoints, ndims) array_like
        2-D array of data point coordinates.
    d : (npoints, ...) array_like
        N-D array of data values at `y`. The length of `d` along the first
        axis must be equal to the length of `y`. Unlike some interpolators, the
        interpolation axis cannot be changed.
    neighbors : int, optional
        If specified, the value of the interpolant at each evaluation point
        will be computed using only this many nearest data points. All the data
        points are used by default.
    smoothing : float or (npoints, ) array_like, optional
        Smoothing parameter. The interpolant perfectly fits the data when this
        is set to 0. For large values, the interpolant approaches a least
        squares fit of a polynomial with the specified degree. Default is 0.
    kernel : str, optional
        Type of RBF. This should be one of

            - 'linear'               : ``-r``
            - 'thin_plate_spline'    : ``r**2 * log(r)``
            - 'cubic'                : ``r**3``
            - 'quintic'              : ``-r**5``
            - 'multiquadric'         : ``-sqrt(1 + r**2)``
            - 'inverse_multiquadric' : ``1/sqrt(1 + r**2)``
            - 'inverse_quadratic'    : ``1/(1 + r**2)``
            - 'gaussian'             : ``exp(-r**2)``

        Default is 'thin_plate_spline'.
    epsilon : float, optional
        Shape parameter that scales the input to the RBF. If `kernel` is
        'linear', 'thin_plate_spline', 'cubic', or 'quintic', this defaults to
        1 and can be ignored because it has the same effect as scaling the
        smoothing parameter. Otherwise, this must be specified.
    degree : int, optional
        Degree of the added polynomial. For some RBFs the interpolant may not
        be well-posed if the polynomial degree is too small. Those RBFs and
        their corresponding minimum degrees are

            - 'multiquadric'      : 0
            - 'linear'            : 0
            - 'thin_plate_spline' : 1
            - 'cubic'             : 1
            - 'quintic'           : 2

        The default value is the minimum degree for `kernel` or 0 if there is
        no minimum degree. Set this to -1 for no added polynomial.

    Notes
    -----
    An RBF is a scalar valued function in N-dimensional space whose value at
    :math:`x` can be expressed in terms of :math:`r=||x - c||`, where :math:`c`
    is the center of the RBF.

    An RBF interpolant for the vector of data values :math:`d`, which are from
    locations :math:`y`, is a linear combination of RBFs centered at :math:`y`
    plus a polynomial with a specified degree. The RBF interpolant is written
    as

    .. math::
        f(x) = K(x, y) a + P(x) b,

    where :math:`K(x, y)` is a matrix of RBFs with centers at :math:`y`
    evaluated at the points :math:`x`, and :math:`P(x)` is a matrix of
    monomials, which span polynomials with the specified degree, evaluated at
    :math:`x`. The coefficients :math:`a` and :math:`b` are the solution to the
    linear equations

    .. math::
        (K(x, y) + \\lambda I) a + P(y) b = d

    and

    .. math::
        P(y)^T a = 0,

    where :math:`\\lambda` is a non-negative smoothing parameter that controls
    how well we want to fit the data. The data are fit exactly when the
    smoothing parameter is 0.

    The above system is uniquely solvable if the following requirements are
    met:

        - :math:`P(y)` must have full column rank. :math:`P(y)` always has full
          column rank when `degree` is -1 or 0. When `degree` is 1,
          :math:`P(y)` has full column rank if the data point locations are not
          all collinear (N=2), coplanar (N=3), etc.
        - If `kernel` is 'multiquadric', 'linear', 'thin_plate_spline',
          'cubic', or 'quintic', then `degree` must not be lower than the
          minimum value listed above.
        - If `smoothing` is 0, then each data point location must be distinct.

    When using an RBF that is not scale invariant ('multiquadric',
    'inverse_multiquadric', 'inverse_quadratic', or 'gaussian'), an appropriate
    shape parameter must be chosen (e.g., through cross validation). Smaller
    values for the shape parameter correspond to wider RBFs. The problem can
    become ill-conditioned or singular when the shape parameter is too small.

    The memory required to solve for the RBF interpolation coefficients
    increases quadratically with the number of data points, which can become
    impractical when interpolating more than about a thousand data points.
    To overcome memory limitations for large interpolation problems, the
    `neighbors` argument can be specified to compute an RBF interpolant for
    each evaluation point using only the nearest data points.

    .. versionadded:: 1.7.0

    See Also
    --------
    NearestNDInterpolator
    LinearNDInterpolator
    CloughTocher2DInterpolator

    References
    ----------
    .. [1] Fasshauer, G., 2007. Meshfree Approximation Methods with Matlab.
        World Scientific Publishing Co.

    .. [2] http://amadeus.math.iit.edu/~fass/603_ch3.pdf

    .. [3] Wahba, G., 1990. Spline Models for Observational Data. SIAM.

    .. [4] http://pages.stat.wisc.edu/~wahba/stat860public/lect/lect8/lect8.pdf

    Examples
    --------
    Demonstrate interpolating scattered data to a grid in 2-D.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import RBFInterpolator
    >>> from scipy.stats.qmc import Halton

    >>> rng = np.random.default_rng()
    >>> xobs = 2*Halton(2, seed=rng).random(100) - 1
    >>> yobs = np.sum(xobs, axis=1)*np.exp(-6*np.sum(xobs**2, axis=1))

    >>> x1 = np.linspace(-1, 1, 50)
    >>> xgrid = np.asarray(np.meshgrid(x1, x1, indexing='ij'))
    >>> xflat = xgrid.reshape(2, -1).T     # make it a 2-D array
    >>> yflat = RBFInterpolator(xobs, yobs)(xflat)
    >>> ygrid = yflat.reshape(50, 50)

    >>> fig, ax = plt.subplots()
    >>> ax.pcolormesh(*xgrid, ygrid, vmin=-0.25, vmax=0.25, shading='gouraud')
    >>> p = ax.scatter(*xobs.T, c=yobs, s=50, ec='k', vmin=-0.25, vmax=0.25)
    >>> fig.colorbar(p)
    >>> plt.show()

    """

    # generic type compatibility with scipy-stubs
    __class_getitem__ = classmethod(GenericAlias)

    def __init__(self, y, d,
                 neighbors=None,
                 smoothing=0.0,
                 kernel="thin_plate_spline",
                 epsilon=None,
                 degree=None):
        xp = array_namespace(y, d, smoothing)
        _backend = _get_backend(xp)

        if neighbors is not None:
            if not is_numpy(xp):
                raise NotImplementedError(
                    "neighbors not None is numpy-only because it relies on KDTree"
                )

        y = _asarray(y, dtype=xp.float64, order="C", xp=xp)
        if y.ndim != 2:
            raise ValueError("`y` must be a 2-dimensional array.")

        ny, ndim = y.shape

        d = xp.asarray(d)
        if xp.isdtype(d.dtype, 'complex floating'):
            d_dtype = xp.complex128
        else:
            d_dtype = xp.float64
        d = _asarray(d, dtype=d_dtype, order="C", xp=xp)

        if d.shape[0] != ny:
            raise ValueError(
                f"Expected the first axis of `d` to have length {ny}."
                )

        d_shape = d.shape[1:]
        d = xp.reshape(d, (ny, -1))
        # If `d` is complex, convert it to a float array with twice as many
        # columns. Otherwise, the LHS matrix would need to be converted to
        # complex and take up 2x more memory than necessary.
        d = d.view(float)     # NB not Array API compliant (and jax copies)

        if isinstance(smoothing, int | float) or smoothing.shape == ():
            smoothing = xp.full(ny, smoothing, dtype=xp.float64)
        else:
            smoothing = _asarray(smoothing, dtype=float, order="C", xp=xp)
            if smoothing.shape != (ny,):
                raise ValueError(
                    "Expected `smoothing` to be a scalar or have shape "
                    f"({ny},)."
                    )

        kernel = kernel.lower()
        if kernel not in _AVAILABLE:
            raise ValueError(f"`kernel` must be one of {_AVAILABLE}.")

        if epsilon is None:
            if kernel in _SCALE_INVARIANT:
                epsilon = 1.0
            else:
                raise ValueError(
                    "`epsilon` must be specified if `kernel` is not one of "
                    f"{_SCALE_INVARIANT}."
                    )
        else:
            epsilon = float(epsilon)

        min_degree = _NAME_TO_MIN_DEGREE.get(kernel, -1)
        if degree is None:
            degree = max(min_degree, 0)
        else:
            degree = int(degree)
            if degree < -1:
                raise ValueError("`degree` must be at least -1.")
            elif -1 < degree < min_degree:
                warnings.warn(
                    f"`degree` should not be below {min_degree} except -1 "
                    f"when `kernel` is '{kernel}'."
                    f"The interpolant may not be uniquely "
                    f"solvable, and the smoothing parameter may have an "
                    f"unintuitive effect.",
                    UserWarning, stacklevel=2
                )

        if neighbors is None:
            nobs = ny
        else:
            # Make sure the number of nearest neighbors used for interpolation
            # does not exceed the number of observations.
            neighbors = int(min(neighbors, ny))
            nobs = neighbors

        powers = _backend._monomial_powers(ndim, degree, xp)
        # The polynomial matrix must have full column rank in order for the
        # interpolant to be well-posed, which is not possible if there are
        # fewer observations than monomials.
        if powers.shape[0] > nobs:
            raise ValueError(
                f"At least {powers.shape[0]} data points are required when "
                f"`degree` is {degree} and the number of dimensions is {ndim}."
                )

        if neighbors is None:
            shift, scale, coeffs = _backend._build_and_solve_system(
                y, d, smoothing, kernel, epsilon, powers,
                xp
            )

            # Make these attributes private since they do not always exist.
            self._shift = shift
            self._scale = scale
            self._coeffs = coeffs

        else:
            self._tree = KDTree(y)

        self.y = y
        self.d = d
        self.d_shape = d_shape
        self.d_dtype = d_dtype
        self.neighbors = neighbors
        self.smoothing = smoothing
        self.kernel = kernel
        self.epsilon = epsilon
        self.powers = powers
        self._xp = xp

    def __setstate__(self, state):
        tpl1, tpl2 = state
        (self.y, self.d, self.d_shape, self.d_dtype, self.neighbors,
         self.smoothing, self.kernel, self.epsilon, self.powers) = tpl1

        if self.neighbors is None:
            self._shift, self._scale, self._coeffs = tpl2
        else:
            self._tree, = tpl2

        self._xp = array_namespace(self.y, self.d, self.smoothing)

    def __getstate__(self):
        tpl = (self.y, self.d, self.d_shape, self.d_dtype, self.neighbors,
               self.smoothing, self.kernel, self.epsilon, self.powers
        )

        if self.neighbors is None:
            tpl2 = (self._shift, self._scale, self._coeffs)
        else:
            tpl2 = (self._tree,)

        return (tpl, tpl2)

    def _chunk_evaluator(
            self,
            x,
            y,
            shift,
            scale,
            coeffs,
            memory_budget=1000000
    ):
        """
        Evaluate the interpolation while controlling memory consumption.
        We chunk the input if we need more memory than specified.

        Parameters
        ----------
        x : (Q, N) float ndarray
            array of points on which to evaluate
        y: (P, N) float ndarray
            array of points on which we know function values
        shift: (N, ) ndarray
            Domain shift used to create the polynomial matrix.
        scale : (N,) float ndarray
            Domain scaling used to create the polynomial matrix.
        coeffs: (P+R, S) float ndarray
            Coefficients in front of basis functions
        memory_budget: int
            Total amount of memory (in units of sizeof(float)) we wish
            to devote for storing the array of coefficients for
            interpolated points. If we need more memory than that, we
            chunk the input.

        Returns
        -------
        (Q, S) float ndarray
        Interpolated array
        """
        _backend = _get_backend(self._xp)

        nx, ndim = x.shape
        if self.neighbors is None:
            nnei = y.shape[0]
        else:
            nnei = self.neighbors
        # in each chunk we consume the same space we already occupy
        chunksize = memory_budget // (self.powers.shape[0] + nnei) + 1
        if chunksize <= nx:
            out = self._xp.empty((nx, self.d.shape[1]), dtype=self._xp.float64)
            for i in range(0, nx, chunksize):
                chunk = _backend.compute_interpolation(
                    x[i:i + chunksize, :],
                    y,
                    self.kernel,
                    self.epsilon,
                    self.powers,
                    shift,
                    scale,
                    coeffs,
                    self._xp
                )
                out = xpx.at(out, (slice(i, i + chunksize), slice(None,))).set(chunk)
        else:
            out = _backend.compute_interpolation(
                x,
                y,
                self.kernel,
                self.epsilon,
                self.powers,
                shift,
                scale,
                coeffs,
                self._xp
            )
        return out

    def __call__(self, x):
        """Evaluate the interpolant at `x`.

        Parameters
        ----------
        x : (npts, ndim) array_like
            Evaluation point coordinates.

        Returns
        -------
        ndarray, shape (npts, )
            Values of the interpolant at `x`.

        """
        x = _asarray(x, dtype=self._xp.float64, order="C", xp=self._xp)
        if x.ndim != 2:
            raise ValueError("`x` must be a 2-dimensional array.")

        nx, ndim = x.shape
        if ndim != self.y.shape[1]:
            raise ValueError("Expected the second axis of `x` to have length "
                             f"{self.y.shape[1]}.")

        # Our memory budget for storing RBF coefficients is
        # based on how many floats in memory we already occupy
        # If this number is below 1e6 we just use 1e6
        # This memory budget is used to decide how we chunk
        # the inputs
        memory_budget = max(xp_size(x) + xp_size(self.y) + xp_size(self.d), 1_000_000)

        if self.neighbors is None:
            out = self._chunk_evaluator(
                x,
                self.y,
                self._shift,
                self._scale,
                self._coeffs,
                memory_budget=memory_budget)
        else:
            # XXX: this relies on KDTree, hence is numpy-only until KDTree is converted
            _build_and_solve_system = _get_backend(np)._build_and_solve_system

            # Get the indices of the k nearest observation points to each
            # evaluation point.
            _, yindices = self._tree.query(x, self.neighbors)
            if self.neighbors == 1:
                # `KDTree` squeezes the output when neighbors=1.
                yindices = yindices[:, None]

            # Multiple evaluation points may have the same neighborhood of
            # observation points. Make the neighborhoods unique so that we only
            # compute the interpolation coefficients once for each
            # neighborhood.
            yindices = np.sort(yindices, axis=1)
            yindices, inv = np.unique(yindices, return_inverse=True, axis=0)
            inv = np.reshape(inv, (-1,))  # flatten, we need 1-D indices
            # `inv` tells us which neighborhood will be used by each evaluation
            # point. Now we find which evaluation points will be using each
            # neighborhood.
            xindices = [[] for _ in range(len(yindices))]
            for i, j in enumerate(inv):
                xindices[j].append(i)

            out = np.empty((nx, self.d.shape[1]), dtype=float)
            for xidx, yidx in zip(xindices, yindices):
                # `yidx` are the indices of the observations in this
                # neighborhood. `xidx` are the indices of the evaluation points
                # that are using this neighborhood.
                xnbr = x[xidx]
                ynbr = self.y[yidx]
                dnbr = self.d[yidx]
                snbr = self.smoothing[yidx]
                shift, scale, coeffs = _build_and_solve_system(
                    ynbr,
                    dnbr,
                    snbr,
                    self.kernel,
                    self.epsilon,
                    self.powers,
                    np
                )
                out[xidx] = self._chunk_evaluator(
                    xnbr,
                    ynbr,
                    shift,
                    scale,
                    coeffs,
                    memory_budget=memory_budget)

        out = out.view(self.d_dtype)    # NB not Array API compliant (and jax copies)
        out = self._xp.reshape(out, (nx, ) + self.d_shape)
        return out
