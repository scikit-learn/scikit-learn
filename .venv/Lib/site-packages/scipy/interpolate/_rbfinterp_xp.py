"""
'Generic' Array API backend for RBF interpolation.

The general logic is this: `_rbfinterp.py` implements the user API and calls
into either `_rbfinterp_np` (the "numpy backend"), or `_rbfinterp_xp` (the
"generic backend".

The numpy backend offloads performance-critical computations to the
pythran-compiled `_rbfinterp_pythran` extension. This way, the call chain is

    _rbfinterp.py <-- _rbfinterp_np.py <-- _rbfinterp_pythran.py

The "generic" backend here is a drop-in replacement of the API of
`_rbfinterp_np.py` for use in `_rbfinterp.py` with non-numpy arrays.

The implementation closely follows `_rbfinterp_np + _rbfinterp_pythran`, with
the following differences:

  -  We used vectorized code not explicit loops in `_build_system` and
     `_build_evaluation_coefficients`; this is more torch/jax friendly;
  - RBF kernels are also "vectorized" and not scalar: they receive an
    array of norms not a single norm;
  - RBF kernels accept an extra xp= argument;

In general, we would prefer less code duplication. The main blocker ATM is
that pythran cannot compile functions with an xp= argument where xp is numpy.
"""
from numpy.linalg import LinAlgError
from ._rbfinterp_common import _monomial_powers_impl


def _monomial_powers(ndim, degree, xp):
    out = _monomial_powers_impl(ndim, degree)
    out = xp.asarray(out)
    if out.shape[0] == 0:
        out = xp.reshape(out, (0, ndim))
    return out


def _build_and_solve_system(y, d, smoothing, kernel, epsilon, powers, xp):
    """Build and solve the RBF interpolation system of equations.

    Parameters
    ----------
    y : (P, N) float ndarray
        Data point coordinates.
    d : (P, S) float ndarray
        Data values at `y`.
    smoothing : (P,) float ndarray
        Smoothing parameter for each data point.
    kernel : str
        Name of the RBF.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.

    Returns
    -------
    coeffs : (P + R, S) float ndarray
        Coefficients for each RBF and monomial.
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.

    """
    lhs, rhs, shift, scale = _build_system(
        y, d, smoothing, kernel, epsilon, powers, xp
        )
    try:
        coeffs = xp.linalg.solve(lhs, rhs)
    except Exception:
        # Best-effort attempt to emit a helpful message.
        # `_rbfinterp_np` backend gives better diagnostics; it is hard to
        # match it in a backend-agnostic way: e.g. jax emits no error at all,
        # and instead returns an array of nans for a singular `lhs`.
        msg = "Singular matrix"
        nmonos = powers.shape[0]
        if nmonos > 0:
            pmat = polynomial_matrix((y - shift)/scale, powers, xp=xp)
            rank = xp.linalg.matrix_rank(pmat)
            if rank < nmonos:
                msg = (
                    "Singular matrix. The matrix of monomials evaluated at "
                    "the data point coordinates does not have full column "
                    f"rank ({rank}/{nmonos})."
                    )
        raise LinAlgError(msg)

    return shift, scale, coeffs


def linear(r, xp):
    return -r


def thin_plate_spline(r, xp):
    # NB: changed w.r.t. pythran, vectorized
    return xp.where(r == 0, 0, r**2 * xp.log(r))


def cubic(r, xp):
    return r**3


def quintic(r, xp):
    return -r**5


def multiquadric(r, xp):
    return -xp.sqrt(r**2 + 1)


def inverse_multiquadric(r, xp):
    return 1.0 / xp.sqrt(r**2 + 1.0)


def inverse_quadratic(r, xp):
    return 1.0 / (r**2 + 1.0)


def gaussian(r, xp):
    return xp.exp(-r**2)


NAME_TO_FUNC = {
   "linear": linear,
   "thin_plate_spline": thin_plate_spline,
   "cubic": cubic,
   "quintic": quintic,
   "multiquadric": multiquadric,
   "inverse_multiquadric": inverse_multiquadric,
   "inverse_quadratic": inverse_quadratic,
   "gaussian": gaussian
   }


def kernel_matrix(x, kernel_func, xp):
    """Evaluate RBFs, with centers at `x`, at `x`."""
    return kernel_func(
        xp.linalg.vector_norm(x[None, :, :] - x[:, None, :], axis=-1), xp
    )


def polynomial_matrix(x, powers, xp):
    """Evaluate monomials, with exponents from `powers`, at `x`."""
    return xp.prod(x[:, None, :] ** powers, axis=-1)


def _build_system(y, d, smoothing, kernel, epsilon, powers, xp):
    """Build the system used to solve for the RBF interpolant coefficients.

    Parameters
    ----------
    y : (P, N) float ndarray
        Data point coordinates.
    d : (P, S) float ndarray
        Data values at `y`.
    smoothing : (P,) float ndarray
        Smoothing parameter for each data point.
    kernel : str
        Name of the RBF.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.

    Returns
    -------
    lhs : (P + R, P + R) float ndarray
        Left-hand side matrix.
    rhs : (P + R, S) float ndarray
        Right-hand side matrix.
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.

    """
    s = d.shape[1]
    r = powers.shape[0]
    kernel_func = NAME_TO_FUNC[kernel]

    # Shift and scale the polynomial domain to be between -1 and 1
    mins = xp.min(y, axis=0)
    maxs = xp.max(y, axis=0)
    shift = (maxs + mins)/2
    scale = (maxs - mins)/2
    # The scale may be zero if there is a single point or all the points have
    # the same value for some dimension. Avoid division by zero by replacing
    # zeros with ones.
    scale = xp.where(scale == 0.0, 1.0, scale)

    yeps = y*epsilon
    yhat = (y - shift)/scale

    out_kernels  = kernel_matrix(yeps, kernel_func, xp)
    out_poly = polynomial_matrix(yhat, powers, xp)

    lhs = xp.concat(
        [
         xp.concat((out_kernels, out_poly), axis=1),
         xp.concat((out_poly.T, xp.zeros((r, r))), axis=1)
        ]
    , axis=0) + xp.diag(xp.concat([smoothing, xp.zeros(r)]))

    rhs = xp.concat([d, xp.zeros((r, s))], axis=0)

    return lhs, rhs, shift, scale


def _build_evaluation_coefficients(
    x, y, kernel, epsilon, powers, shift, scale, xp
):
    """Construct the coefficients needed to evaluate
    the RBF.

    Parameters
    ----------
    x : (Q, N) float ndarray
        Evaluation point coordinates.
    y : (P, N) float ndarray
        Data point coordinates.
    kernel : str
        Name of the RBF.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.
    shift : (N,) float ndarray
        Shifts the polynomial domain for numerical stability.
    scale : (N,) float ndarray
        Scales the polynomial domain for numerical stability.

    Returns
    -------
    (Q, P + R) float ndarray

    """
    kernel_func = NAME_TO_FUNC[kernel]

    yeps = y*epsilon
    xeps = x*epsilon
    xhat = (x - shift)/scale

    # NB: changed w.r.t. pythran
    vec = xp.concat(
        [
            kernel_func(
                xp.linalg.vector_norm(
                    xeps[:, None, :] - yeps[None, :, :], axis=-1
                ), xp
            ),
            xp.prod(xhat[:, None, :] ** powers, axis=-1)
        ], axis=-1
    )

    return vec


def compute_interpolation(x, y, kernel, epsilon, powers, shift, scale, coeffs, xp):
    vec = _build_evaluation_coefficients(
        x, y, kernel, epsilon, powers, shift, scale, xp
    )
    return vec @ coeffs
