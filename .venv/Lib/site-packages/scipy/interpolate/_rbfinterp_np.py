import numpy as np
from numpy.linalg import LinAlgError
from scipy.linalg.lapack import dgesv  # type: ignore[attr-defined]
from ._rbfinterp_common import _monomial_powers_impl

from ._rbfinterp_pythran import (
    _build_system as _pythran_build_system,
    _build_evaluation_coefficients as _pythran_build_evaluation_coefficients,
    _polynomial_matrix as _pythran_polynomial_matrix
)


# trampolines for pythran-compiled functions to drop the `xp` argument
def _build_evaluation_coefficients(
    x, y, kernel, epsilon, powers, shift, scale, xp
):
    return _pythran_build_evaluation_coefficients(
        x, y, kernel, epsilon, powers, shift, scale
    )

def polynomial_matrix(x, powers, xp):
    return _pythran_polynomial_matrix(x, powers)


def _monomial_powers(ndim, degree, xp):
    out = _monomial_powers_impl(ndim, degree)
    out = np.asarray(out, dtype=np.int64)
    if len(out) == 0:
        out = out.reshape(0, ndim)
    return out


def _build_system(y, d, smoothing, kernel, epsilon, powers, xp):
    return _pythran_build_system(y, d, smoothing, kernel, epsilon, powers)


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
    _, _, coeffs, info = dgesv(lhs, rhs, overwrite_a=True, overwrite_b=True)
    if info < 0:
        raise ValueError(f"The {-info}-th argument had an illegal value.")
    elif info > 0:
        msg = "Singular matrix."
        nmonos = powers.shape[0]
        if nmonos > 0:
            pmat = polynomial_matrix((y - shift)/scale, powers, xp)
            rank = np.linalg.matrix_rank(pmat)
            if rank < nmonos:
                msg = (
                    "Singular matrix. The matrix of monomials evaluated at "
                    "the data point coordinates does not have full column "
                    f"rank ({rank}/{nmonos})."
                    )

        raise LinAlgError(msg)

    return shift, scale, coeffs

def compute_interpolation(x, y, kernel, epsilon, powers, shift, scale, coeffs, xp):
    vec = _build_evaluation_coefficients(
        x, y, kernel, epsilon, powers, shift, scale, xp
    )
    return vec @ coeffs
