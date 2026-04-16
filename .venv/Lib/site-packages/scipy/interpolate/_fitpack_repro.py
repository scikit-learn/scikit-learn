""" Replicate FITPACK's logic for constructing smoothing spline functions and curves.

    Currently provides analogs of splrep and splprep python routines, i.e.
    curfit.f and parcur.f routines (the drivers are fpcurf.f and fppara.f, respectively)

    The Fortran sources are from
    https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/

    .. [1] P. Dierckx, "Algorithms for smoothing data with periodic and
        parametric splines, Computer Graphics and Image Processing",
        20 (1982) 171-184.
        :doi:`10.1016/0146-664X(82)90043-0`.
    .. [2] P. Dierckx, "Curve and surface fitting with splines", Monographs on
         Numerical Analysis, Oxford University Press, 1993.
    .. [3] P. Dierckx, "An algorithm for smoothing, differentiation and integration
         of experimental data using spline functions",
         Journal of Computational and Applied Mathematics, vol. I, no 3, p. 165 (1975).
         https://doi.org/10.1016/0771-050X(75)90034-0
"""
import warnings
import operator
import numpy as np

from scipy._lib._array_api import array_namespace, concat_1d, xp_capabilities

from ._bsplines import (
    _not_a_knot, make_interp_spline, BSpline, fpcheck, _lsq_solve_qr,
    _lsq_solve_qr_for_root_rati_periodic, _periodic_knots
)
from . import _dierckx      # type: ignore[attr-defined]


#    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#    c  part 1: determination of the number of knots and their position     c
#    c  **************************************************************      c
#
# https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L31


# Hardcoded in curfit.f
TOL = 0.001
MAXIT = 20


def _get_residuals(x, y, t, k, w, periodic=False):
    # inline the relevant part of
    # >>> spl = make_lsq_spline(x, y, w=w2, t=t, k=k)
    # NB:
    #     1. y is assumed to be 2D here. For 1D case (parametric=False),
    #        the call must have been preceded by y = y[:, None] (cf _validate_inputs)
    #     2. We always sum the squares across axis=1:
    #         * For 1D (parametric=False), the last dimension has size one,
    #           so the summation is a no-op.
    #         * For 2D (parametric=True), the summation is actually how the
    #           'residuals' are defined, see Eq. (42) in Dierckx1982
    #           (the reference is in the docstring of `class F`) below.
    _, _, _, fp, residuals = _lsq_solve_qr(x, y, t, k, w, periodic=periodic)
    if np.isnan(residuals.sum()):
        raise ValueError(_iermesg[1])
    return residuals, fp

def _validate_bc_type(bc_type):
    if bc_type is None:
        return "not-a-knot"

    if bc_type not in ("not-a-knot", "periodic"):
        raise ValueError("Only 'not-a-knot' and 'periodic' "
                         "boundary conditions are recognised, "
                         f"found {bc_type}")

    return bc_type


def add_knot(x, t, k, residuals):
    """Add a new knot.

    (Approximately) replicate FITPACK's logic:
      1. split the `x` array into knot intervals, ``t(j+k) <= x(i) <= t(j+k+1)``
      2. find the interval with the maximum sum of residuals
      3. insert a new knot into the middle of that interval.

    NB: a new knot is in fact an `x` value at the middle of the interval.
    So *the knots are a subset of `x`*.

    This routine is an analog of
    https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpcurf.f#L190-L215
    (cf _split function)

    and https://github.com/scipy/scipy/blob/v1.11.4/scipy/interpolate/fitpack/fpknot.f
    """
    new_knot = _dierckx.fpknot(x, t, k, residuals)

    idx_t = np.searchsorted(t, new_knot)
    t_new = np.r_[t[:idx_t], new_knot, t[idx_t:]]
    return t_new


def _validate_inputs(x, y, w, k, s, xb, xe, parametric, periodic=False):
    """Common input validations for generate_knots and make_splrep.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if not x.flags.c_contiguous:
        x = x.copy()
    if not y.flags.c_contiguous:
        y = y.copy()

    if w is None:
        w = np.ones_like(x, dtype=float)
    else:
        w = np.asarray(w, dtype=float)
        if not w.flags.c_contiguous:
            w = w.copy()
        if w.ndim != 1:
            raise ValueError(f"{w.ndim = } not implemented yet.")
        if (w < 0).any():
            raise ValueError("Weights must be non-negative")
        if w.sum() == 0:
            raise ValueError("All weights are zero.")


    if y.ndim == 0 or y.ndim > 2:
        raise ValueError(f"{y.ndim = } not supported (must be 1 or 2.)")

    parametric = bool(parametric)
    if parametric:
        if y.ndim != 2:
            raise ValueError(f"{y.ndim = } != 2 not supported with {parametric =}.")
    else:
        if y.ndim != 1:
            raise ValueError(f"{y.ndim = } != 1 not supported with {parametric =}.")
        # all _impl functions expect y.ndim = 2
        y = y[:, None]

    if w.shape[0] != x.shape[0]:
        raise ValueError(f"Weights is incompatible: {w.shape =} != {x.shape}.")

    if x.shape[0] != y.shape[0]:
        raise ValueError(f"Data is incompatible: {x.shape = } and {y.shape = }.")
    if x.ndim != 1 or (x[1:] < x[:-1]).any():
        raise ValueError("Expect `x` to be an ordered 1D sequence.")

    k = operator.index(k)

    if s < 0:
        raise ValueError(f"`s` must be non-negative. Got {s = }")

    if xb is None:
        xb = min(x)
    if xe is None:
        xe = max(x)

    if periodic and not np.allclose(y[0], y[-1], atol=1e-15):
        raise ValueError("First and last points does not match which is required "
                         "for `bc_type='periodic'`.")

    return x, y, w, k, s, xb, xe


@xp_capabilities(cpu_only=True, jax_jit=False, allow_dask_compute=True)
def generate_knots(x, y, *, w=None, xb=None, xe=None,
                   k=3, s=0, nest=None, bc_type=None):
    """Generate knot vectors until the Least SQuares (LSQ) criterion is satified.

    Parameters
    ----------
    x, y : array_like
        The data points defining the curve ``y = f(x)``.
    w : array_like, optional
        Weights.
    xb : float, optional
        The boundary of the approximation interval. If None (default),
        is set to ``x[0]``.
    xe : float, optional
        The boundary of the approximation interval. If None (default),
        is set to ``x[-1]``.
    k : int, optional
        The spline degree. Default is cubic, ``k = 3``.
    s : float, optional
        The smoothing factor. Default is ``s = 0``.
    nest : int, optional
        Stop when at least this many knots are placed.
          ends are equivalent.
    bc_type : str, optional
        Boundary conditions.
        Default is `"not-a-knot"`.
        The following boundary conditions are recognized:

        * ``"not-a-knot"`` (default): The first and second segments are the
          same polynomial. This is equivalent to having ``bc_type=None``.
        * ``"periodic"``: The values and the first ``k-1`` derivatives at the
          ends are equivalent.

    Yields
    ------
    t : ndarray
        Knot vectors with an increasing number of knots.
        The generator is finite: it stops when the smoothing critetion is
        satisfied, or when then number of knots exceeds the maximum value:
        the user-provided `nest` or `x.size + k + 1` --- which is the knot vector
        for the interpolating spline.

    Examples
    --------
    Generate some noisy data and fit a sequence of LSQ splines:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import make_lsq_spline, generate_knots
    >>> rng = np.random.default_rng(12345)
    >>> x = np.linspace(-3, 3, 50)
    >>> y = np.exp(-x**2) + 0.1 * rng.standard_normal(size=50)

    >>> knots = list(generate_knots(x, y, s=1e-10))
    >>> for t in knots[::3]:
    ...     spl = make_lsq_spline(x, y, t)
    ...     xs = xs = np.linspace(-3, 3, 201)
    ...     plt.plot(xs, spl(xs), '-', label=f'n = {len(t)}', lw=3, alpha=0.7)
    >>> plt.plot(x, y, 'o', label='data')
    >>> plt.plot(xs, np.exp(-xs**2), '--')
    >>> plt.legend()

    Note that increasing the number of knots make the result follow the data
    more and more closely.

    Also note that a step of the generator may add multiple knots:

    >>> [len(t) for t in knots]
    [8, 9, 10, 12, 16, 24, 40, 48, 52, 54]

    Notes
    -----
    The routine generates successive knots vectors of increasing length, starting
    from ``2*(k+1)`` to ``len(x) + k + 1``, trying to make knots more dense
    in the regions where the deviation of the LSQ spline from data is large.

    When the maximum number of knots, ``len(x) + k + 1`` is reached
    (this happens when ``s`` is small and ``nest`` is large), the generator
    stops, and the last output is the knots for the interpolation with the
    not-a-knot boundary condition.

    Knots are located at data sites, unless ``k`` is even and the number of knots
    is ``len(x) + k + 1``. In that case, the last output of the generator
    has internal knots at Greville sites, ``(x[1:] + x[:-1]) / 2``.

    .. versionadded:: 1.15.0

    """
    xp = array_namespace(x, y, w)
    bc_type = _validate_bc_type(bc_type)
    periodic = bc_type == 'periodic'

    if s == 0:
        if nest is not None or w is not None:
            raise ValueError("s == 0 is interpolation only")
        # For special-case k=1 (e.g., Lyche and Morken, Eq.(2.16)),
        # _not_a_knot produces desired knot vector
        if periodic and k > 1:
            t = _periodic_knots(x, k)
        else:
            t = _not_a_knot(x, k)
        yield xp.asarray(t)
        return

    x, y, w, k, s, xb, xe = _validate_inputs(
        x, y, w, k, s, xb, xe, parametric=np.ndim(y) == 2,
        periodic=periodic
    )

    yield from _generate_knots_impl(x, y, w, xb, xe, k, s, nest, periodic, xp=xp)


def _generate_knots_impl(x, y, w, xb, xe, k, s, nest, periodic, xp=np):

    acc = s * TOL
    m = x.size    # the number of data points

    if nest is None:
        # the max number of knots. This is set in _fitpack_impl.py line 274
        # and fitpack.pyf line 198
        # Ref: https://github.com/scipy/scipy/blob/596b586e25e34bd842b575bac134b4d6924c6556/scipy/interpolate/_fitpack_impl.py#L260-L263
        if periodic:
            nest = max(m + 2*k, 2*k + 3)
        else:
            nest = max(m + k + 1, 2*k + 3)
    else:
        if nest < 2*(k + 1):
            raise ValueError(f"`nest` too small: {nest = } < 2*(k+1) = {2*(k+1)}.")

    if not periodic:
        nmin = 2*(k + 1)    # the number of knots for an LSQ polynomial approximation
        nmax = m + k + 1  # the number of knots for the spline interpolation
    else:
        # Ref: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L54
        nmin = 2*(k + 1)    # the number of knots for an LSQ polynomial approximation
        # Ref: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L61
        nmax = m + 2*k  # the number of knots for the spline interpolation

    per = xe - xb
    # Ref: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L107-L123
    # Computes fp0 for constant function
    if periodic:
        t = np.zeros(nmin, dtype=float)
        for i in range(0, k + 1):
            t[i] = x[0] - (k - i) * per
            t[i + k + 1] = x[m - 1] + i * per
        _, fp = _get_residuals(x, y, t, k, w, periodic=periodic)
        # For periodic splines, check whether constant function
        # satisfies accuracy criterion
        # Also if maximal number of nodes is equal to the minimal
        # then constant function is the direct solution
        # Ref: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L600-L610
        if fp - s < acc or nmax == nmin:
            yield t
            return
    else:
        fp = 0.0
        fpold = 0.0

    # start from no internal knots
    if not periodic:
        t = np.asarray([xb]*(k+1) + [xe]*(k+1), dtype=float)
    else:
        # Ref: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L131
        # Initialize knot vector `t` of size (2k + 3) with zeros.
        # The central knot `t[k + 1]` is seeded with the midpoint value from `x`.
        # Note that, in the `if periodic:` block (in the main loop below),
        # the boundary knots `t[k]` and `t[n - k - 1]` are set to the endpoints `xb` and
        # `xe`. Then, the surrounding knots on both ends are updated to ensure
        # periodicity.
        # Specifically:
        # - Left-side knots are mirrored from the right end minus the period (`per`).
        # - Right-side knots are mirrored from the left end plus the period.
        # These updates ensure that the knot vector wraps around correctly for periodic
        # B-spline fitting.
        t = np.zeros(2*k + 3, dtype=float)
        t[k + 1] = x[(m + 1)//2 - 1]
        nplus = 1
    n = t.shape[0]

    # c  main loop for the different sets of knots. m is a safe upper bound
    # c  for the number of trials.
    for iter in range(m):
        # Ref: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L147-L158
        if periodic:
            n = t.shape[0]
            t[k] = xb
            t[n - k - 1] = xe
            for j in range(1, k + 1):
                t[k - j] = t[n - k - j - 1] - per
                t[n - k + j - 1] = t[k + j] + per
        yield xp.asarray(t)

        # construct the LSQ spline with this set of knots
        fpold = fp
        residuals, fp = _get_residuals(x, y, t, k, w=w,
                                        periodic=periodic)
        fpms = fp - s

        # c  test whether the approximation sinf(x) is an acceptable solution.
        # c  if f(p=inf) < s accept the choice of knots.
        if (abs(fpms) < acc) or (fpms < 0):
            return

        # ### c  increase the number of knots. ###

        # c  determine the number of knots nplus we are going to add.
        if n == nmin:
            # the first iteration
            nplus = 1
        else:
            delta = fpold - fp
            npl1 = int(nplus * fpms / delta) if delta > acc else nplus*2
            nplus = min(nplus*2, max(npl1, nplus//2, 1))

        # actually add knots
        for j in range(nplus):
            t = add_knot(x, t, k, residuals)

            # check if we have enough knots already

            n = t.shape[0]
            # c  if n = nmax, sinf(x) is an interpolating spline.
            # c  if n=nmax we locate the knots as for interpolation.
            if n >= nmax:
                if not periodic:
                    t = _not_a_knot(x, k)
                else:
                    t = _periodic_knots(x, k)
                yield xp.asarray(t)
                return

            # c  if n=nest we cannot increase the number of knots because of
            # c  the storage capacity limitation.
            if n >= nest:
                yield xp.asarray(t)
                return

            # recompute if needed
            if j < nplus - 1:
                residuals, _ = _get_residuals(x, y, t, k, w=w, periodic=periodic)

    # this should never be reached
    return


#   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#   c  part 2: determination of the smoothing spline sp(x).                c
#   c  ***************************************************                 c
#   c  we have determined the number of knots and their position.          c
#   c  we now compute the b-spline coefficients of the smoothing spline    c
#   c  sp(x). the observation matrix a is extended by the rows of matrix   c
#   c  b expressing that the kth derivative discontinuities of sp(x) at    c
#   c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
#   c  ponding weights of these additional rows are set to 1/p.            c
#   c  iteratively we then have to determine the value of p such that      c
#   c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
#   c  the least-squares kth degree polynomial corresponds to p=0, and     c
#   c  that the least-squares spline corresponds to p=infinity. the        c
#   c  iteration process which is proposed here, makes use of rational     c
#   c  interpolation. since f(p) is a convex and strictly decreasing       c
#   c  function of p, it can be approximated by a rational function        c
#   c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
#   c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
#   c  to calculate the new value of p such that r(p)=s. convergence is    c
#   c  guaranteed by taking f1>0 and f3<0.                                 c
#   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


def prodd(t, i, j, k):
    res = 1.0
    for s in range(k+2):
        if i + s != j:
            res *= (t[j] - t[i+s])
    return res


def disc(t, k):
    """Discontinuity matrix: jumps of k-th derivatives of b-splines at internal knots.

    See Eqs. (9)-(10) of Ref. [1], or, equivalently, Eq. (3.43) of Ref. [2].

    This routine assumes internal knots are all simple (have multiplicity =1).

    Parameters
    ----------
    t : ndarray, 1D, shape(n,)
        Knots.
    k : int
        The spline degree

    Returns
    -------
    disc : ndarray, shape(n-2*k-1, k+2)
        The jumps of the k-th derivatives of b-splines at internal knots,
        ``t[k+1], ...., t[n-k-1]``.
    offset : ndarray, shape(2-2*k-1,)
        Offsets
    nc : int

    Notes
    -----

    The normalization here follows FITPACK:
    (https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpdisc.f#L36)

    The k-th derivative jumps are multiplied by a factor::

        (delta / nrint)**k / k!

    where ``delta`` is the length of the interval spanned by internal knots, and
    ``nrint`` is one less the number of internal knots (i.e., the number of
    subintervals between them).

    References
    ----------
    .. [1] Paul Dierckx, Algorithms for smoothing data with periodic and parametric
           splines, Computer Graphics and Image Processing, vol. 20, p. 171 (1982).
           :doi:`10.1016/0146-664X(82)90043-0`

    .. [2] Tom Lyche and Knut Morken, Spline methods,
        http://www.uio.no/studier/emner/matnat/ifi/INF-MAT5340/v05/undervisningsmateriale/

    """
    n = t.shape[0]

    # the length of the base interval spanned by internal knots & the number
    # of subintervas between these internal knots
    delta = t[n - k - 1] - t[k]
    nrint = n - 2*k - 1

    matr = np.empty((nrint - 1, k + 2), dtype=float)
    for jj in range(nrint - 1):
        j = jj + k + 1
        for ii in range(k + 2):
            i = jj + ii
            matr[jj, ii] = (t[i + k + 1] - t[i]) / prodd(t, i, j, k)
        # NB: equivalent to
        # row = [(t[i + k + 1] - t[i]) / prodd(t, i, j, k) for i in range(j-k-1, j+1)]
        # assert (matr[j-k-1, :] == row).all()

    # follow FITPACK
    matr *= (delta/ nrint)**k

    # make it packed
    offset = np.array([i for i in range(nrint-1)], dtype=np.int64)
    nc = n - k - 1
    return matr, offset, nc


class F:
    """ The r.h.s. of ``f(p) = s``.

    Given scalar `p`, we solve the system of equations in the LSQ sense:

        | A     |  @ | c | = | y |
        | B / p |    | 0 |   | 0 |

    where `A` is the matrix of b-splines and `b` is the discontinuity matrix
    (the jumps of the k-th derivatives of b-spline basis elements at knots).

    Since we do that repeatedly while minimizing over `p`, we QR-factorize
    `A` only once and update the QR factorization only of the `B` rows of the
    augmented matrix |A, B/p|.

    The system of equations is Eq. (15) Ref. [1]_, the strategy and implementation
    follows that of FITPACK, see specific links below.

    References
    ----------
    [1] P. Dierckx, Algorithms for Smoothing Data with Periodic and Parametric Splines,
        COMPUTER GRAPHICS AND IMAGE PROCESSING vol. 20, pp 171-184 (1982.)
        https://doi.org/10.1016/0146-664X(82)90043-0

    """
    def __init__(self, x, y, t, k, s, w=None, *, R=None, Y=None):
        self.x = x
        self.y = y
        self.t = t
        self.k = k
        w = np.ones_like(x, dtype=float) if w is None else w
        if w.ndim != 1:
            raise ValueError(f"{w.ndim = } != 1.")
        self.w = w
        self.s = s

        if y.ndim != 2:
            raise ValueError(f"F: expected y.ndim == 2, got {y.ndim = } instead.")

        # ### precompute what we can ###

        # https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L250
        # c  evaluate the discontinuity jump of the kth derivative of the
        # c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
        b, b_offset, b_nc = disc(t, k)

        # the QR factorization of the data matrix, if not provided
        # NB: otherwise, must be consistent with x,y & s, but this is not checked
        if R is None and Y is None:
            R, Y, _, _, _ = _lsq_solve_qr(x, y, t, k, w)

        # prepare to combine R and the discontinuity matrix (AB); also r.h.s. (YY)
        # https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L269
        # c  the rows of matrix b with weight 1/p are rotated into the
        # c  triangularised observation matrix a which is stored in g.
        nc = t.shape[0] - k - 1
        nz = k + 1
        if R.shape[1] != nz:
            raise ValueError(f"Internal error: {R.shape[1] =} != {k+1 =}.")

        # r.h.s. of the augmented system
        z = np.zeros((b.shape[0], Y.shape[1]), dtype=float)
        self.YY = np.r_[Y[:nc], z]

        # l.h.s. of the augmented system
        AA = np.zeros((nc + b.shape[0], self.k+2), dtype=float)
        AA[:nc, :nz] = R[:nc, :]
        # AA[nc:, :] = b.a / p  # done in __call__(self, p)
        self.AA  = AA
        self.offset = np.r_[np.arange(nc, dtype=np.int64), b_offset]

        self.nc = nc
        self.b = b

    def __call__(self, p):
        # https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L279
        # c  the row of matrix b is rotated into triangle by givens transformation

        # copy the precomputed matrices over for in-place work
        # R = PackedMatrix(self.AB.a.copy(), self.AB.offset.copy(), nc)
        AB = self.AA.copy()
        offset = self.offset.copy()
        nc = self.nc

        AB[nc:, :] = self.b / p
        QY = self.YY.copy()

        # heavy lifting happens here, in-place
        _dierckx.qr_reduce(AB, offset, nc, QY, startrow=nc)

        # solve for the coefficients
        c, _, fp = _dierckx.fpback(AB, nc, self.x, self.y, self.t, self.k, self.w, QY)

        spl = BSpline(self.t, c, self.k)
        self.spl = spl   # store it

        return fp - self.s

class Fperiodic:
    """
    Fit a smooth periodic B-spline curve to given data points.

    This class fits a periodic B-spline curve S(t) of degree k through data points
    (x, y) with knots t. The spline is smooth and repeats itself at the start and
    end, meaning the function and its derivatives up to order k-1 are equal
    at the boundaries.

    We want to find spline coefficients c that minimize the difference between
    the spline and the data, while also keeping the spline smooth. This is done
    by solving:

        minimize || W^{1/2} (Y - B c) ||^2 + s * c^T @ R @ c
        subject to periodic constraints on c.

    where:
      - Y is the data values,
      - B is the matrix of B-spline basis functions at points x,
      - W is a weighting matrix for the data points,
      - s is the smoothing parameter (larger s means smoother curve),
      - R is a matrix that penalizes wiggliness of the spline,
      - c spline coefficients to be solved for
      - periodic constraints ensure the spline repeats smoothly.

    The solution is obtained by forming augmented matrices and performing
    a QR factorization that incorporates these constraints, following the
    approach in FITPACK's `fpperi.f`.

    Parameters:
    -----------
    x : array_like, shape (n,)
    y : array_like, shape (n, m)
    t : array_like, shape (nt,)
        Knot vector for the spline
    k : int
        Degree of the spline.
    s : float
        Controls smoothness: bigger s means smoother curve, smaller s fits data closer.
    w : array_like, shape (n,), optional
        Weights for data points. Defaults to all ones.
    R, Y, A1, A2, Z : arrays, optional
        Precomputed matrices from least squares and QR factorization steps to speed up
        repeated fits with the same knots and data.

    Attributes:
    -----------
    G1, G2 : arrays
        Augmented matrices combining the original QR factors and constraints related to
        the spline basis and data. G1 is roughly the "upper-triangular" part;
        G2 contains additional constraint information for periodicity.

    H1, H2 : arrays
        Matrices associated with the discontinuity jump constraints of the k-th
        derivative of B-splines at the knots. These encode the periodicity
        conditions and are scaled by the smoothing parameter.

    offset : array
        Offset indices used for efficient indexing during QR reduction.

    Methods:
    --------
    __call__(p):
        Perform QR reduction of augmented matrices scaled by 1/p, solve for spline
        coefficients, and return residual difference fp - s.

    References:
    -----------
    - FITPACK's fpperi.f and fpcurf.f Fortran routines for periodic spline fitting.
    """
    def __init__(self, x, y, t, k, s, w=None, *,
                 R=None, Y=None, A1=None, A2=None, Z=None):
        # Initialize the class with input data points x, y,
        # knot vector t, spline degree k, smoothing factor s,
        # optional weights w, and optionally precomputed matrices.
        self.x = x
        self.y = y
        self.t = t
        self.k = k

        # If weights not provided, default to uniform weights (all ones)
        w = np.ones_like(x, dtype=float) if w is None else w

        if w.ndim != 1:
            raise ValueError(f"{w.ndim = } != 1.")
        self.w = w

        self.s = s

        if y.ndim != 2:
            raise ValueError(f"F: expected y.ndim == 2, got {y.ndim = } instead.")

        # ### precompute matrices and factors needed for spline fitting ###

        # Compute the discontinuity jump vector 'b' for the k-th derivative
        # of B-splines at internal knots. This is needed for enforcing
        # periodicity constraints. The jump discontinuities are stored in b.
        # Refer: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpcurf.f#L252
        b, b_offset, b_nc = disc(t, k)

        # If QR factorization or auxiliary matrices (A1, A2, Z) are not provided,
        # compute them via least squares on the data (x,y) with weights w.
        # These matrices come from fitting B-spline basis to data.
        if ((A1 is None or A2 is None or Z is None) or
            (R is None and Y is None)):
            # _lsq_solve_qr_for_root_rati_periodic computes QR factorization of
            # data matrix and returns related matrices.
            # Refer: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L171-L215
            R, A1, A2, Z, _, _, _, _ = _lsq_solve_qr_for_root_rati_periodic(
                x, y, t, k, w)

        # Initialize augmented matrices G1, G2, H1, H2 and offset used
        # for periodic B-spline fitting with constraints.
        # This calls the C++ function `init_augmented_matrices` to set
        # up these matrices based on A1, A2, and discontinuity vector b.
        # Refer: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L441-L493
        G1, G2, H1, H2, offset = _dierckx.init_augmented_matrices(A1, A2, b, len(t), k)

        # Store these matrices as class attributes for use in evaluation
        self.G1_ = G1
        self.G2_ = G2
        self.H1_ = H1
        self.H2_ = H2
        self.Z_ = Z
        self.offset_ = offset

    def __call__(self, p):
        # Create copies of the augmented matrices to avoid overwriting originals
        G1 = self.G1_.copy()
        G2 = self.G2_.copy()
        H1 = self.H1_.copy()
        H2 = self.H2_.copy()
        Z = self.Z_.copy()

        # Scale H1 and H2 by the inverse of p (1/p), applying the pinv multiplier
        # This scaling is required before QR reduction step to control smoothing.
        pinv = 1/p
        H1 = H1 * pinv
        H2 = H2 * pinv

        # Initialize vector c for coefficients with shape compatible with matrices.
        # The first part of c is set from Z, which is related to least squares fit.
        c = np.empty((len(self.t) - self.k - 1, Z.shape[1]), dtype=Z.dtype)
        c[:len(self.t) - 2*self.k - 1, :] = Z[:len(self.t) - 2*self.k - 1, :]

        # Perform QR factorization reduction on the augmented matrices
        # This corresponds to the C++ function qr_reduce_augmented_matrices.
        # It applies Givens rotations to eliminate entries and
        # reduces the problem dimension by accounting for constraints.
        # Refer: https://github.com/scipy/scipy/blob/maintenance/1.16.x/scipy/interpolate/fitpack/fpperi.f#L498-L534
        _dierckx.qr_reduce_augmented_matrices(
            G1, G2, H1, H2, c, self.offset_, len(self.t), self.k)

        # Solve for the B-spline coefficients by backward substitution
        # using the reduced matrices. This corresponds to the fpbacp routine in fitpack.
        c, _, fp = _dierckx.fpbacp(
            G1, G2, c,
            self.k, self.k + 1, self.x[:-1], self.y[:-1, :],
            self.t, self.w[:-1])

        # Construct a BSpline object using knot vector t, coefficients c, and degree k
        spl = BSpline(self.t, c, self.k)
        self.spl = spl  # store it in the object

        # Return the difference between the final residual fp and smoothing parameter s
        # This quantity is used to assess fit quality in root-finding procedures.
        return fp - self.s


def fprati(p1, f1, p2, f2, p3, f3):
    """The root of r(p) = (u*p + v) / (p + w) given three points and values,
    (p1, f2), (p2, f2) and (p3, f3).

    The FITPACK analog adjusts the bounds, and we do not
    https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fprati.f

    NB: FITPACK uses p < 0 to encode p=infinity. We just use the infinity itself.
    Since the bracket is ``p1 <= p2 <= p3``, ``p3`` can be infinite (in fact,
    this is what the minimizer starts with, ``p3=inf``).
    """
    h1 = f1 * (f2 - f3)
    h2 = f2 * (f3 - f1)
    h3 = f3 * (f1 - f2)
    if p3 == np.inf:
        return -(p2*h1 + p1*h2) / h3
    return -(p1*p2*h3 + p2*p3*h1 + p1*p3*h2) / (p1*h1 + p2*h2 + p3*h3)


class Bunch:
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)


_iermesg1 = """error. a theoretically impossible result was found during
the iteration process for finding a smoothing spline with
fp = s. probably causes : s too small.
"""

_iermesg = {
1: _iermesg1 + """the weighted sum of squared residuals is becoming NaN
""",
2: _iermesg1 + """there is an approximation returned but the corresponding
weighted sum of squared residuals does not satisfy the
condition abs(fp-s)/s < tol.
""",
3: """error. the maximal number of iterations maxit (set to 20
by the program) allowed for finding a smoothing spline
with fp=s has been reached. probably causes : s too small
there is an approximation returned but the corresponding
weighted sum of squared residuals does not satisfy the
condition abs(fp-s)/s < tol.
"""
}


def root_rati(f, p0, bracket, acc):
    """Solve `f(p) = 0` using a rational function approximation.

    In a nutshell, since the function f(p) is known to be monotonically decreasing, we
       - maintain the bracket (p1, f1), (p2, f2) and (p3, f3)
       - at each iteration step, approximate f(p) by a rational function
         r(p) = (u*p + v) / (p + w)
         and make a step to p_new to the root of f(p): r(p_new) = 0.
         The coefficients u, v and w are found from the bracket values p1..3 and f1...3

    The algorithm and implementation follows
    https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L229
    and
    https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fppara.f#L290

    Note that the latter is for parametric splines and the former is for 1D spline
    functions. The minimization is indentical though [modulo a summation over the
    dimensions in the computation of f(p)], so we reuse the minimizer for both
    d=1 and d>1.
    """
    # Magic values from
    # https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L27
    con1 = 0.1
    con9 = 0.9
    con4 = 0.04

    # bracketing flags (follow FITPACK)
    # https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fppara.f#L365
    ich1, ich3 = 0, 0

    (p1, f1), (p3, f3)  = bracket
    p = p0

    for it in range(MAXIT):
        p2, f2 = p, f(p)

        # c  test whether the approximation sp(x) is an acceptable solution.
        if abs(f2) < acc:
            ier, converged = 0, True
            break

        # c  carry out one more step of the iteration process.
        if ich3 == 0:
            if f2 - f3 <= acc:
                # c  our initial choice of p is too large.
                p3 = p2
                f3 = f2
                p = p*con4
                if p <= p1:
                     p = p1*con9 + p2*con1
                continue
            else:
                if f2 < 0:
                    ich3 = 1

        if ich1 == 0:
            if f1 - f2 <= acc:
                # c  our initial choice of p is too small
                p1 = p2
                f1 = f2
                p = p/con4
                if p3 != np.inf and p <= p3:
                     p = p2*con1 + p3*con9
                continue
            else:
                if f2 > 0:
                    ich1 = 1

        # c  test whether the iteration process proceeds as theoretically expected.
        # [f(p) should be monotonically decreasing]
        if f1 <= f2 or f2 <= f3:
            ier, converged = 2, False
            break

        # actually make the iteration step
        p = fprati(p1, f1, p2, f2, p3, f3)

        # c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
        if f2 < 0:
            p3, f3 = p2, f2
        else:
            p1, f1 = p2, f2

    else:
        # not converged in MAXIT iterations
        ier, converged = 3, False

    if ier != 0:
        warnings.warn(RuntimeWarning(_iermesg[ier]), stacklevel=2)

    return Bunch(converged=converged, root=p, iterations=it, ier=ier)

def _make_splrep_impl(x, y, w, xb, xe, k, s, t, nest, periodic, xp=np):
    """Shared infra for make_splrep and make_splprep.
    """
    acc = s * TOL
    m = x.size    # the number of data points

    if nest is None:
        # the max number of knots. This is set in _fitpack_impl.py line 274
        # and fitpack.pyf line 198
        # Ref: https://github.com/scipy/scipy/blob/596b586e25e34bd842b575bac134b4d6924c6556/scipy/interpolate/_fitpack_impl.py#L260-L263
        if periodic:
            nest = max(m + 2*k, 2*k + 3)
        else:
            nest = max(m + k + 1, 2*k + 3)
    else:
        if nest < 2*(k + 1):
            raise ValueError(f"`nest` too small: {nest = } < 2*(k+1) = {2*(k+1)}.")
        if t is not None:
            raise ValueError("Either supply `t` or `nest`.")

    if t is None:
        gen = _generate_knots_impl(x, y, w, xb, xe, k, s, nest, periodic)
        t = list(gen)[-1]
    else:
        fpcheck(x, t, k, periodic=periodic)

    if t.shape[0] == 2 * (k + 1):
        # nothing to optimize
        _, _, c, _, _ = _lsq_solve_qr(x, y, t, k, w, periodic=periodic)
        t, c = xp.asarray(t), xp.asarray(c)
        return BSpline(t, c, k)

    ### solve ###

    # c  initial value for p.
    # https://github.com/scipy/scipy/blob/maintenance/1.11.x/scipy/interpolate/fitpack/fpcurf.f#L253
    if periodic:
        # N.B. - Check _lsq_solve_qr computation
        # of p for periodic splines
        R, A1, A2, Z, Y, _, p, _ = _lsq_solve_qr_for_root_rati_periodic(x, y, t, k, w)
    else:
        R, Y, _, _, _ = _lsq_solve_qr(x, y, t, k, w, periodic=periodic)
    nc = t.shape[0] -k -1
    if not periodic:
        p = nc / R[:, 0].sum()

    # ### bespoke solver ####
    # initial conditions
    # f(p=inf) : LSQ spline with knots t   (XXX: reuse R, c)
    # N.B. - Check _lsq_solve_qr which is called
    # via _get_residuals for logic behind
    # computation of fp for periodic splines
    _, fp = _get_residuals(x, y, t, k, w, periodic=periodic)
    fpinf = fp - s

    # f(p=0): LSQ spline without internal knots
    if not periodic:
        _, fp0 = _get_residuals(x, y, np.array([xb]*(k+1) + [xe]*(k+1)), k, w)
        fp0 = fp0 - s
    else:
        # f(p=0) is fp for constant function
        # in case of periodic splines
        per = xe - xb
        tc = np.zeros(2*(k + 1), dtype=float)
        for i in range(0, k + 1):
            tc[i] = x[0] - (k - i) * per
            tc[i + k + 1] = x[m - 1] + i * per
        _, fp0 = _get_residuals(x, y, tc, k, w, periodic=periodic)
        fp0 = fp0 - s

    # solve
    bracket = (0, fp0), (np.inf, fpinf)
    if not periodic:
        f = F(x, y, t, k=k, s=s, w=w, R=R, Y=Y)
    else:
        f = Fperiodic(x, y, t, k=k, s=s, w=w, R=R, Y=Y, A1=A1, A2=A2, Z=Z)
    _ = root_rati(f, p, bracket, acc)

    # solve ALTERNATIVE: is roughly equivalent, gives slightly different results
    # starting from scratch, that would have probably been tolerable;
    # backwards compatibility dictates that we replicate the FITPACK minimizer though.
 #   f = F(x, y, t, k=k, s=s, w=w, R=R, Y=Y)
 #   from scipy.optimize import root_scalar
 #   res_ = root_scalar(f, x0=p, rtol=acc)
 #   assert res_.converged

    # f.spl is the spline corresponding to the found `p` value
    t, c, k = f.spl.tck
    axis, extrap = f.spl.axis, f.spl.extrapolate
    t, c = xp.asarray(t), xp.asarray(c)
    spl = BSpline.construct_fast(t, c, k, axis=axis, extrapolate=extrap)
    return spl


@xp_capabilities(cpu_only=True, jax_jit=False, allow_dask_compute=True)
def make_splrep(x, y, *, w=None, xb=None, xe=None,
                k=3, s=0, t=None, nest=None, bc_type=None):
    r"""Create a smoothing B-spline function with bounded error, minimizing derivative jumps.

    Given the set of data points ``(x[i], y[i])``, determine a smooth spline
    approximation of degree ``k`` on the interval ``xb <= x <= xe``.

    Parameters
    ----------
    x, y : array_like, shape (m,)
        The data points defining a curve ``y = f(x)``.
    w : array_like, shape (m,), optional
        Strictly positive 1D array of weights, of the same length as `x` and `y`.
        The weights are used in computing the weighted least-squares spline
        fit. If the errors in the y values have standard-deviation given by the
        vector ``d``, then `w` should be ``1/d``.
        Default is ``np.ones(m)``.
    xb, xe : float, optional
        The interval to fit.  If None, these default to ``x[0]`` and ``x[-1]``,
        respectively.
    k : int, optional
        The degree of the spline fit. It is recommended to use cubic splines,
        ``k=3``, which is the default. Even values of `k` should be avoided,
        especially with small `s` values.
    s : float, optional
        The smoothing condition. The amount of smoothness is determined by
        satisfying the LSQ (least-squares) constraint::

            sum((w * (g(x)  - y))**2 ) <= s

        where ``g(x)`` is the smoothed fit to ``(x, y)``. The user can use `s`
        to control the tradeoff between closeness to data and smoothness of fit.
        Larger `s` means more smoothing while smaller values of `s` indicate less
        smoothing.
        Recommended values of `s` depend on the weights, `w`. If the weights
        represent the inverse of the standard deviation of `y`, then a good `s`
        value should be found in the range ``(m-sqrt(2*m), m+sqrt(2*m))`` where
        ``m`` is the number of datapoints in `x`, `y`, and `w`.
        Default is ``s = 0.0``, i.e. interpolation.
    t : array_like, optional
        The spline knots. If None (default), the knots will be constructed
        automatically.
        There must be at least ``2*k + 2`` and at most ``m + k + 1`` knots.
    nest : int, optional
        The target length of the knot vector. Should be between ``2*(k + 1)``
        (the minimum number of knots for a degree-``k`` spline), and
        ``m + k + 1`` (the number of knots of the interpolating spline).
        The actual number of knots returned by this routine may be slightly
        larger than `nest`.
        Default is None (no limit, add up to ``m + k + 1`` knots).
    periodic : bool, optional
        If True, data points are considered periodic with period ``x[m-1]`` -
        ``x[0]`` and a smooth periodic spline approximation is returned. Values of
        ``y[m-1]`` and ``w[m-1]`` are not used.
        The default is False, corresponding to boundary condition 'not-a-knot'.
    bc_type : str, optional
        Boundary conditions.
        Default is `"not-a-knot"`.
        The following boundary conditions are recognized:

        * ``"not-a-knot"`` (default): The first and second segments are the
          same polynomial. This is equivalent to having ``bc_type=None``.
        * ``"periodic"``: The values and the first ``k-1`` derivatives at the
          ends are equivalent.

    Returns
    -------
    spl : a `BSpline` instance
        For `s=0`,  ``spl(x) == y``.
        For non-zero values of `s` the `spl` represents the smoothed approximation
        to `(x, y)`, generally with fewer knots.

    See Also
    --------
    generate_knots : is used under the hood for generating the knots
    make_splprep : the analog of this routine for parametric curves
    make_interp_spline : construct an interpolating spline (``s = 0``)
    make_lsq_spline : construct the least-squares spline given the knot vector
    splrep : a FITPACK analog of this routine

    References
    ----------
    .. [1] P. Dierckx, "Algorithms for smoothing data with periodic and
        parametric splines, Computer Graphics and Image Processing",
        20 (1982) 171-184.
    .. [2] P. Dierckx, "Curve and surface fitting with splines", Monographs on
        Numerical Analysis, Oxford University Press, 1993.

    Notes
    -----
    This routine constructs the smoothing spline function, :math:`g(x)`, to
    minimize the sum of jumps, :math:`D_j`, of the ``k``-th derivative at the
    internal knots (:math:`x_b < t_i < x_e`), where

    .. math::

        D_i = g^{(k)}(t_i + 0) - g^{(k)}(t_i - 0)

    Specifically, the routine constructs the spline function :math:`g(x)` which
    minimizes

    .. math::

            \sum_i | D_i |^2 \to \mathrm{min}

    provided that

    .. math::

           \sum_{j=1}^m (w_j \times (g(x_j) - y_j))^2 \leqslant s ,

    where :math:`s > 0` is the input parameter.

    In other words, we balance maximizing the smoothness (measured as the jumps
    of the derivative, the first criterion), and the deviation of :math:`g(x_j)`
    from the data :math:`y_j` (the second criterion).

    Note that the summation in the second criterion is over all data points,
    and in the first criterion it is over the internal spline knots (i.e.
    those with ``xb < t[i] < xe``). The spline knots are in general a subset
    of data, see `generate_knots` for details.

    Also note the difference of this routine to `make_lsq_spline`: the latter
    routine does not consider smoothness and simply solves a least-squares
    problem

    .. math::

        \sum w_j \times (g(x_j) - y_j)^2 \to \mathrm{min}

    for a spline function :math:`g(x)` with a _fixed_ knot vector ``t``.

    .. versionadded:: 1.15.0
    """  # noqa:E501
    xp = array_namespace(x, y, w, t)
    if t is not None:
        t = xp.asarray(t)
    bc_type = _validate_bc_type(bc_type)

    if s == 0:
        if t is not None or w is not None or nest is not None:
            raise ValueError("s==0 is for interpolation only")
        return make_interp_spline(x, y, k=k, bc_type=bc_type)

    periodic = (bc_type == 'periodic')

    x, y, w, k, s, xb, xe = _validate_inputs(x, y, w, k, s, xb, xe,
                                             parametric=False, periodic=periodic)

    spl = _make_splrep_impl(x, y, w, xb, xe, k, s, t, nest, periodic, xp=xp)

    # postprocess: squeeze out the last dimension: was added to simplify the internals.
    spl.c = spl.c[:, 0]
    return spl


@xp_capabilities(cpu_only=True, jax_jit=False, allow_dask_compute=True)
def make_splprep(x, *, w=None, u=None, ub=None, ue=None,
                 k=3, s=0, t=None, nest=None, bc_type=None):
    r"""
    Create a smoothing parametric B-spline curve with bounded error, minimizing derivative jumps.

    Given a list of N 1D arrays, `x`, which represent a curve in
    N-dimensional space parametrized by `u`, find a smooth approximating
    spline curve ``g(u)``.

    Parameters
    ----------
    x : array_like, shape (ndim, m)
        Sampled data points representing the curve in ``ndim`` dimensions.
        The typical use is a list of 1D arrays, each of length ``m``.
    w : array_like, shape(m,), optional
        Strictly positive 1D array of weights.
        The weights are used in computing the weighted least-squares spline
        fit. If the errors in the `x` values have standard deviation given by
        the vector d, then `w` should be 1/d. Default is ``np.ones(m)``.
    u : array_like, optional
        An array of parameter values for the curve in the parametric form.
        If not given, these values are calculated automatically, according to::

            v[0] = 0
            v[i] = v[i-1] + distance(x[i], x[i-1])
            u[i] = v[i] / v[-1]

    ub, ue : float, optional
        The end-points of the parameters interval. Default to ``u[0]`` and ``u[-1]``.
    k : int, optional
        Degree of the spline. Cubic splines, ``k=3``, are recommended.
        Even values of `k` should be avoided especially with a small ``s`` value.
        Default is ``k=3``
    s : float, optional
        A smoothing condition.  The amount of smoothness is determined by
        satisfying the conditions::

            sum((w * (g(u) - x))**2) <= s,

        where ``g(u)`` is the smoothed approximation to ``x``.  The user can
        use `s` to control the trade-off between closeness and smoothness
        of fit.  Larger ``s`` means more smoothing while smaller values of ``s``
        indicate less smoothing.
        Recommended values of ``s`` depend on the weights, ``w``.  If the weights
        represent the inverse of the standard deviation of ``x``, then a good
        ``s`` value should be found in the range ``(m - sqrt(2*m), m + sqrt(2*m))``,
        where ``m`` is the number of data points in ``x`` and ``w``.
    t : array_like, optional
        The spline knots. If None (default), the knots will be constructed
        automatically.
        There must be at least ``2*k + 2`` and at most ``m + k + 1`` knots.
    nest : int, optional
        The target length of the knot vector. Should be between ``2*(k + 1)``
        (the minimum number of knots for a degree-``k`` spline), and
        ``m + k + 1`` (the number of knots of the interpolating spline).
        The actual number of knots returned by this routine may be slightly
        larger than `nest`.
        Default is None (no limit, add up to ``m + k + 1`` knots).
    bc_type : str, optional
        Boundary conditions.
        Default is `"not-a-knot"`.
        The following boundary conditions are recognized:

        * ``"not-a-knot"`` (default): The first and second segments are the
          same polynomial. This is equivalent to having ``bc_type=None``.
        * ``"periodic"``: The values and the first ``k-1`` derivatives at the
          ends are equivalent.

    Returns
    -------
    spl : a `BSpline` instance
        For `s=0`,  ``spl(u) == x``.
        For non-zero values of ``s``, `spl` represents the smoothed approximation
        to ``x``, generally with fewer knots.
    u : ndarray
        The values of the parameters

    See Also
    --------
    generate_knots : is used under the hood for generating the knots
    make_splrep : the analog of this routine 1D functions
    make_interp_spline : construct an interpolating spline (``s = 0``)
    make_lsq_spline : construct the least-squares spline given the knot vector
    splprep : a FITPACK analog of this routine

    Notes
    -----
    Given a set of :math:`m` data points in :math:`D` dimensions, :math:`\vec{x}_j`,
    with :math:`j=1, ..., m` and :math:`\vec{x}_j = (x_{j; 1}, ..., x_{j; D})`,
    this routine constructs the parametric spline curve :math:`g_a(u)` with
    :math:`a=1, ..., D`, to minimize the sum of jumps, :math:`D_{i; a}`, of the
    ``k``-th derivative at the internal knots (:math:`u_b < t_i < u_e`), where

    .. math::

        D_{i; a} = g_a^{(k)}(t_i + 0) - g_a^{(k)}(t_i - 0)

    Specifically, the routine constructs the spline function :math:`g(u)` which
    minimizes

    .. math::

            \sum_i \sum_{a=1}^D | D_{i; a} |^2 \to \mathrm{min}

    provided that

    .. math::

        \sum_{j=1}^m \sum_{a=1}^D (w_j \times (g_a(u_j) - x_{j; a}))^2 \leqslant s

    where :math:`u_j` is the value of the parameter corresponding to the data point
    :math:`(x_{j; 1}, ..., x_{j; D})`, and :math:`s > 0` is the input parameter.

    In other words, we balance maximizing the smoothness (measured as the jumps
    of the derivative, the first criterion), and the deviation of :math:`g(u_j)`
    from the data :math:`x_j` (the second criterion).

    Note that the summation in the second criterion is over all data points,
    and in the first criterion it is over the internal spline knots (i.e.
    those with ``ub < t[i] < ue``). The spline knots are in general a subset
    of data, see `generate_knots` for details.

    .. versionadded:: 1.15.0

    References
    ----------
    .. [1] P. Dierckx, "Algorithms for smoothing data with periodic and
        parametric splines, Computer Graphics and Image Processing",
        20 (1982) 171-184.
    .. [2] P. Dierckx, "Curve and surface fitting with splines", Monographs on
        Numerical Analysis, Oxford University Press, 1993.
    """  # noqa:E501

    bc_type = _validate_bc_type(bc_type)

    periodic = (bc_type == 'periodic')

    # x can be either a 2D array or a list of 1D arrays
    if isinstance(x, list):
        xp = array_namespace(*x, w, u, t)
    else:
        xp = array_namespace(x, w, u, t)

    x = xp.stack(x, axis=1)

    # construct the default parametrization of the curve
    if u is None:
        dp = (x[1:, :] - x[:-1, :])**2
        u = xp.cumulative_sum( xp.sqrt(xp.sum(dp, axis=1)) )
        u = concat_1d(xp, 0., u / u[-1])

    if s == 0:
        if t is not None or w is not None or nest is not None:
            raise ValueError("s==0 is for interpolation only")
        return make_interp_spline(u, x.T, k=k, axis=1), u

    u, x, w, k, s, ub, ue = _validate_inputs(u, x, w, k, s, ub, ue,
                                             parametric=True, periodic=periodic)

    if t is not None:
        t = xp.asarray(t)

    spl = _make_splrep_impl(u, x, w, ub, ue, k, s, t,
                            nest, periodic=periodic, xp=xp)

    # posprocess: `axis=1` so that spl(u).shape == np.shape(x)
    # when `x` is a list of 1D arrays (cf original splPrep)
    cc = spl.c.T
    spl1 = BSpline(spl.t, cc, spl.k, axis=1)

    return spl1, xp.asarray(u)
