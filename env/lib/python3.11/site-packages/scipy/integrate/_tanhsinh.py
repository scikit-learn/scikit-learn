# mypy: disable-error-code="attr-defined"
import numpy as np
from scipy import special
import scipy._lib._elementwise_iterative_method as eim
from scipy._lib._util import _RichResult

# todo:
#  figure out warning situation
#  address https://github.com/scipy/scipy/pull/18650#discussion_r1233032521
#  without `minweight`, we are also suppressing infinities within the interval.
#    Is that OK? If so, we can probably get rid of `status=3`.
#  Add heuristic to stop when improvement is too slow / antithrashing
#  support singularities? interval subdivision? this feature will be added
#    eventually, but do we adjust the interface now?
#  When doing log-integration, should the tolerances control the error of the
#    log-integral or the error of the integral?  The trouble is that `log`
#    inherently looses some precision so it may not be possible to refine
#    the integral further. Example: 7th moment of stats.f(15, 20)
#  respect function evaluation limit?
#  make public?


def _tanhsinh(f, a, b, *, args=(), log=False, maxfun=None, maxlevel=None,
              minlevel=2, atol=None, rtol=None, preserve_shape=False,
              callback=None):
    """Evaluate a convergent integral numerically using tanh-sinh quadrature.

    In practice, tanh-sinh quadrature achieves quadratic convergence for
    many integrands: the number of accurate *digits* scales roughly linearly
    with the number of function evaluations [1]_.

    Either or both of the limits of integration may be infinite, and
    singularities at the endpoints are acceptable. Divergent integrals and
    integrands with non-finite derivatives or singularities within an interval
    are out of scope, but the latter may be evaluated be calling `_tanhsinh` on
    each sub-interval separately.

    Parameters
    ----------
    f : callable
        The function to be integrated. The signature must be::
            func(x: ndarray, *fargs) -> ndarray
         where each element of ``x`` is a finite real and ``fargs`` is a tuple,
         which may contain an arbitrary number of arrays that are broadcastable
         with `x`. ``func`` must be an elementwise-scalar function; see
         documentation of parameter `preserve_shape` for details.
         If ``func`` returns a value with complex dtype when evaluated at
         either endpoint, subsequent arguments ``x`` will have complex dtype
         (but zero imaginary part).
    a, b : array_like
        Real lower and upper limits of integration. Must be broadcastable.
        Elements may be infinite.
    args : tuple, optional
        Additional positional arguments to be passed to `func`. Must be arrays
        broadcastable with `a` and `b`. If the callable to be integrated
        requires arguments that are not broadcastable with `a` and `b`, wrap
        that callable with `f`. See Examples.
    log : bool, default: False
        Setting to True indicates that `f` returns the log of the integrand
        and that `atol` and `rtol` are expressed as the logs of the absolute
        and relative errors. In this case, the result object will contain the
        log of the integral and error. This is useful for integrands for which
        numerical underflow or overflow would lead to inaccuracies.
        When ``log=True``, the integrand (the exponential of `f`) must be real,
        but it may be negative, in which case the log of the integrand is a
        complex number with an imaginary part that is an odd multiple of π.
    maxlevel : int, default: 10
        The maximum refinement level of the algorithm.

        At the zeroth level, `f` is called once, performing 16 function
        evaluations. At each subsequent level, `f` is called once more,
        approximately doubling the number of function evaluations that have
        been performed. Accordingly, for many integrands, each successive level
        will double the number of accurate digits in the result (up to the
        limits of floating point precision).

        The algorithm will terminate after completing level `maxlevel` or after
        another termination condition is satisfied, whichever comes first.
    minlevel : int, default: 2
        The level at which to begin iteration (default: 2). This does not
        change the total number of function evaluations or the abscissae at
        which the function is evaluated; it changes only the *number of times*
        `f` is called. If ``minlevel=k``, then the integrand is evaluated at
        all abscissae from levels ``0`` through ``k`` in a single call.
        Note that if `minlevel` exceeds `maxlevel`, the provided `minlevel` is
        ignored, and `minlevel` is set equal to `maxlevel`.
    atol, rtol : float, optional
        Absolute termination tolerance (default: 0) and relative termination
        tolerance (default: ``eps**0.75``, where ``eps`` is the precision of
        the result dtype), respectively. The error estimate is as
        described in [1]_ Section 5. While not theoretically rigorous or
        conservative, it is said to work well in practice. Must be non-negative
        and finite if `log` is False, and must be expressed as the log of a
        non-negative and finite number if `log` is True.
    preserve_shape : bool, default: False
        In the following, "arguments of `f`" refers to the array ``x`` and
        any arrays within ``fargs``. Let ``shape`` be the broadcasted shape
        of `a`, `b`, and all elements of `args` (which is conceptually
        distinct from ``fargs`` passed into `f`).

        - When ``preserve_shape=False`` (default), `f` must accept arguments
          of *any* broadcastable shapes.

        - When ``preserve_shape=True``, `f` must accept arguments of shape
          ``shape`` *or* ``shape + (n,)``, where ``(n,)`` is the number of
          abscissae at which the function is being evaluated.

        In either case, for each scalar element ``xi`` within `x`, the array
        returned by `f` must include the scalar ``f(xi)`` at the same index.
        Consequently, the shape of the output is always the shape of the input
        ``x``.

        See Examples.

    callback : callable, optional
        An optional user-supplied function to be called before the first
        iteration and after each iteration.
        Called as ``callback(res)``, where ``res`` is a ``_RichResult``
        similar to that returned by `_differentiate` (but containing the
        current iterate's values of all variables). If `callback` raises a
        ``StopIteration``, the algorithm will terminate immediately and
        `_tanhsinh` will return a result object.

    Returns
    -------
    res : _RichResult
        An instance of `scipy._lib._util._RichResult` with the following
        attributes. (The descriptions are written as though the values will be
        scalars; however, if `func` returns an array, the outputs will be
        arrays of the same shape.)
        success : bool
            ``True`` when the algorithm terminated successfully (status ``0``).
        status : int
            An integer representing the exit status of the algorithm.
            ``0`` : The algorithm converged to the specified tolerances.
            ``-1`` : (unused)
            ``-2`` : The maximum number of iterations was reached.
            ``-3`` : A non-finite value was encountered.
            ``-4`` : Iteration was terminated by `callback`.
            ``1`` : The algorithm is proceeding normally (in `callback` only).
        integral : float
            An estimate of the integral
        error : float
            An estimate of the error. Only available if level two or higher
            has been completed; otherwise NaN.
        maxlevel : int
            The maximum refinement level used.
        nfev : int
            The number of points at which `func` was evaluated.

    See Also
    --------
    quad, quadrature

    Notes
    -----
    Implements the algorithm as described in [1]_ with minor adaptations for
    finite-precision arithmetic, including some described by [2]_ and [3]_. The
    tanh-sinh scheme was originally introduced in [4]_.

    Due to floating-point error in the abscissae, the function may be evaluated
    at the endpoints of the interval during iterations. The values returned by
    the function at the endpoints will be ignored.

    References
    ----------
    [1] Bailey, David H., Karthik Jeyabalan, and Xiaoye S. Li. "A comparison of
        three high-precision quadrature schemes." Experimental Mathematics 14.3
        (2005): 317-329.
    [2] Vanherck, Joren, Bart Sorée, and Wim Magnus. "Tanh-sinh quadrature for
        single and multiple integration using floating-point arithmetic."
        arXiv preprint arXiv:2007.15057 (2020).
    [3] van Engelen, Robert A.  "Improving the Double Exponential Quadrature
        Tanh-Sinh, Sinh-Sinh and Exp-Sinh Formulas."
        https://www.genivia.com/files/qthsh.pdf
    [4] Takahasi, Hidetosi, and Masatake Mori. "Double exponential formulas for
        numerical integration." Publications of the Research Institute for
        Mathematical Sciences 9.3 (1974): 721-741.

    Example
    -------
    Evaluate the Gaussian integral:

    >>> import numpy as np
    >>> from scipy.integrate._tanhsinh import _tanhsinh
    >>> def f(x):
    ...     return np.exp(-x**2)
    >>> res = _tanhsinh(f, -np.inf, np.inf)
    >>> res.integral  # true value is np.sqrt(np.pi), 1.7724538509055159
     1.7724538509055159
    >>> res.error  # actual error is 0
    4.0007963937534104e-16

    The value of the Gaussian function (bell curve) is nearly zero for
    arguments sufficiently far from zero, so the value of the integral
    over a finite interval is nearly the same.

    >>> _tanhsinh(f, -20, 20).integral
    1.772453850905518

    However, with unfavorable integration limits, the integration scheme
    may not be able to find the important region.

    >>> _tanhsinh(f, -np.inf, 1000).integral
    4.500490856620352

    In such cases, or when there are singularities within the interval,
    break the integral into parts with endpoints at the important points.

    >>> _tanhsinh(f, -np.inf, 0).integral + _tanhsinh(f, 0, 1000).integral
    1.772453850905404

    For integration involving very large or very small magnitudes, use
    log-integration. (For illustrative purposes, the following example shows a
    case in which both regular and log-integration work, but for more extreme
    limits of integration, log-integration would avoid the underflow
    experienced when evaluating the integral normally.)

    >>> res = _tanhsinh(f, 20, 30, rtol=1e-10)
    >>> res.integral, res.error
    4.7819613911309014e-176, 4.670364401645202e-187
    >>> def log_f(x):
    ...     return -x**2
    >>> np.exp(res.integral), np.exp(res.error)
    4.7819613911306924e-176, 4.670364401645093e-187

    The limits of integration and elements of `args` may be broadcastable
    arrays, and integration is performed elementwise.

    >>> from scipy import stats
    >>> dist = stats.gausshyper(13.8, 3.12, 2.51, 5.18)
    >>> a, b = dist.support()
    >>> x = np.linspace(a, b, 100)
    >>> res = _tanhsinh(dist.pdf, a, x)
    >>> ref = dist.cdf(x)
    >>> np.allclose(res.integral, ref)

    By default, `preserve_shape` is False, and therefore the callable
    `f` may be called with arrays of any broadcastable shapes.
    For example:

    >>> shapes = []
    >>> def f(x, c):
    ...    shape = np.broadcast_shapes(x.shape, c.shape)
    ...    shapes.append(shape)
    ...    return np.sin(c*x)
    >>>
    >>> c = [1, 10, 30, 100]
    >>> res = _tanhsinh(f, 0, 1, args=(c,), minlevel=1)
    >>> shapes
    [(4,), (4, 66), (3, 64), (2, 128), (1, 256)]

    To understand where these shapes are coming from - and to better
    understand how `_tanhsinh` computes accurate results - note that
    higher values of ``c`` correspond with higher frequency sinusoids.
    The higher frequency sinusoids make the integrand more complicated,
    so more function evaluations are required to achieve the target
    accuracy:

    >>> res.nfev
    array([ 67, 131, 259, 515])

    The initial ``shape``, ``(4,)``, corresponds with evaluating the
    integrand at a single abscissa and all four frequencies; this is used
    for input validation and to determine the size and dtype of the arrays
    that store results. The next shape corresponds with evaluating the
    integrand at an initial grid of abscissae and all four frequencies.
    Successive calls to the function double the total number of abscissae at
    which the function has been evaluated. However, in later function
    evaluations, the integrand is evaluated at fewer frequencies because
    the corresponding integral has already converged to the required
    tolerance. This saves function evaluations to improve performance, but
    it requires the function to accept arguments of any shape.

    "Vector-valued" integrands, such as those written for use with
    `scipy.integrate.quad_vec`, are unlikely to satisfy this requirement.
    For example, consider

    >>> def f(x):
    ...    return [x, np.sin(10*x), np.cos(30*x), x*np.sin(100*x)**2]

    This integrand is not compatible with `_tanhsinh` as written; for instance,
    the shape of the output will not be the same as the shape of ``x``. Such a
    function *could* be converted to a compatible form with the introduction of
    additional parameters, but this would be inconvenient. In such cases,
    a simpler solution would be to use `preserve_shape`.

    >>> shapes = []
    >>> def f(x):
    ...     shapes.append(x.shape)
    ...     x0, x1, x2, x3 = x
    ...     return [x0, np.sin(10*x1), np.cos(30*x2), x3*np.sin(100*x3)]
    >>>
    >>> a = np.zeros(4)
    >>> res = _tanhsinh(f, a, 1, preserve_shape=True)
    >>> shapes
    [(4,), (4, 66), (4, 64), (4, 128), (4, 256)]

    Here, the broadcasted shape of `a` and `b` is ``(4,)``. With
    ``preserve_shape=True``, the function may be called with argument
    ``x`` of shape ``(4,)`` or ``(4, n)``, and this is what we observe.

    """
    (f, a, b, log, maxfun, maxlevel, minlevel,
     atol, rtol, args, preserve_shape, callback) = _tanhsinh_iv(
        f, a, b, log, maxfun, maxlevel, minlevel, atol,
        rtol, args, preserve_shape, callback)

    # Initialization
    # `eim._initialize` does several important jobs, including
    # ensuring that limits, each of the `args`, and the output of `f`
    # broadcast correctly and are of consistent types. To save a function
    # evaluation, I pass the midpoint of the integration interval. This comes
    # at a cost of some gymnastics to ensure that the midpoint has the right
    # shape and dtype. Did you know that 0d and >0d arrays follow different
    # type promotion rules?
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        c = ((a.ravel() + b.ravel())/2).reshape(a.shape)
        inf_a, inf_b = np.isinf(a), np.isinf(b)
        c[inf_a] = b[inf_a] - 1  # takes care of infinite a
        c[inf_b] = a[inf_b] + 1  # takes care of infinite b
        c[inf_a & inf_b] = 0  # takes care of infinite a and b
        temp = eim._initialize(f, (c,), args, complex_ok=True,
                               preserve_shape=preserve_shape)
    f, xs, fs, args, shape, dtype, xp = temp
    a = np.broadcast_to(a, shape).astype(dtype).ravel()
    b = np.broadcast_to(b, shape).astype(dtype).ravel()

    # Transform improper integrals
    a, b, a0, negative, abinf, ainf, binf = _transform_integrals(a, b)

    # Define variables we'll need
    nit, nfev = 0, 1  # one function evaluation performed above
    zero = -np.inf if log else 0
    pi = dtype.type(np.pi)
    maxiter = maxlevel - minlevel + 1
    eps = np.finfo(dtype).eps
    if rtol is None:
        rtol = 0.75*np.log(eps) if log else eps**0.75

    Sn = np.full(shape, zero, dtype=dtype).ravel()  # latest integral estimate
    Sn[np.isnan(a) | np.isnan(b) | np.isnan(fs[0])] = np.nan
    Sk = np.empty_like(Sn).reshape(-1, 1)[:, 0:0]  # all integral estimates
    aerr = np.full(shape, np.nan, dtype=dtype).ravel()  # absolute error
    status = np.full(shape, eim._EINPROGRESS, dtype=int).ravel()
    h0 = np.real(_get_base_step(dtype=dtype))  # base step

    # For term `d4` of error estimate ([1] Section 5), we need to keep the
    # most extreme abscissae and corresponding `fj`s, `wj`s in Euler-Maclaurin
    # sum. Here, we initialize these variables.
    xr0 = np.full(shape, -np.inf, dtype=dtype).ravel()
    fr0 = np.full(shape, np.nan, dtype=dtype).ravel()
    wr0 = np.zeros(shape, dtype=dtype).ravel()
    xl0 = np.full(shape, np.inf, dtype=dtype).ravel()
    fl0 = np.full(shape, np.nan, dtype=dtype).ravel()
    wl0 = np.zeros(shape, dtype=dtype).ravel()
    d4 = np.zeros(shape, dtype=dtype).ravel()

    work = _RichResult(
        Sn=Sn, Sk=Sk, aerr=aerr, h=h0, log=log, dtype=dtype, pi=pi, eps=eps,
        a=a.reshape(-1, 1), b=b.reshape(-1, 1),  # integration limits
        n=minlevel, nit=nit, nfev=nfev, status=status,  # iter/eval counts
        xr0=xr0, fr0=fr0, wr0=wr0, xl0=xl0, fl0=fl0, wl0=wl0, d4=d4,  # err est
        ainf=ainf, binf=binf, abinf=abinf, a0=a0.reshape(-1, 1))  # transforms
    # Constant scalars don't need to be put in `work` unless they need to be
    # passed outside `tanhsinh`. Examples: atol, rtol, h0, minlevel.

    # Correspondence between terms in the `work` object and the result
    res_work_pairs = [('status', 'status'), ('integral', 'Sn'),
                      ('error', 'aerr'), ('nit', 'nit'), ('nfev', 'nfev')]

    def pre_func_eval(work):
        # Determine abscissae at which to evaluate `f`
        work.h = h0 / 2**work.n
        xjc, wj = _get_pairs(work.n, h0, dtype=work.dtype,
                             inclusive=(work.n == minlevel))
        work.xj, work.wj = _transform_to_limits(xjc, wj, work.a, work.b)

        # Perform abscissae substitutions for infinite limits of integration
        xj = work.xj.copy()
        xj[work.abinf] = xj[work.abinf] / (1 - xj[work.abinf]**2)
        xj[work.binf] = 1/xj[work.binf] - 1 + work.a0[work.binf]
        xj[work.ainf] *= -1
        return xj

    def post_func_eval(x, fj, work):
        # Weight integrand as required by substitutions for infinite limits
        if work.log:
            fj[work.abinf] += (np.log(1 + work.xj[work.abinf] ** 2)
                               - 2*np.log(1 - work.xj[work.abinf] ** 2))
            fj[work.binf] -= 2 * np.log(work.xj[work.binf])
        else:
            fj[work.abinf] *= ((1 + work.xj[work.abinf]**2) /
                               (1 - work.xj[work.abinf]**2)**2)
            fj[work.binf] *= work.xj[work.binf]**-2.

        # Estimate integral with Euler-Maclaurin Sum
        fjwj, Sn = _euler_maclaurin_sum(fj, work)
        if work.Sk.shape[-1]:
            Snm1 = work.Sk[:, -1]
            Sn = (special.logsumexp([Snm1 - np.log(2), Sn], axis=0) if log
                  else Snm1 / 2 + Sn)

        work.fjwj = fjwj
        work.Sn = Sn

    def check_termination(work):
        """Terminate due to convergence or encountering non-finite values"""
        stop = np.zeros(work.Sn.shape, dtype=bool)

        # Terminate before first iteration if integration limits are equal
        if work.nit == 0:
            i = (work.a == work.b).ravel()  # ravel singleton dimension
            zero = -np.inf if log else 0
            work.Sn[i] = zero
            work.aerr[i] = zero
            work.status[i] = eim._ECONVERGED
            stop[i] = True
        else:
            # Terminate if convergence criterion is met
            work.rerr, work.aerr = _estimate_error(work)
            i = ((work.rerr < rtol) | (work.rerr + np.real(work.Sn) < atol) if log
                 else (work.rerr < rtol) | (work.rerr * abs(work.Sn) < atol))
            work.status[i] = eim._ECONVERGED
            stop[i] = True

        # Terminate if integral estimate becomes invalid
        if log:
            i = (np.isposinf(np.real(work.Sn)) | np.isnan(work.Sn)) & ~stop
        else:
            i = ~np.isfinite(work.Sn) & ~stop
        work.status[i] = eim._EVALUEERR
        stop[i] = True

        return stop

    def post_termination_check(work):
        work.n += 1
        work.Sk = np.concatenate((work.Sk, work.Sn[:, np.newaxis]), axis=-1)
        return

    def customize_result(res, shape):
        # If the integration limits were such that b < a, we reversed them
        # to perform the calculation, and the final result needs to be negated.
        if log and np.any(negative):
            pi = res['integral'].dtype.type(np.pi)
            j = np.complex64(1j)  # minimum complex type
            res['integral'] = res['integral'] + negative*pi*j
        else:
            res['integral'][negative] *= -1

        # For this algorithm, it seems more appropriate to report the maximum
        # level rather than the number of iterations in which it was performed.
        res['maxlevel'] = minlevel + res['nit'] - 1
        res['maxlevel'][res['nit'] == 0] = -1
        del res['nit']
        return shape

    # Suppress all warnings initially, since there are many places in the code
    # for which this is expected behavior.
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        res = eim._loop(work, callback, shape, maxiter, f, args, dtype, pre_func_eval,
                        post_func_eval, check_termination, post_termination_check,
                        customize_result, res_work_pairs, xp, preserve_shape)
    return res


def _get_base_step(dtype=np.float64):
    # Compute the base step length for the provided dtype. Theoretically, the
    # Euler-Maclaurin sum is infinite, but it gets cut off when either the
    # weights underflow or the abscissae cannot be distinguished from the
    # limits of integration. The latter happens to occur first for float32 and
    # float64, and it occurs when `xjc` (the abscissa complement)
    # in `_compute_pair` underflows. We can solve for the argument `tmax` at
    # which it will underflow using [2] Eq. 13.
    fmin = 4*np.finfo(dtype).tiny  # stay a little away from the limit
    tmax = np.arcsinh(np.log(2/fmin - 1) / np.pi)

    # Based on this, we can choose a base step size `h` for level 0.
    # The number of function evaluations will be `2 + m*2^(k+1)`, where `k` is
    # the level and `m` is an integer we get to choose. I choose
    # m = _N_BASE_STEPS = `8` somewhat arbitrarily, but a rationale is that a
    # power of 2 makes floating point arithmetic more predictable. It also
    # results in a base step size close to `1`, which is what [1] uses (and I
    # used here until I found [2] and these ideas settled).
    h0 = tmax / _N_BASE_STEPS
    return h0.astype(dtype)


_N_BASE_STEPS = 8


def _compute_pair(k, h0):
    # Compute the abscissa-weight pairs for each level k. See [1] page 9.

    # For now, we compute and store in 64-bit precision. If higher-precision
    # data types become better supported, it would be good to compute these
    # using the highest precision available. Or, once there is an Array API-
    # compatible arbitrary precision array, we can compute at the required
    # precision.

    # "....each level k of abscissa-weight pairs uses h = 2 **-k"
    # We adapt to floating point arithmetic using ideas of [2].
    h = h0 / 2**k
    max = _N_BASE_STEPS * 2**k

    # For iterations after the first, "....the integrand function needs to be
    # evaluated only at the odd-indexed abscissas at each level."
    j = np.arange(max+1) if k == 0 else np.arange(1, max+1, 2)
    jh = j * h

    # "In this case... the weights wj = u1/cosh(u2)^2, where..."
    pi_2 = np.pi / 2
    u1 = pi_2*np.cosh(jh)
    u2 = pi_2*np.sinh(jh)
    # Denominators get big here. Overflow then underflow doesn't need warning.
    # with np.errstate(under='ignore', over='ignore'):
    wj = u1 / np.cosh(u2)**2
    # "We actually store 1-xj = 1/(...)."
    xjc = 1 / (np.exp(u2) * np.cosh(u2))  # complement of xj = np.tanh(u2)

    # When level k == 0, the zeroth xj corresponds with xj = 0. To simplify
    # code, the function will be evaluated there twice; each gets half weight.
    wj[0] = wj[0] / 2 if k == 0 else wj[0]

    return xjc, wj  # store at full precision


def _pair_cache(k, h0):
    # Cache the abscissa-weight pairs up to a specified level.
    # Abscissae and weights of consecutive levels are concatenated.
    # `index` records the indices that correspond with each level:
    # `xjc[index[k]:index[k+1]` extracts the level `k` abscissae.
    if h0 != _pair_cache.h0:
        _pair_cache.xjc = np.empty(0)
        _pair_cache.wj = np.empty(0)
        _pair_cache.indices = [0]

    xjcs = [_pair_cache.xjc]
    wjs = [_pair_cache.wj]

    for i in range(len(_pair_cache.indices)-1, k + 1):
        xjc, wj = _compute_pair(i, h0)
        xjcs.append(xjc)
        wjs.append(wj)
        _pair_cache.indices.append(_pair_cache.indices[-1] + len(xjc))

    _pair_cache.xjc = np.concatenate(xjcs)
    _pair_cache.wj = np.concatenate(wjs)
    _pair_cache.h0 = h0

_pair_cache.xjc = np.empty(0)
_pair_cache.wj = np.empty(0)
_pair_cache.indices = [0]
_pair_cache.h0 = None


def _get_pairs(k, h0, inclusive=False, dtype=np.float64):
    # Retrieve the specified abscissa-weight pairs from the cache
    # If `inclusive`, return all up to and including the specified level
    if len(_pair_cache.indices) <= k+2 or h0 != _pair_cache.h0:
        _pair_cache(k, h0)

    xjc = _pair_cache.xjc
    wj = _pair_cache.wj
    indices = _pair_cache.indices

    start = 0 if inclusive else indices[k]
    end = indices[k+1]

    return xjc[start:end].astype(dtype), wj[start:end].astype(dtype)


def _transform_to_limits(xjc, wj, a, b):
    # Transform integral according to user-specified limits. This is just
    # math that follows from the fact that the standard limits are (-1, 1).
    # Note: If we had stored xj instead of xjc, we would have
    # xj = alpha * xj + beta, where beta = (a + b)/2
    alpha = (b - a) / 2
    xj = np.concatenate((-alpha * xjc + b, alpha * xjc + a), axis=-1)
    wj = wj*alpha  # arguments get broadcasted, so we can't use *=
    wj = np.concatenate((wj, wj), axis=-1)

    # Points at the boundaries can be generated due to finite precision
    # arithmetic, but these function values aren't supposed to be included in
    # the Euler-Maclaurin sum. Ideally we wouldn't evaluate the function at
    # these points; however, we can't easily filter out points since this
    # function is vectorized. Instead, zero the weights.
    invalid = (xj <= a) | (xj >= b)
    wj[invalid] = 0
    return xj, wj


def _euler_maclaurin_sum(fj, work):
    # Perform the Euler-Maclaurin Sum, [1] Section 4

    # The error estimate needs to know the magnitude of the last term
    # omitted from the Euler-Maclaurin sum. This is a bit involved because
    # it may have been computed at a previous level. I sure hope it's worth
    # all the trouble.
    xr0, fr0, wr0 = work.xr0, work.fr0, work.wr0
    xl0, fl0, wl0 = work.xl0, work.fl0, work.wl0

    # It is much more convenient to work with the transposes of our work
    # variables here.
    xj, fj, wj = work.xj.T, fj.T, work.wj.T
    n_x, n_active = xj.shape  # number of abscissae, number of active elements

    # We'll work with the left and right sides separately
    xr, xl = xj.reshape(2, n_x // 2, n_active).copy()  # this gets modified
    fr, fl = fj.reshape(2, n_x // 2, n_active)
    wr, wl = wj.reshape(2, n_x // 2, n_active)

    invalid_r = ~np.isfinite(fr) | (wr == 0)
    invalid_l = ~np.isfinite(fl) | (wl == 0)

    # integer index of the maximum abscissa at this level
    xr[invalid_r] = -np.inf
    ir = np.argmax(xr, axis=0, keepdims=True)
    # abscissa, function value, and weight at this index
    xr_max = np.take_along_axis(xr, ir, axis=0)[0]
    fr_max = np.take_along_axis(fr, ir, axis=0)[0]
    wr_max = np.take_along_axis(wr, ir, axis=0)[0]
    # boolean indices at which maximum abscissa at this level exceeds
    # the incumbent maximum abscissa (from all previous levels)
    j = xr_max > xr0
    # Update record of the incumbent abscissa, function value, and weight
    xr0[j] = xr_max[j]
    fr0[j] = fr_max[j]
    wr0[j] = wr_max[j]

    # integer index of the minimum abscissa at this level
    xl[invalid_l] = np.inf
    il = np.argmin(xl, axis=0, keepdims=True)
    # abscissa, function value, and weight at this index
    xl_min = np.take_along_axis(xl, il, axis=0)[0]
    fl_min = np.take_along_axis(fl, il, axis=0)[0]
    wl_min = np.take_along_axis(wl, il, axis=0)[0]
    # boolean indices at which minimum abscissa at this level is less than
    # the incumbent minimum abscissa (from all previous levels)
    j = xl_min < xl0
    # Update record of the incumbent abscissa, function value, and weight
    xl0[j] = xl_min[j]
    fl0[j] = fl_min[j]
    wl0[j] = wl_min[j]
    fj = fj.T

    # Compute the error estimate `d4` - the magnitude of the leftmost or
    # rightmost term, whichever is greater.
    flwl0 = fl0 + np.log(wl0) if work.log else fl0 * wl0  # leftmost term
    frwr0 = fr0 + np.log(wr0) if work.log else fr0 * wr0  # rightmost term
    magnitude = np.real if work.log else np.abs
    work.d4 = np.maximum(magnitude(flwl0), magnitude(frwr0))

    # There are two approaches to dealing with function values that are
    # numerically infinite due to approaching a singularity - zero them, or
    # replace them with the function value at the nearest non-infinite point.
    # [3] pg. 22 suggests the latter, so let's do that given that we have the
    # information.
    fr0b = np.broadcast_to(fr0[np.newaxis, :], fr.shape)
    fl0b = np.broadcast_to(fl0[np.newaxis, :], fl.shape)
    fr[invalid_r] = fr0b[invalid_r]
    fl[invalid_l] = fl0b[invalid_l]

    # When wj is zero, log emits a warning
    # with np.errstate(divide='ignore'):
    fjwj = fj + np.log(work.wj) if work.log else fj * work.wj

    # update integral estimate
    Sn = (special.logsumexp(fjwj + np.log(work.h), axis=-1) if work.log
          else np.sum(fjwj, axis=-1) * work.h)

    work.xr0, work.fr0, work.wr0 = xr0, fr0, wr0
    work.xl0, work.fl0, work.wl0 = xl0, fl0, wl0

    return fjwj, Sn


def _estimate_error(work):
    # Estimate the error according to [1] Section 5

    if work.n == 0 or work.nit == 0:
        # The paper says to use "one" as the error before it can be calculated.
        # NaN seems to be more appropriate.
        nan = np.full_like(work.Sn, np.nan)
        return nan, nan

    indices = _pair_cache.indices

    n_active = len(work.Sn)  # number of active elements
    axis_kwargs = dict(axis=-1, keepdims=True)

    # With a jump start (starting at level higher than 0), we haven't
    # explicitly calculated the integral estimate at lower levels. But we have
    # all the function value-weight products, so we can compute the
    # lower-level estimates.
    if work.Sk.shape[-1] == 0:
        h = 2 * work.h  # step size at this level
        n_x = indices[work.n]  # number of abscissa up to this level
        # The right and left fjwj terms from all levels are concatenated along
        # the last axis. Get out only the terms up to this level.
        fjwj_rl = work.fjwj.reshape(n_active, 2, -1)
        fjwj = fjwj_rl[:, :, :n_x].reshape(n_active, 2*n_x)
        # Compute the Euler-Maclaurin sum at this level
        Snm1 = (special.logsumexp(fjwj, **axis_kwargs) + np.log(h) if work.log
                else np.sum(fjwj, **axis_kwargs) * h)
        work.Sk = np.concatenate((Snm1, work.Sk), axis=-1)

    if work.n == 1:
        nan = np.full_like(work.Sn, np.nan)
        return nan, nan

    # The paper says not to calculate the error for n<=2, but it's not clear
    # about whether it starts at level 0 or level 1. We start at level 0, so
    # why not compute the error beginning in level 2?
    if work.Sk.shape[-1] < 2:
        h = 4 * work.h  # step size at this level
        n_x = indices[work.n-1]  # number of abscissa up to this level
        # The right and left fjwj terms from all levels are concatenated along
        # the last axis. Get out only the terms up to this level.
        fjwj_rl = work.fjwj.reshape(len(work.Sn), 2, -1)
        fjwj = fjwj_rl[..., :n_x].reshape(n_active, 2*n_x)
        # Compute the Euler-Maclaurin sum at this level
        Snm2 = (special.logsumexp(fjwj, **axis_kwargs) + np.log(h) if work.log
                else np.sum(fjwj, **axis_kwargs) * h)
        work.Sk = np.concatenate((Snm2, work.Sk), axis=-1)

    Snm2 = work.Sk[..., -2]
    Snm1 = work.Sk[..., -1]

    e1 = work.eps

    if work.log:
        log_e1 = np.log(e1)
        # Currently, only real integrals are supported in log-scale. All
        # complex values have imaginary part in increments of pi*j, which just
        # carries sign information of the original integral, so use of
        # `np.real` here is equivalent to absolute value in real scale.
        d1 = np.real(special.logsumexp([work.Sn, Snm1 + work.pi*1j], axis=0))
        d2 = np.real(special.logsumexp([work.Sn, Snm2 + work.pi*1j], axis=0))
        d3 = log_e1 + np.max(np.real(work.fjwj), axis=-1)
        d4 = work.d4
        aerr = np.max([d1 ** 2 / d2, 2 * d1, d3, d4], axis=0)
        rerr = np.maximum(log_e1, aerr - np.real(work.Sn))
    else:
        # Note: explicit computation of log10 of each of these is unnecessary.
        d1 = np.abs(work.Sn - Snm1)
        d2 = np.abs(work.Sn - Snm2)
        d3 = e1 * np.max(np.abs(work.fjwj), axis=-1)
        d4 = work.d4
        # If `d1` is 0, no need to warn. This does the right thing.
        # with np.errstate(divide='ignore'):
        aerr = np.max([d1**(np.log(d1)/np.log(d2)), d1**2, d3, d4], axis=0)
        rerr = np.maximum(e1, aerr/np.abs(work.Sn))
    return rerr, aerr.reshape(work.Sn.shape)


def _transform_integrals(a, b):
    # Transform integrals to a form with finite a < b
    # For b < a, we reverse the limits and will multiply the final result by -1
    # For infinite limit on the right, we use the substitution x = 1/t - 1 + a
    # For infinite limit on the left, we substitute x = -x and treat as above
    # For infinite limits, we substitute x = t / (1-t**2)

    negative = b < a
    a[negative], b[negative] = b[negative], a[negative]

    abinf = np.isinf(a) & np.isinf(b)
    a[abinf], b[abinf] = -1, 1

    ainf = np.isinf(a)
    a[ainf], b[ainf] = -b[ainf], -a[ainf]

    binf = np.isinf(b)
    a0 = a.copy()
    a[binf], b[binf] = 0, 1

    return a, b, a0, negative, abinf, ainf, binf


def _tanhsinh_iv(f, a, b, log, maxfun, maxlevel, minlevel,
                 atol, rtol, args, preserve_shape, callback):
    # Input validation and standardization

    message = '`f` must be callable.'
    if not callable(f):
        raise ValueError(message)

    message = 'All elements of `a` and `b` must be real numbers.'
    a, b = np.broadcast_arrays(a, b)
    if np.any(np.iscomplex(a)) or np.any(np.iscomplex(b)):
        raise ValueError(message)

    message = '`log` must be True or False.'
    if log not in {True, False}:
        raise ValueError(message)
    log = bool(log)

    if atol is None:
        atol = -np.inf if log else 0

    rtol_temp = rtol if rtol is not None else 0.

    params = np.asarray([atol, rtol_temp, 0.])
    message = "`atol` and `rtol` must be real numbers."
    if not np.issubdtype(params.dtype, np.floating):
        raise ValueError(message)

    if log:
        message = '`atol` and `rtol` may not be positive infinity.'
        if np.any(np.isposinf(params)):
            raise ValueError(message)
    else:
        message = '`atol` and `rtol` must be non-negative and finite.'
        if np.any(params < 0) or np.any(np.isinf(params)):
            raise ValueError(message)
    atol = params[0]
    rtol = rtol if rtol is None else params[1]

    BIGINT = float(2**62)
    if maxfun is None and maxlevel is None:
        maxlevel = 10

    maxfun = BIGINT if maxfun is None else maxfun
    maxlevel = BIGINT if maxlevel is None else maxlevel

    message = '`maxfun`, `maxlevel`, and `minlevel` must be integers.'
    params = np.asarray([maxfun, maxlevel, minlevel])
    if not (np.issubdtype(params.dtype, np.number)
            and np.all(np.isreal(params))
            and np.all(params.astype(np.int64) == params)):
        raise ValueError(message)
    message = '`maxfun`, `maxlevel`, and `minlevel` must be non-negative.'
    if np.any(params < 0):
        raise ValueError(message)
    maxfun, maxlevel, minlevel = params.astype(np.int64)
    minlevel = min(minlevel, maxlevel)

    if not np.iterable(args):
        args = (args,)

    message = '`preserve_shape` must be True or False.'
    if preserve_shape not in {True, False}:
        raise ValueError(message)

    if callback is not None and not callable(callback):
        raise ValueError('`callback` must be callable.')

    return (f, a, b, log, maxfun, maxlevel, minlevel,
            atol, rtol, args, preserve_shape, callback)


def _logsumexp(x, axis=0):
    # logsumexp raises with empty array
    x = np.asarray(x)
    shape = list(x.shape)
    if shape[axis] == 0:
        shape.pop(axis)
        return np.full(shape, fill_value=-np.inf, dtype=x.dtype)
    else:
        return special.logsumexp(x, axis=axis)


def _nsum_iv(f, a, b, step, args, log, maxterms, atol, rtol):
    # Input validation and standardization

    message = '`f` must be callable.'
    if not callable(f):
        raise ValueError(message)

    message = 'All elements of `a`, `b`, and `step` must be real numbers.'
    a, b, step = np.broadcast_arrays(a, b, step)
    dtype = np.result_type(a.dtype, b.dtype, step.dtype)
    if not np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.complexfloating):
        raise ValueError(message)

    valid_a = np.isfinite(a)
    valid_b = b >= a  # NaNs will be False
    valid_step = np.isfinite(step) & (step > 0)
    valid_abstep = valid_a & valid_b & valid_step

    message = '`log` must be True or False.'
    if log not in {True, False}:
        raise ValueError(message)

    if atol is None:
        atol = -np.inf if log else 0

    rtol_temp = rtol if rtol is not None else 0.

    params = np.asarray([atol, rtol_temp, 0.])
    message = "`atol` and `rtol` must be real numbers."
    if not np.issubdtype(params.dtype, np.floating):
        raise ValueError(message)

    if log:
        message = '`atol`, `rtol` may not be positive infinity or NaN.'
        if np.any(np.isposinf(params) | np.isnan(params)):
            raise ValueError(message)
    else:
        message = '`atol`, and `rtol` must be non-negative and finite.'
        if np.any((params < 0) | (~np.isfinite(params))):
            raise ValueError(message)
    atol = params[0]
    rtol = rtol if rtol is None else params[1]

    maxterms_int = int(maxterms)
    if maxterms_int != maxterms or maxterms < 0:
        message = "`maxterms` must be a non-negative integer."
        raise ValueError(message)

    if not np.iterable(args):
        args = (args,)

    return f, a, b, step, valid_abstep, args, log, maxterms_int, atol, rtol


def _nsum(f, a, b, step=1, args=(), log=False, maxterms=int(2**20), atol=None,
          rtol=None):
    r"""Evaluate a convergent sum.

    For finite `b`, this evaluates::

        f(a + np.arange(n)*step).sum()

    where ``n = int((b - a) / step) + 1``. If `f` is smooth, positive, and
    monotone decreasing, `b` may be infinite, in which case the infinite sum
    is approximated using integration.

    Parameters
    ----------
    f : callable
        The function that evaluates terms to be summed. The signature must be::

            f(x: ndarray, *args) -> ndarray

         where each element of ``x`` is a finite real and ``args`` is a tuple,
         which may contain an arbitrary number of arrays that are broadcastable
         with `x`. `f` must represent a smooth, positive, and monotone decreasing
         function of `x`; `_nsum` performs no checks to verify that these conditions
         are met and may return erroneous results if they are violated.
    a, b : array_like
        Real lower and upper limits of summed terms. Must be broadcastable.
        Each element of `a` must be finite and less than the corresponding
        element in `b`, but elements of `b` may be infinite.
    step : array_like
        Finite, positive, real step between summed terms. Must be broadcastable
        with `a` and `b`.
    args : tuple, optional
        Additional positional arguments to be passed to `f`. Must be arrays
        broadcastable with `a`, `b`, and `step`. If the callable to be summed
        requires arguments that are not broadcastable with `a`, `b`, and `step`,
        wrap that callable with `f`. See Examples.
    log : bool, default: False
        Setting to True indicates that `f` returns the log of the terms
        and that `atol` and `rtol` are expressed as the logs of the absolute
        and relative errors. In this case, the result object will contain the
        log of the sum and error. This is useful for summands for which
        numerical underflow or overflow would lead to inaccuracies.
    maxterms : int, default: 2**32
        The maximum number of terms to evaluate when summing directly. 
        Additional function evaluations may be performed for input
        validation and integral evaluation. 
    atol, rtol : float, optional
        Absolute termination tolerance (default: 0) and relative termination
        tolerance (default: ``eps**0.5``, where ``eps`` is the precision of
        the result dtype), respectively. Must be non-negative
        and finite if `log` is False, and must be expressed as the log of a
        non-negative and finite number if `log` is True.

    Returns
    -------
    res : _RichResult
        An instance of `scipy._lib._util._RichResult` with the following
        attributes. (The descriptions are written as though the values will be
        scalars; however, if `func` returns an array, the outputs will be

        arrays of the same shape.)
        success : bool
            ``True`` when the algorithm terminated successfully (status ``0``).
        status : int
            An integer representing the exit status of the algorithm.
            ``0`` : The algorithm converged to the specified tolerances.
            ``-1`` : Element(s) of `a`, `b`, or `step` are invalid
            ``-2`` : Numerical integration reached its iteration limit; the sum may be divergent.
            ``-3`` : A non-finite value was encountered.
        sum : float
            An estimate of the sum.
        error : float
            An estimate of the absolute error, assuming all terms are non-negative.
        nfev : int
            The number of points at which `func` was evaluated.

    See Also
    --------
    tanhsinh

    Notes
    -----
    The method implemented for infinite summation is related to the integral
    test for convergence of an infinite series: assuming `step` size 1 for
    simplicity of exposition, the sum of a monotone decreasing function is bounded by

    .. math::

        \int_u^\infty f(x) dx \leq \sum_{k=u}^\infty f(k) \leq \int_u^\infty f(x) dx + f(u)

    Let :math:`a` represent  `a`, :math:`n` represent `maxterms`, :math:`\epsilon_a`
    represent `atol`, and :math:`\epsilon_r` represent `rtol`.
    The implementation first evaluates the integral :math:`S_l=\int_a^\infty f(x) dx`
    as a lower bound of the infinite sum. Then, it seeks a value :math:`c > a` such
    that :math:`f(c) < \epsilon_a + S_l \epsilon_r`, if it exists; otherwise,
    let :math:`c = a + n`. Then the infinite sum is approximated as
    
    .. math::

        \sum_{k=a}^{c-1} f(k) + \int_c^\infty f(x) dx + f(c)/2

    and the reported error is :math:`f(c)/2` plus the error estimate of
    numerical integration. The approach described above is generalized for non-unit
    `step` and finite `b` that is too large for direct evaluation of the sum,
    i.e. ``b - a + 1 > maxterms``.

    References
    ----------
    [1] Wikipedia. "Integral test for convergence."
    https://en.wikipedia.org/wiki/Integral_test_for_convergence

    Examples
    --------
    Compute the infinite sum of the reciprocals of squared integers.
    
    >>> import numpy as np
    >>> from scipy.integrate._tanhsinh import _nsum
    >>> res = _nsum(lambda k: 1/k**2, 1, np.inf, maxterms=1e3)
    >>> ref = np.pi**2/6  # true value
    >>> res.error  # estimated error
    4.990014980029223e-07
    >>> (res.sum - ref)/ref  # true error
    -1.0101760641302586e-10
    >>> res.nfev  # number of points at which callable was evaluated
    1142
    
    Compute the infinite sums of the reciprocals of integers raised to powers ``p``.
    
    >>> from scipy import special
    >>> p = np.arange(2, 10)
    >>> res = _nsum(lambda k, p: 1/k**p, 1, np.inf, maxterms=1e3, args=(p,))
    >>> ref = special.zeta(p, 1)
    >>> np.allclose(res.sum, ref)
    True
    
    """ # noqa: E501
    # Potential future work:
    # - more careful testing of when `b` is slightly less than `a` plus an
    #   integer multiple of step (needed before this is public)
    # - improve error estimate of `_direct` sum
    # - add other methods for convergence acceleration (Richardson, epsilon)
    # - support infinite lower limit?
    # - support negative monotone increasing functions?
    # - b < a / negative step?
    # - complex-valued function?
    # - check for violations of monotonicity?

    # Function-specific input validation / standardization
    tmp = _nsum_iv(f, a, b, step, args, log, maxterms, atol, rtol)
    f, a, b, step, valid_abstep, args, log, maxterms, atol, rtol = tmp

    # Additional elementwise algorithm input validation / standardization
    tmp = eim._initialize(f, (a,), args, complex_ok=False)
    f, xs, fs, args, shape, dtype, xp = tmp

    # Finish preparing `a`, `b`, and `step` arrays
    a = xs[0]
    b = np.broadcast_to(b, shape).ravel().astype(dtype)
    step = np.broadcast_to(step, shape).ravel().astype(dtype)
    valid_abstep = np.broadcast_to(valid_abstep, shape).ravel()
    nterms = np.floor((b - a) / step)
    b = a + nterms*step

    # Define constants
    eps = np.finfo(dtype).eps
    zero = np.asarray(-np.inf if log else 0, dtype=dtype)[()]
    if rtol is None:
        rtol = 0.5*np.log(eps) if log else eps**0.5
    constants = (dtype, log, eps, zero, rtol, atol, maxterms)

    # Prepare result arrays
    S = np.empty_like(a)
    E = np.empty_like(a)
    status = np.zeros(len(a), dtype=int)
    nfev = np.ones(len(a), dtype=int)  # one function evaluation above

    # Branch for direct sum evaluation / integral approximation / invalid input
    i1 = (nterms + 1 <= maxterms) & valid_abstep
    i2 = (nterms + 1 > maxterms) & valid_abstep
    i3 = ~valid_abstep

    if np.any(i1):
        args_direct = [arg[i1] for arg in args]
        tmp = _direct(f, a[i1], b[i1], step[i1], args_direct, constants)
        S[i1], E[i1] = tmp[:-1]
        nfev[i1] += tmp[-1]
        status[i1] = -3 * (~np.isfinite(S[i1]))

    if np.any(i2):
        args_indirect = [arg[i2] for arg in args]
        tmp = _integral_bound(f, a[i2], b[i2], step[i2], args_indirect, constants)
        S[i2], E[i2], status[i2] = tmp[:-1]
        nfev[i2] += tmp[-1]

    if np.any(i3):
        S[i3], E[i3] = np.nan, np.nan
        status[i3] = -1

    # Return results
    S, E = S.reshape(shape)[()], E.reshape(shape)[()]
    status, nfev = status.reshape(shape)[()], nfev.reshape(shape)[()]
    return _RichResult(sum=S, error=E, status=status, success=status == 0,
                       nfev=nfev)


def _direct(f, a, b, step, args, constants, inclusive=True):
    # Directly evaluate the sum.

    # When used in the context of distributions, `args` would contain the
    # distribution parameters. We have broadcasted for simplicity, but we could
    # reduce function evaluations when distribution parameters are the same but
    # sum limits differ. Roughly:
    # - compute the function at all points between min(a) and max(b),
    # - compute the cumulative sum,
    # - take the difference between elements of the cumulative sum
    #   corresponding with b and a.
    # This is left to future enhancement

    dtype, log, eps, zero, _, _, _ = constants

    # To allow computation in a single vectorized call, find the maximum number
    # of points (over all slices) at which the function needs to be evaluated.
    # Note: if `inclusive` is `True`, then we want `1` more term in the sum.
    # I didn't think it was great style to use `True` as `1` in Python, so I
    # explicitly converted it to an `int` before using it.
    inclusive_adjustment = int(inclusive)
    steps = np.round((b - a) / step) + inclusive_adjustment
    # Equivalently, steps = np.round((b - a) / step) + inclusive
    max_steps = int(np.max(steps))

    # In each slice, the function will be evaluated at the same number of points,
    # but excessive points (those beyond the right sum limit `b`) are replaced
    # with NaN to (potentially) reduce the time of these unnecessary calculations.
    # Use a new last axis for these calculations for consistency with other
    # elementwise algorithms.
    a2, b2, step2 = a[:, np.newaxis], b[:, np.newaxis], step[:, np.newaxis]
    args2 = [arg[:, np.newaxis] for arg in args]
    ks = a2 + np.arange(max_steps, dtype=dtype) * step2
    i_nan = ks >= (b2 + inclusive_adjustment*step2/2)
    ks[i_nan] = np.nan
    fs = f(ks, *args2)

    # The function evaluated at NaN is NaN, and NaNs are zeroed in the sum.
    # In some cases it may be faster to loop over slices than to vectorize
    # like this. This is an optimization that can be added later.
    fs[i_nan] = zero
    nfev = max_steps - i_nan.sum(axis=-1)
    S = _logsumexp(fs, axis=-1) if log else np.sum(fs, axis=-1)
    # Rough, non-conservative error estimate. See gh-19667 for improvement ideas.
    E = np.real(S) + np.log(eps) if log else eps * abs(S)
    return S, E, nfev


def _integral_bound(f, a, b, step, args, constants):
    # Estimate the sum with integral approximation
    dtype, log, _, _, rtol, atol, maxterms = constants
    log2 = np.log(2, dtype=dtype)

    # Get a lower bound on the sum and compute effective absolute tolerance
    lb = _tanhsinh(f, a, b, args=args, atol=atol, rtol=rtol, log=log)
    tol = np.broadcast_to(atol, lb.integral.shape)
    tol = _logsumexp((tol, rtol + lb.integral)) if log else tol + rtol*lb.integral
    i_skip = lb.status < 0  # avoid unnecessary f_evals if integral is divergent
    tol[i_skip] = np.nan
    status = lb.status

    # As in `_direct`, we'll need a temporary new axis for points
    # at which to evaluate the function. Append axis at the end for
    # consistency with other elementwise algorithms.
    a2 = a[..., np.newaxis]
    step2 = step[..., np.newaxis]
    args2 = [arg[..., np.newaxis] for arg in args]

    # Find the location of a term that is less than the tolerance (if possible)
    log2maxterms = np.floor(np.log2(maxterms)) if maxterms else 0
    n_steps = np.concatenate([2**np.arange(0, log2maxterms), [maxterms]], dtype=dtype)
    nfev = len(n_steps)
    ks = a2 + n_steps * step2
    fks = f(ks, *args2)
    nt = np.minimum(np.sum(fks > tol[:, np.newaxis], axis=-1),  n_steps.shape[-1]-1)
    n_steps = n_steps[nt]

    # Directly evaluate the sum up to this term
    k = a + n_steps * step
    left, left_error, left_nfev = _direct(f, a, k, step, args,
                                          constants, inclusive=False)
    i_skip |= np.isposinf(left)  # if sum is not finite, no sense in continuing
    status[np.isposinf(left)] = -3
    k[i_skip] = np.nan

    # Use integration to estimate the remaining sum
    # Possible optimization for future work: if there were no terms less than
    # the tolerance, there is no need to compute the integral to better accuracy.
    # Something like:
    # atol = np.maximum(atol, np.minimum(fk/2 - fb/2))
    # rtol = np.maximum(rtol, np.minimum((fk/2 - fb/2)/left))
    # where `fk`/`fb` are currently calculated below.
    right = _tanhsinh(f, k, b, args=args, atol=atol, rtol=rtol, log=log)

    # Calculate the full estimate and error from the pieces
    fk = fks[np.arange(len(fks)), nt]
    fb = f(b, *args)
    nfev += 1
    if log:
        log_step = np.log(step)
        S_terms = (left, right.integral - log_step, fk - log2, fb - log2)
        S = _logsumexp(S_terms, axis=0)
        E_terms = (left_error, right.error - log_step, fk-log2, fb-log2+np.pi*1j)
        E = _logsumexp(E_terms, axis=0).real
    else:
        S = left + right.integral/step + fk/2 + fb/2
        E = left_error + right.error/step + fk/2 - fb/2
    status[~i_skip] = right.status[~i_skip]
    return S, E, status, left_nfev + right.nfev + nfev + lb.nfev
