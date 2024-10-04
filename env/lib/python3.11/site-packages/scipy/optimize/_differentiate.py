# mypy: disable-error-code="attr-defined"
import numpy as np
import scipy._lib._elementwise_iterative_method as eim
from scipy._lib._util import _RichResult

_EERRORINCREASE = -1  # used in _differentiate

def _differentiate_iv(func, x, args, atol, rtol, maxiter, order, initial_step,
                      step_factor, step_direction, preserve_shape, callback):
    # Input validation for `_differentiate`

    if not callable(func):
        raise ValueError('`func` must be callable.')

    # x has more complex IV that is taken care of during initialization
    x = np.asarray(x)
    dtype = x.dtype if np.issubdtype(x.dtype, np.inexact) else np.float64

    if not np.iterable(args):
        args = (args,)

    if atol is None:
        atol = np.finfo(dtype).tiny

    if rtol is None:
        rtol = np.sqrt(np.finfo(dtype).eps)

    message = 'Tolerances and step parameters must be non-negative scalars.'
    tols = np.asarray([atol, rtol, initial_step, step_factor])
    if (not np.issubdtype(tols.dtype, np.number)
            or np.any(tols < 0)
            or tols.shape != (4,)):
        raise ValueError(message)
    initial_step, step_factor = tols[2:].astype(dtype)

    maxiter_int = int(maxiter)
    if maxiter != maxiter_int or maxiter <= 0:
        raise ValueError('`maxiter` must be a positive integer.')

    order_int = int(order)
    if order_int != order or order <= 0:
        raise ValueError('`order` must be a positive integer.')

    step_direction = np.sign(step_direction).astype(dtype)
    x, step_direction = np.broadcast_arrays(x, step_direction)
    x, step_direction = x[()], step_direction[()]

    message = '`preserve_shape` must be True or False.'
    if preserve_shape not in {True, False}:
        raise ValueError(message)

    if callback is not None and not callable(callback):
        raise ValueError('`callback` must be callable.')

    return (func, x, args, atol, rtol, maxiter_int, order_int, initial_step,
            step_factor, step_direction, preserve_shape, callback)


def _differentiate(func, x, *, args=(), atol=None, rtol=None, maxiter=10,
                   order=8, initial_step=0.5, step_factor=2.0,
                   step_direction=0, preserve_shape=False, callback=None):
    """Evaluate the derivative of an elementwise scalar function numerically.

    Parameters
    ----------
    func : callable
        The function whose derivative is desired. The signature must be::

            func(x: ndarray, *fargs) -> ndarray

         where each element of ``x`` is a finite real number and ``fargs`` is a tuple,
         which may contain an arbitrary number of arrays that are broadcastable
         with `x`. ``func`` must be an elementwise function: each element
         ``func(x)[i]`` must equal ``func(x[i])`` for all indices ``i``.
    x : array_like
        Abscissae at which to evaluate the derivative.
    args : tuple, optional
        Additional positional arguments to be passed to `func`. Must be arrays
        broadcastable with `x`. If the callable to be differentiated requires
        arguments that are not broadcastable with `x`, wrap that callable with
        `func`. See Examples.
    atol, rtol : float, optional
        Absolute and relative tolerances for the stopping condition: iteration
        will stop when ``res.error < atol + rtol * abs(res.df)``. The default
        `atol` is the smallest normal number of the appropriate dtype, and
        the default `rtol` is the square root of the precision of the
        appropriate dtype.
    order : int, default: 8
        The (positive integer) order of the finite difference formula to be
        used. Odd integers will be rounded up to the next even integer.
    initial_step : float, default: 0.5
        The (absolute) initial step size for the finite difference derivative
        approximation.
    step_factor : float, default: 2.0
        The factor by which the step size is *reduced* in each iteration; i.e.
        the step size in iteration 1 is ``initial_step/step_factor``. If
        ``step_factor < 1``, subsequent steps will be greater than the initial
        step; this may be useful if steps smaller than some threshold are
        undesirable (e.g. due to subtractive cancellation error).
    maxiter : int, default: 10
        The maximum number of iterations of the algorithm to perform. See
        notes.
    step_direction : array_like
        An array representing the direction of the finite difference steps (for
        use when `x` lies near to the boundary of the domain of the function.)
        Must be broadcastable with `x` and all `args`.
        Where 0 (default), central differences are used; where negative (e.g.
        -1), steps are non-positive; and where positive (e.g. 1), all steps are
        non-negative.
    preserve_shape : bool, default: False
        In the following, "arguments of `func`" refers to the array ``x`` and
        any arrays within ``fargs``. Let ``shape`` be the broadcasted shape
        of `x` and all elements of `args` (which is conceptually
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
        `_differentiate` will return a result.

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
            ``-1`` : The error estimate increased, so iteration was terminated.
            ``-2`` : The maximum number of iterations was reached.
            ``-3`` : A non-finite value was encountered.
            ``-4`` : Iteration was terminated by `callback`.
            ``1`` : The algorithm is proceeding normally (in `callback` only).
        df : float
            The derivative of `func` at `x`, if the algorithm terminated
            successfully.
        error : float
            An estimate of the error: the magnitude of the difference between
            the current estimate of the derivative and the estimate in the
            previous iteration.
        nit : int
            The number of iterations performed.
        nfev : int
            The number of points at which `func` was evaluated.
        x : float
            The value at which the derivative of `func` was evaluated
            (after broadcasting with `args` and `step_direction`).

    Notes
    -----
    The implementation was inspired by jacobi [1]_, numdifftools [2]_, and
    DERIVEST [3]_, but the implementation follows the theory of Taylor series
    more straightforwardly (and arguably naively so).
    In the first iteration, the derivative is estimated using a finite
    difference formula of order `order` with maximum step size `initial_step`.
    Each subsequent iteration, the maximum step size is reduced by
    `step_factor`, and the derivative is estimated again until a termination
    condition is reached. The error estimate is the magnitude of the difference
    between the current derivative approximation and that of the previous
    iteration.

    The stencils of the finite difference formulae are designed such that
    abscissae are "nested": after `func` is evaluated at ``order + 1``
    points in the first iteration, `func` is evaluated at only two new points
    in each subsequent iteration; ``order - 1`` previously evaluated function
    values required by the finite difference formula are reused, and two
    function values (evaluations at the points furthest from `x`) are unused.

    Step sizes are absolute. When the step size is small relative to the
    magnitude of `x`, precision is lost; for example, if `x` is ``1e20``, the
    default initial step size of ``0.5`` cannot be resolved. Accordingly,
    consider using larger initial step sizes for large magnitudes of `x`.

    The default tolerances are challenging to satisfy at points where the
    true derivative is exactly zero. If the derivative may be exactly zero,
    consider specifying an absolute tolerance (e.g. ``atol=1e-16``) to
    improve convergence.

    References
    ----------
    [1]_ Hans Dembinski (@HDembinski). jacobi.
         https://github.com/HDembinski/jacobi
    [2]_ Per A. Brodtkorb and John D'Errico. numdifftools.
         https://numdifftools.readthedocs.io/en/latest/
    [3]_ John D'Errico. DERIVEST: Adaptive Robust Numerical Differentiation.
         https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation
    [4]_ Numerical Differentition. Wikipedia.
         https://en.wikipedia.org/wiki/Numerical_differentiation

    Examples
    --------
    Evaluate the derivative of ``np.exp`` at several points ``x``.

    >>> import numpy as np
    >>> from scipy.optimize._differentiate import _differentiate
    >>> f = np.exp
    >>> df = np.exp  # true derivative
    >>> x = np.linspace(1, 2, 5)
    >>> res = _differentiate(f, x)
    >>> res.df  # approximation of the derivative
    array([2.71828183, 3.49034296, 4.48168907, 5.75460268, 7.3890561 ])
    >>> res.error  # estimate of the error
    array(
        [7.12940817e-12, 9.16688947e-12, 1.17594823e-11, 1.50972568e-11, 1.93942640e-11]
    )
    >>> abs(res.df - df(x))  # true error
    array(
        [3.06421555e-14, 3.01980663e-14, 5.06261699e-14, 6.30606678e-14, 8.34887715e-14]
    )

    Show the convergence of the approximation as the step size is reduced.
    Each iteration, the step size is reduced by `step_factor`, so for
    sufficiently small initial step, each iteration reduces the error by a
    factor of ``1/step_factor**order`` until finite precision arithmetic
    inhibits further improvement.

    >>> iter = list(range(1, 12))  # maximum iterations
    >>> hfac = 2  # step size reduction per iteration
    >>> hdir = [-1, 0, 1]  # compare left-, central-, and right- steps
    >>> order = 4  # order of differentiation formula
    >>> x = 1
    >>> ref = df(x)
    >>> errors = []  # true error
    >>> for i in iter:
    ...     res = _differentiate(f, x, maxiter=i, step_factor=hfac,
    ...                          step_direction=hdir, order=order,
    ...                          atol=0, rtol=0)  # prevent early termination
    ...     errors.append(abs(res.df - ref))
    >>> errors = np.array(errors)
    >>> plt.semilogy(iter, errors[:, 0], label='left differences')
    >>> plt.semilogy(iter, errors[:, 1], label='central differences')
    >>> plt.semilogy(iter, errors[:, 2], label='right differences')
    >>> plt.xlabel('iteration')
    >>> plt.ylabel('error')
    >>> plt.legend()
    >>> plt.show()
    >>> (errors[1, 1] / errors[0, 1], 1 / hfac**order)
    (0.06215223140159822, 0.0625)

    The implementation is vectorized over `x`, `step_direction`, and `args`.
    The function is evaluated once before the first iteration to perform input
    validation and standardization, and once per iteration thereafter.

    >>> def f(x, p):
    ...     print('here')
    ...     f.nit += 1
    ...     return x**p
    >>> f.nit = 0
    >>> def df(x, p):
    ...     return p*x**(p-1)
    >>> x = np.arange(1, 5)
    >>> p = np.arange(1, 6).reshape((-1, 1))
    >>> hdir = np.arange(-1, 2).reshape((-1, 1, 1))
    >>> res = _differentiate(f, x, args=(p,), step_direction=hdir, maxiter=1)
    >>> np.allclose(res.df, df(x, p))
    True
    >>> res.df.shape
    (3, 5, 4)
    >>> f.nit
    2

    By default, `preserve_shape` is False, and therefore the callable
    `f` may be called with arrays of any broadcastable shapes.
    For example:

    >>> shapes = []
    >>> def f(x, c):
    ...    shape = np.broadcast_shapes(x.shape, c.shape)
    ...    shapes.append(shape)
    ...    return np.sin(c*x)
    >>>
    >>> c = [1, 5, 10, 20]
    >>> res = _differentiate(f, 0, args=(c,))
    >>> shapes
    [(4,), (4, 8), (4, 2), (3, 2), (2, 2), (1, 2)]

    To understand where these shapes are coming from - and to better
    understand how `_differentiate` computes accurate results - note that
    higher values of ``c`` correspond with higher frequency sinusoids.
    The higher frequency sinusoids make the function's derivative change
    faster, so more function evaluations are required to achieve the target
    accuracy:

    >>> res.nfev
    array([11, 13, 15, 17])

    The initial ``shape``, ``(4,)``, corresponds with evaluating the
    function at a single abscissa and all four frequencies; this is used
    for input validation and to determine the size and dtype of the arrays
    that store results. The next shape corresponds with evaluating the
    function at an initial grid of abscissae and all four frequencies.
    Successive calls to the function evaluate the function at two more
    abscissae, increasing the effective order of the approximation by two.
    However, in later function evaluations, the function is evaluated at
    fewer frequencies because the corresponding derivative has already
    converged to the required tolerance. This saves function evaluations to
    improve performance, but it requires the function to accept arguments of
    any shape.

    "Vector-valued" functions are unlikely to satisfy this requirement.
    For example, consider

    >>> def f(x):
    ...    return [x, np.sin(3*x), x+np.sin(10*x), np.sin(20*x)*(x-1)**2]

    This integrand is not compatible with `_differentiate` as written; for instance,
    the shape of the output will not be the same as the shape of ``x``. Such a
    function *could* be converted to a compatible form with the introduction of
    additional parameters, but this would be inconvenient. In such cases,
    a simpler solution would be to use `preserve_shape`.

    >>> shapes = []
    >>> def f(x):
    ...     shapes.append(x.shape)
    ...     x0, x1, x2, x3 = x
    ...     return [x0, np.sin(3*x1), x2+np.sin(10*x2), np.sin(20*x3)*(x3-1)**2]
    >>>
    >>> x = np.zeros(4)
    >>> res = _differentiate(f, x, preserve_shape=True)
    >>> shapes
    [(4,), (4, 8), (4, 2), (4, 2), (4, 2), (4, 2)]

    Here, the shape of ``x`` is ``(4,)``. With ``preserve_shape=True``, the
    function may be called with argument ``x`` of shape ``(4,)`` or ``(4, n)``,
    and this is what we observe.

    """
    # TODO (followup):
    #  - investigate behavior at saddle points
    #  - array initial_step / step_factor?
    #  - multivariate functions?

    res = _differentiate_iv(func, x, args, atol, rtol, maxiter, order, initial_step,
                            step_factor, step_direction, preserve_shape, callback)
    (func, x, args, atol, rtol, maxiter, order,
     h0, fac, hdir, preserve_shape, callback) = res

    # Initialization
    # Since f(x) (no step) is not needed for central differences, it may be
    # possible to eliminate this function evaluation. However, it's useful for
    # input validation and standardization, and everything else is designed to
    # reduce function calls, so let's keep it simple.
    temp = eim._initialize(func, (x,), args, preserve_shape=preserve_shape)
    func, xs, fs, args, shape, dtype, xp = temp
    x, f = xs[0], fs[0]
    df = np.full_like(f, np.nan)
    # Ideally we'd broadcast the shape of `hdir` in `_elementwise_algo_init`, but
    # it's simpler to do it here than to generalize `_elementwise_algo_init` further.
    # `hdir` and `x` are already broadcasted in `_differentiate_iv`, so we know
    # that `hdir` can be broadcasted to the final shape.
    hdir = np.broadcast_to(hdir, shape).flatten()

    status = np.full_like(x, eim._EINPROGRESS, dtype=int)  # in progress
    nit, nfev = 0, 1  # one function evaluations performed above
    # Boolean indices of left, central, right, and (all) one-sided steps
    il = hdir < 0
    ic = hdir == 0
    ir = hdir > 0
    io = il | ir

    # Most of these attributes are reasonably obvious, but:
    # - `fs` holds all the function values of all active `x`. The zeroth
    #   axis corresponds with active points `x`, the first axis corresponds
    #   with the different steps (in the order described in
    #   `_differentiate_weights`).
    # - `terms` (which could probably use a better name) is half the `order`,
    #   which is always even.
    work = _RichResult(x=x, df=df, fs=f[:, np.newaxis], error=np.nan, h=h0,
                       df_last=np.nan, error_last=np.nan, h0=h0, fac=fac,
                       atol=atol, rtol=rtol, nit=nit, nfev=nfev,
                       status=status, dtype=dtype, terms=(order+1)//2,
                       hdir=hdir, il=il, ic=ic, ir=ir, io=io)
    # This is the correspondence between terms in the `work` object and the
    # final result. In this case, the mapping is trivial. Note that `success`
    # is prepended automatically.
    res_work_pairs = [('status', 'status'), ('df', 'df'), ('error', 'error'),
                      ('nit', 'nit'), ('nfev', 'nfev'), ('x', 'x')]

    def pre_func_eval(work):
        """Determine the abscissae at which the function needs to be evaluated.

        See `_differentiate_weights` for a description of the stencil (pattern
        of the abscissae).

        In the first iteration, there is only one stored function value in
        `work.fs`, `f(x)`, so we need to evaluate at `order` new points. In
        subsequent iterations, we evaluate at two new points. Note that
        `work.x` is always flattened into a 1D array after broadcasting with
        all `args`, so we add a new axis at the end and evaluate all point
        in one call to the function.

        For improvement:
        - Consider measuring the step size actually taken, since `(x + h) - x`
          is not identically equal to `h` with floating point arithmetic.
        - Adjust the step size automatically if `x` is too big to resolve the
          step.
        - We could probably save some work if there are no central difference
          steps or no one-sided steps.
        """
        n = work.terms  # half the order
        h = work.h  # step size
        c = work.fac  # step reduction factor
        d = c**0.5  # square root of step reduction factor (one-sided stencil)
        # Note - no need to be careful about dtypes until we allocate `x_eval`

        if work.nit == 0:
            hc = h / c**np.arange(n)
            hc = np.concatenate((-hc[::-1], hc))
        else:
            hc = np.asarray([-h, h]) / c**(n-1)

        if work.nit == 0:
            hr = h / d**np.arange(2*n)
        else:
            hr = np.asarray([h, h/d]) / c**(n-1)

        n_new = 2*n if work.nit == 0 else 2  # number of new abscissae
        x_eval = np.zeros((len(work.hdir), n_new), dtype=work.dtype)
        il, ic, ir = work.il, work.ic, work.ir
        x_eval[ir] = work.x[ir, np.newaxis] + hr
        x_eval[ic] = work.x[ic, np.newaxis] + hc
        x_eval[il] = work.x[il, np.newaxis] - hr
        return x_eval

    def post_func_eval(x, f, work):
        """ Estimate the derivative and error from the function evaluations

        As in `pre_func_eval`: in the first iteration, there is only one stored
        function value in `work.fs`, `f(x)`, so we need to add the `order` new
        points. In subsequent iterations, we add two new points. The tricky
        part is getting the order to match that of the weights, which is
        described in `_differentiate_weights`.

        For improvement:
        - Change the order of the weights (and steps in `pre_func_eval`) to
          simplify `work_fc` concatenation and eliminate `fc` concatenation.
        - It would be simple to do one-step Richardson extrapolation with `df`
          and `df_last` to increase the order of the estimate and/or improve
          the error estimate.
        - Process the function evaluations in a more numerically favorable
          way. For instance, combining the pairs of central difference evals
          into a second-order approximation and using Richardson extrapolation
          to produce a higher order approximation seemed to retain accuracy up
          to very high order.
        - Alternatively, we could use `polyfit` like Jacobi. An advantage of
          fitting polynomial to more points than necessary is improved noise
          tolerance.
        """
        n = work.terms
        n_new = n if work.nit == 0 else 1
        il, ic, io = work.il, work.ic, work.io

        # Central difference
        # `work_fc` is *all* the points at which the function has been evaluated
        # `fc` is the points we're using *this iteration* to produce the estimate
        work_fc = (f[ic, :n_new], work.fs[ic, :], f[ic, -n_new:])
        work_fc = np.concatenate(work_fc, axis=-1)
        if work.nit == 0:
            fc = work_fc
        else:
            fc = (work_fc[:, :n], work_fc[:, n:n+1], work_fc[:, -n:])
            fc = np.concatenate(fc, axis=-1)

        # One-sided difference
        work_fo = np.concatenate((work.fs[io, :], f[io, :]), axis=-1)
        if work.nit == 0:
            fo = work_fo
        else:
            fo = np.concatenate((work_fo[:, 0:1], work_fo[:, -2*n:]), axis=-1)

        work.fs = np.zeros((len(ic), work.fs.shape[-1] + 2*n_new))
        work.fs[ic] = work_fc
        work.fs[io] = work_fo

        wc, wo = _differentiate_weights(work, n)
        work.df_last = work.df.copy()
        work.df[ic] = fc @ wc / work.h
        work.df[io] = fo @ wo / work.h
        work.df[il] *= -1

        work.h /= work.fac
        work.error_last = work.error
        # Simple error estimate - the difference in derivative estimates between
        # this iteration and the last. This is typically conservative because if
        # convergence has begin, the true error is much closer to the difference
        # between the current estimate and the *next* error estimate. However,
        # we could use Richarson extrapolation to produce an error estimate that
        # is one order higher, and take the difference between that and
        # `work.df` (which would just be constant factor that depends on `fac`.)
        work.error = abs(work.df - work.df_last)

    def check_termination(work):
        """Terminate due to convergence, non-finite values, or error increase"""
        stop = np.zeros_like(work.df).astype(bool)

        i = work.error < work.atol + work.rtol*abs(work.df)
        work.status[i] = eim._ECONVERGED
        stop[i] = True

        if work.nit > 0:
            i = ~((np.isfinite(work.x) & np.isfinite(work.df)) | stop)
            work.df[i], work.status[i] = np.nan, eim._EVALUEERR
            stop[i] = True

        # With infinite precision, there is a step size below which
        # all smaller step sizes will reduce the error. But in floating point
        # arithmetic, catastrophic cancellation will begin to cause the error
        # to increase again. This heuristic tries to avoid step sizes that are
        # too small. There may be more theoretically sound approaches for
        # detecting a step size that minimizes the total error, but this
        # heuristic seems simple and effective.
        i = (work.error > work.error_last*10) & ~stop
        work.status[i] = _EERRORINCREASE
        stop[i] = True

        return stop

    def post_termination_check(work):
        return

    def customize_result(res, shape):
        return shape

    return eim._loop(work, callback, shape, maxiter, func, args, dtype,
                     pre_func_eval, post_func_eval, check_termination,
                     post_termination_check, customize_result, res_work_pairs,
                     xp, preserve_shape)


def _differentiate_weights(work, n):
    # This produces the weights of the finite difference formula for a given
    # stencil. In experiments, use of a second-order central difference formula
    # with Richardson extrapolation was more accurate numerically, but it was
    # more complicated, and it would have become even more complicated when
    # adding support for one-sided differences. However, now that all the
    # function evaluation values are stored, they can be processed in whatever
    # way is desired to produce the derivative estimate. We leave alternative
    # approaches to future work. To be more self-contained, here is the theory
    # for deriving the weights below.
    #
    # Recall that the Taylor expansion of a univariate, scalar-values function
    # about a point `x` may be expressed as:
    #      f(x + h)  =     f(x) + f'(x)*h + f''(x)/2!*h**2  + O(h**3)
    # Suppose we evaluate f(x), f(x+h), and f(x-h).  We have:
    #      f(x)      =     f(x)
    #      f(x + h)  =     f(x) + f'(x)*h + f''(x)/2!*h**2  + O(h**3)
    #      f(x - h)  =     f(x) - f'(x)*h + f''(x)/2!*h**2  + O(h**3)
    # We can solve for weights `wi` such that:
    #   w1*f(x)      = w1*(f(x))
    # + w2*f(x + h)  = w2*(f(x) + f'(x)*h + f''(x)/2!*h**2) + O(h**3)
    # + w3*f(x - h)  = w3*(f(x) - f'(x)*h + f''(x)/2!*h**2) + O(h**3)
    #                =     0    + f'(x)*h + 0               + O(h**3)
    # Then
    #     f'(x) ~ (w1*f(x) + w2*f(x+h) + w3*f(x-h))/h
    # is a finite difference derivative approximation with error O(h**2),
    # and so it is said to be a "second-order" approximation. Under certain
    # conditions (e.g. well-behaved function, `h` sufficiently small), the
    # error in the approximation will decrease with h**2; that is, if `h` is
    # reduced by a factor of 2, the error is reduced by a factor of 4.
    #
    # By default, we use eighth-order formulae. Our central-difference formula
    # uses abscissae:
    #   x-h/c**3, x-h/c**2, x-h/c, x-h, x, x+h, x+h/c, x+h/c**2, x+h/c**3
    # where `c` is the step factor. (Typically, the step factor is greater than
    # one, so the outermost points - as written above - are actually closest to
    # `x`.) This "stencil" is chosen so that each iteration, the step can be
    # reduced by the factor `c`, and most of the function evaluations can be
    # reused with the new step size. For example, in the next iteration, we
    # will have:
    #   x-h/c**4, x-h/c**3, x-h/c**2, x-h/c, x, x+h/c, x+h/c**2, x+h/c**3, x+h/c**4
    # We do not reuse `x-h` and `x+h` for the new derivative estimate.
    # While this would increase the order of the formula and thus the
    # theoretical convergence rate, it is also less stable numerically.
    # (As noted above, there are other ways of processing the values that are
    # more stable. Thus, even now we store `f(x-h)` and `f(x+h)` in `work.fs`
    # to simplify future development of this sort of improvement.)
    #
    # The (right) one-sided formula is produced similarly using abscissae
    #   x, x+h, x+h/d, x+h/d**2, ..., x+h/d**6, x+h/d**7, x+h/d**7
    # where `d` is the square root of `c`. (The left one-sided formula simply
    # uses -h.) When the step size is reduced by factor `c = d**2`, we have
    # abscissae:
    #   x, x+h/d**2, x+h/d**3..., x+h/d**8, x+h/d**9, x+h/d**9
    # `d` is chosen as the square root of `c` so that the rate of the step-size
    # reduction is the same per iteration as in the central difference case.
    # Note that because the central difference formulas are inherently of even
    # order, for simplicity, we use only even-order formulas for one-sided
    # differences, too.

    # It's possible for the user to specify `fac` in, say, double precision but
    # `x` and `args` in single precision. `fac` gets converted to single
    # precision, but we should always use double precision for the intermediate
    # calculations here to avoid additional error in the weights.
    fac = work.fac.astype(np.float64)

    # Note that if the user switches back to floating point precision with
    # `x` and `args`, then `fac` will not necessarily equal the (lower
    # precision) cached `_differentiate_weights.fac`, and the weights will
    # need to be recalculated. This could be fixed, but it's late, and of
    # low consequence.
    if fac != _differentiate_weights.fac:
        _differentiate_weights.central = []
        _differentiate_weights.right = []
        _differentiate_weights.fac = fac

    if len(_differentiate_weights.central) != 2*n + 1:
        # Central difference weights. Consider refactoring this; it could
        # probably be more compact.
        i = np.arange(-n, n + 1)
        p = np.abs(i) - 1.  # center point has power `p` -1, but sign `s` is 0
        s = np.sign(i)

        h = s / fac ** p
        A = np.vander(h, increasing=True).T
        b = np.zeros(2*n + 1)
        b[1] = 1
        weights = np.linalg.solve(A, b)

        # Enforce identities to improve accuracy
        weights[n] = 0
        for i in range(n):
            weights[-i-1] = -weights[i]

        # Cache the weights. We only need to calculate them once unless
        # the step factor changes.
        _differentiate_weights.central = weights

        # One-sided difference weights. The left one-sided weights (with
        # negative steps) are simply the negative of the right one-sided
        # weights, so no need to compute them separately.
        i = np.arange(2*n + 1)
        p = i - 1.
        s = np.sign(i)

        h = s / np.sqrt(fac) ** p
        A = np.vander(h, increasing=True).T
        b = np.zeros(2 * n + 1)
        b[1] = 1
        weights = np.linalg.solve(A, b)

        _differentiate_weights.right = weights

    return (_differentiate_weights.central.astype(work.dtype, copy=False),
            _differentiate_weights.right.astype(work.dtype, copy=False))
_differentiate_weights.central = []
_differentiate_weights.right = []
_differentiate_weights.fac = None


def _jacobian(func, x, *, atol=None, rtol=None, maxiter=10,
              order=8, initial_step=0.5, step_factor=2.0):
    r"""Evaluate the Jacobian of a function numerically.

    Parameters
    ----------
    func : callable
        The function whose Jacobian is desired. The signature must be::

            func(x: ndarray) -> ndarray

         where each element of ``x`` is a finite real. If the function to be
         differentiated accepts additional, arguments wrap it (e.g. using
         `functools.partial` or ``lambda``) and pass the wrapped callable
         into `_jacobian`. See Notes regarding vectorization and the dimensionality
         of the input and output.
    x : array_like
        Points at which to evaluate the Jacobian. Must have at least one dimension.
        See Notes regarding the dimensionality and vectorization.
    atol, rtol : float, optional
        Absolute and relative tolerances for the stopping condition: iteration
        will stop for each element of the Jacobian when
        ``res.error < atol + rtol * abs(res.df)``. The default `atol` is the
        smallest normal number of the appropriate dtype, and the default `rtol`
        is the square root of the precision of the appropriate dtype.
    order : int, default: 8
        The (positive integer) order of the finite difference formula to be
        used. Odd integers will be rounded up to the next even integer.
    initial_step : float, default: 0.5
        The (absolute) initial step size for the finite difference derivative
        approximation.
    step_factor : float, default: 2.0
        The factor by which the step size is *reduced* in each iteration; i.e.
        the step size in iteration 1 is ``initial_step/step_factor``. If
        ``step_factor < 1``, subsequent steps will be greater than the initial
        step; this may be useful if steps smaller than some threshold are
        undesirable (e.g. due to subtractive cancellation error).
    maxiter : int, default: 10
        The maximum number of iterations of the algorithm to perform.

    Returns
    -------
    res : _RichResult
        An instance of `scipy._lib._util._RichResult` with the following
        attributes.

        success : bool array
            ``True`` when the algorithm terminated successfully (status ``0``).
        status : int array
            An integer representing the exit status of the algorithm.
            ``0`` : The algorithm converged to the specified tolerances.
            ``-1`` : The error estimate increased, so iteration was terminated.
            ``-2`` : The maximum number of iterations was reached.
            ``-3`` : A non-finite value was encountered.
            ``-4`` : Iteration was terminated by `callback`.
            ``1`` : The algorithm is proceeding normally (in `callback` only).
        df : float array
            The Jacobian of `func` at `x`, if the algorithm terminated
            successfully.
        error : float array
            An estimate of the error: the magnitude of the difference between
            the current estimate of the derivative and the estimate in the
            previous iteration.
        nit : int array
            The number of iterations performed.
        nfev : int array
            The number of points at which `func` was evaluated.
        x : float array
            The value at which the derivative of `func` was evaluated.

    See Also
    --------
    _differentiate

    Notes
    -----
    Suppose we wish to evaluate the Jacobian of a function
    :math:`f: \mathbf{R^m} \rightarrow \mathbf{R^n}`, and assign to variables
    ``m`` and ``n`` the positive integer values of :math:`m` and :math:`n`,
    respectively. If we wish to evaluate the Jacobian at a single point,
    then:

    - argument `x` must be an array of shape ``(m,)``
    - argument `func` must be vectorized to accept an array of shape ``(m, p)``.
      The first axis represents the :math:`m` inputs of :math:`f`; the second
      is for evaluating the function at multiple points in a single call.
    - argument `func` must return an array of shape ``(n, p)``. The first
      axis represents the :math:`n` outputs of :math:`f`; the second
      is for the result of evaluating the function at multiple points.
    - attribute ``df`` of the result object will be an array of shape ``(n, m)``,
      the Jacobian.

    This function is also vectorized in the sense that the Jacobian can be
    evaluated at ``k`` points in a single call. In this case, `x` would be an
    array of shape ``(m, k)``, `func` would accept an array of shape
    ``(m, k, p)`` and return an array of shape ``(n, k, p)``, and the ``df``
    attribute of the result would have shape ``(n, m, k)``.

    References
    ----------
    .. [1] Jacobian matrix and determinant, *Wikipedia*,
           https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant

    Examples
    --------
    The Rosenbrock function maps from :math:`\mathbf{R}^m \righarrow \mathbf{R}`;
    the SciPy implementation `scipy.optimize.rosen` is vectorized to accept an
    array of shape ``(m, p)`` and return an array of shape ``m``. Suppose we wish
    to evaluate the Jacobian (AKA the gradient because the function returns a scalar)
    at ``[0.5, 0.5, 0.5]``.

    >>> import numpy as np
    >>> from scipy.optimize._differentiate import _jacobian as jacobian
    >>> from scipy.optimize import rosen, rosen_der
    >>> m = 3
    >>> x = np.full(m, 0.5)
    >>> res = jacobian(rosen, x)
    >>> ref = rosen_der(x)  # reference value of the gradient
    >>> res.df, ref
    (array([-51.,  -1.,  50.]), array([-51.,  -1.,  50.]))

    As an example of a function with multiple outputs, consider Example 4
    from [1]_.

    >>> def f(x):
    ...     x1, x2, x3 = x    ...
    ...     return [x1, 5*x3, 4*x2**2 - 2*x3, x3*np.sin(x1)]

    The true Jacobian is given by:

    >>> def df(x):
    ...         x1, x2, x3 = x
    ...         one = np.ones_like(x1)
    ...         return [[one, 0*one, 0*one],
    ...                 [0*one, 0*one, 5*one],
    ...                 [0*one, 8*x2, -2*one],
    ...                 [x3*np.cos(x1), 0*one, np.sin(x1)]]

    Evaluate the Jacobian at an arbitrary point.

    >>> rng = np.random.default_rng(389252938452)
    >>> x = rng.random(size=3)
    >>> res = jacobian(f, x)
    >>> ref = df(x)
    >>> res.df.shape == (4, 3)
    True
    >>> np.allclose(res.df, ref)
    True

    Evaluate the Jacobian at 10 arbitrary points in a single call.

    >>> x = rng.random(size=(3, 10))
    >>> res = jacobian(f, x)
    >>> ref = df(x)
    >>> res.df.shape == (4, 3, 10)
    True
    >>> np.allclose(res.df, ref)
    True

    """
    x = np.asarray(x)
    int_dtype = np.issubdtype(x.dtype, np.integer)
    x0 = np.asarray(x, dtype=float) if int_dtype else x

    if x0.ndim < 1:
        message = "Argument `x` must be at least 1-D."
        raise ValueError(message)

    m = x0.shape[0]
    i = np.arange(m)

    def wrapped(x):
        p = () if x.ndim == x0.ndim else (x.shape[-1],)  # number of abscissae
        new_dims = (1,) if x.ndim == x0.ndim else (1, -1)
        new_shape = (m, m) + x0.shape[1:] + p
        xph = np.expand_dims(x0, new_dims)
        xph = np.broadcast_to(xph, new_shape).copy()
        xph[i, i] = x
        return func(xph)

    res = _differentiate(wrapped, x, atol=atol, rtol=rtol,
                         maxiter=maxiter, order=order, initial_step=initial_step,
                         step_factor=step_factor, preserve_shape=True)
    del res.x  # the user knows `x`, and the way it gets broadcasted is meaningless here
    return res
