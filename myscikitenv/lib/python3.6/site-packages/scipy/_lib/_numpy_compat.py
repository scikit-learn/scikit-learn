"""Functions copypasted from newer versions of numpy.

"""
from __future__ import division, print_function, absolute_import

import warnings
from warnings import WarningMessage
import re
from functools import wraps
import numpy as np

from scipy._lib._version import NumpyVersion


if NumpyVersion(np.__version__) > '1.7.0.dev':
    _assert_warns = np.testing.assert_warns
else:
    def _assert_warns(warning_class, func, *args, **kw):
        r"""
        Fail unless the given callable throws the specified warning.

        This definition is copypasted from numpy 1.9.0.dev.
        The version in earlier numpy returns None.

        Parameters
        ----------
        warning_class : class
            The class defining the warning that `func` is expected to throw.
        func : callable
            The callable to test.
        *args : Arguments
            Arguments passed to `func`.
        **kwargs : Kwargs
            Keyword arguments passed to `func`.

        Returns
        -------
        The value returned by `func`.

        """
        with warnings.catch_warnings(record=True) as l:
            warnings.simplefilter('always')
            result = func(*args, **kw)
            if not len(l) > 0:
                raise AssertionError("No warning raised when calling %s"
                        % func.__name__)
            if not l[0].category is warning_class:
                raise AssertionError("First warning for %s is not a "
                        "%s( is %s)" % (func.__name__, warning_class, l[0]))
        return result


if NumpyVersion(np.__version__) >= '1.10.0':
    from numpy import broadcast_to
else:
    # Definition of `broadcast_to` from numpy 1.10.0.

    def _maybe_view_as_subclass(original_array, new_array):
        if type(original_array) is not type(new_array):
            # if input was an ndarray subclass and subclasses were OK,
            # then view the result as that subclass.
            new_array = new_array.view(type=type(original_array))
            # Since we have done something akin to a view from original_array, we
            # should let the subclass finalize (if it has it implemented, i.e., is
            # not None).
            if new_array.__array_finalize__:
                new_array.__array_finalize__(original_array)
        return new_array

    def _broadcast_to(array, shape, subok, readonly):
        shape = tuple(shape) if np.iterable(shape) else (shape,)
        array = np.array(array, copy=False, subok=subok)
        if not shape and array.shape:
            raise ValueError('cannot broadcast a non-scalar to a scalar array')
        if any(size < 0 for size in shape):
            raise ValueError('all elements of broadcast shape must be non-'
                             'negative')
        broadcast = np.nditer(
            (array,), flags=['multi_index', 'refs_ok', 'zerosize_ok'],
            op_flags=['readonly'], itershape=shape, order='C').itviews[0]
        result = _maybe_view_as_subclass(array, broadcast)
        if not readonly and array.flags.writeable:
            result.flags.writeable = True
        return result

    def broadcast_to(array, shape, subok=False):
        return _broadcast_to(array, shape, subok=subok, readonly=True)


if NumpyVersion(np.__version__) >= '1.11.0':
    def get_randint(random_state):
        return random_state.randint
else:
    # In NumPy versions previous to 1.11.0 the randint funtion and the randint
    # method of RandomState does only work with int32 values.
    def get_randint(random_state):
        def randint_patched(low, high, size, dtype=np.int32):
            low = max(low, np.iinfo(dtype).min, np.iinfo(np.int32).min)
            high = min(high, np.iinfo(dtype).max, np.iinfo(np.int32).max)
            integers = random_state.randint(low, high=high, size=size)
            return integers.astype(dtype, copy=False)
        return randint_patched


if NumpyVersion(np.__version__) >= '1.9.0':
    from numpy import unique
else:
    # the return_counts keyword was added in 1.9.0
    def unique(ar, return_index=False, return_inverse=False, return_counts=False):
        """
        Find the unique elements of an array.

        Returns the sorted unique elements of an array. There are three optional
        outputs in addition to the unique elements: the indices of the input array
        that give the unique values, the indices of the unique array that
        reconstruct the input array, and the number of times each unique value
        comes up in the input array.

        Parameters
        ----------
        ar : array_like
            Input array. This will be flattened if it is not already 1-D.
        return_index : bool, optional
            If True, also return the indices of `ar` that result in the unique
            array.
        return_inverse : bool, optional
            If True, also return the indices of the unique array that can be used
            to reconstruct `ar`.
        return_counts : bool, optional
            If True, also return the number of times each unique value comes up
            in `ar`.

            .. versionadded:: 1.9.0

        Returns
        -------
        unique : ndarray
            The sorted unique values.
        unique_indices : ndarray, optional
            The indices of the first occurrences of the unique values in the
            (flattened) original array. Only provided if `return_index` is True.
        unique_inverse : ndarray, optional
            The indices to reconstruct the (flattened) original array from the
            unique array. Only provided if `return_inverse` is True.
        unique_counts : ndarray, optional
            The number of times each of the unique values comes up in the
            original array. Only provided if `return_counts` is True.

            .. versionadded:: 1.9.0

        Notes
        -----
        Taken over from numpy 1.12.0-dev (c8408bf9c).  Omitted examples,
        see numpy documentation for those.

        """
        ar = np.asanyarray(ar).flatten()

        optional_indices = return_index or return_inverse
        optional_returns = optional_indices or return_counts

        if ar.size == 0:
            if not optional_returns:
                ret = ar
            else:
                ret = (ar,)
                if return_index:
                    ret += (np.empty(0, np.bool),)
                if return_inverse:
                    ret += (np.empty(0, np.bool),)
                if return_counts:
                    ret += (np.empty(0, np.intp),)
            return ret

        if optional_indices:
            perm = ar.argsort(kind='mergesort' if return_index else 'quicksort')
            aux = ar[perm]
        else:
            ar.sort()
            aux = ar
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))

        if not optional_returns:
            ret = aux[flag]
        else:
            ret = (aux[flag],)
            if return_index:
                ret += (perm[flag],)
            if return_inverse:
                iflag = np.cumsum(flag) - 1
                inv_idx = np.empty(ar.shape, dtype=np.intp)
                inv_idx[perm] = iflag
                ret += (inv_idx,)
            if return_counts:
                idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
                ret += (np.diff(idx),)
        return ret


if NumpyVersion(np.__version__) > '1.12.0.dev':
    polyvalfromroots = np.polynomial.polynomial.polyvalfromroots
else:
    def polyvalfromroots(x, r, tensor=True):
        r"""
        Evaluate a polynomial specified by its roots at points x.

        This function is copypasted from numpy 1.12.0.dev.

        If `r` is of length `N`, this function returns the value

        .. math:: p(x) = \prod_{n=1}^{N} (x - r_n)

        The parameter `x` is converted to an array only if it is a tuple or a
        list, otherwise it is treated as a scalar. In either case, either `x`
        or its elements must support multiplication and addition both with
        themselves and with the elements of `r`.

        If `r` is a 1-D array, then `p(x)` will have the same shape as `x`.  If
        `r` is multidimensional, then the shape of the result depends on the
        value of `tensor`. If `tensor is ``True`` the shape will be r.shape[1:]
        + x.shape; that is, each polynomial is evaluated at every value of `x`.
        If `tensor` is ``False``, the shape will be r.shape[1:]; that is, each
        polynomial is evaluated only for the corresponding broadcast value of
        `x`. Note that scalars have shape (,).

        Parameters
        ----------
        x : array_like, compatible object
            If `x` is a list or tuple, it is converted to an ndarray, otherwise
            it is left unchanged and treated as a scalar. In either case, `x`
            or its elements must support addition and multiplication with with
            themselves and with the elements of `r`.
        r : array_like
            Array of roots. If `r` is multidimensional the first index is the
            root index, while the remaining indices enumerate multiple
            polynomials. For instance, in the two dimensional case the roots of
            each polynomial may be thought of as stored in the columns of `r`.
        tensor : boolean, optional
            If True, the shape of the roots array is extended with ones on the
            right, one for each dimension of `x`. Scalars have dimension 0 for
            this action. The result is that every column of coefficients in `r`
            is evaluated for every element of `x`. If False, `x` is broadcast
            over the columns of `r` for the evaluation.  This keyword is useful
            when `r` is multidimensional. The default value is True.

        Returns
        -------
        values : ndarray, compatible object
            The shape of the returned array is described above.

        See Also
        --------
        polyroots, polyfromroots, polyval

        Examples
        --------
        >>> from numpy.polynomial.polynomial import polyvalfromroots
        >>> polyvalfromroots(1, [1,2,3])
        0.0
        >>> a = np.arange(4).reshape(2,2)
        >>> a
        array([[0, 1],
               [2, 3]])
        >>> polyvalfromroots(a, [-1, 0, 1])
        array([[ -0.,   0.],
               [  6.,  24.]])
        >>> r = np.arange(-2, 2).reshape(2,2) # multidimensional coefficients
        >>> r # each column of r defines one polynomial
        array([[-2, -1],
               [ 0,  1]])
        >>> b = [-2, 1]
        >>> polyvalfromroots(b, r, tensor=True)
        array([[-0.,  3.],
               [ 3., 0.]])
        >>> polyvalfromroots(b, r, tensor=False)
        array([-0.,  0.])
        """
        r = np.array(r, ndmin=1, copy=0)
        if r.dtype.char in '?bBhHiIlLqQpP':
            r = r.astype(np.double)
        if isinstance(x, (tuple, list)):
            x = np.asarray(x)
        if isinstance(x, np.ndarray):
            if tensor:
                r = r.reshape(r.shape + (1,)*x.ndim)
            elif x.ndim >= r.ndim:
                raise ValueError("x.ndim must be < r.ndim when tensor == "
                                 "False")
        return np.prod(x - r, axis=0)


try:
    from numpy.testing import suppress_warnings
except ImportError:
    class suppress_warnings(object):
        """
        Context manager and decorator doing much the same as
        ``warnings.catch_warnings``.

        However, it also provides a filter mechanism to work around
        https://bugs.python.org/issue4180.

        This bug causes Python before 3.4 to not reliably show warnings again
        after they have been ignored once (even within catch_warnings). It
        means that no "ignore" filter can be used easily, since following
        tests might need to see the warning. Additionally it allows easier
        specificity for testing warnings and can be nested.

        Parameters
        ----------
        forwarding_rule : str, optional
            One of "always", "once", "module", or "location". Analogous to
            the usual warnings module filter mode, it is useful to reduce
            noise mostly on the outmost level. Unsuppressed and unrecorded
            warnings will be forwarded based on this rule. Defaults to "always".
            "location" is equivalent to the warnings "default", match by exact
            location the warning warning originated from.

        Notes
        -----
        Filters added inside the context manager will be discarded again
        when leaving it. Upon entering all filters defined outside a
        context will be applied automatically.

        When a recording filter is added, matching warnings are stored in the
        ``log`` attribute as well as in the list returned by ``record``.

        If filters are added and the ``module`` keyword is given, the
        warning registry of this module will additionally be cleared when
        applying it, entering the context, or exiting it. This could cause
        warnings to appear a second time after leaving the context if they
        were configured to be printed once (default) and were already
        printed before the context was entered.

        Nesting this context manager will work as expected when the
        forwarding rule is "always" (default). Unfiltered and unrecorded
        warnings will be passed out and be matched by the outer level.
        On the outmost level they will be printed (or caught by another
        warnings context). The forwarding rule argument can modify this
        behaviour.

        Like ``catch_warnings`` this context manager is not threadsafe.

        Examples
        --------
        >>> with suppress_warnings() as sup:
        ...     sup.filter(DeprecationWarning, "Some text")
        ...     sup.filter(module=np.ma.core)
        ...     log = sup.record(FutureWarning, "Does this occur?")
        ...     command_giving_warnings()
        ...     # The FutureWarning was given once, the filtered warnings were
        ...     # ignored. All other warnings abide outside settings (may be
        ...     # printed/error)
        ...     assert_(len(log) == 1)
        ...     assert_(len(sup.log) == 1)  # also stored in log attribute

        Or as a decorator:

        >>> sup = suppress_warnings()
        >>> sup.filter(module=np.ma.core)  # module must match exact
        >>> @sup
        >>> def some_function():
        ...     # do something which causes a warning in np.ma.core
        ...     pass
        """
        def __init__(self, forwarding_rule="always"):
            self._entered = False

            # Suppressions are either instance or defined inside one with block:
            self._suppressions = []

            if forwarding_rule not in {"always", "module", "once", "location"}:
                raise ValueError("unsupported forwarding rule.")
            self._forwarding_rule = forwarding_rule

        def _clear_registries(self):
            if hasattr(warnings, "_filters_mutated"):
                # clearing the registry should not be necessary on new pythons,
                # instead the filters should be mutated.
                warnings._filters_mutated()
                return
            # Simply clear the registry, this should normally be harmless,
            # note that on new pythons it would be invalidated anyway.
            for module in self._tmp_modules:
                if hasattr(module, "__warningregistry__"):
                    module.__warningregistry__.clear()

        def _filter(self, category=Warning, message="", module=None, record=False):
            if record:
                record = []  # The log where to store warnings
            else:
                record = None
            if self._entered:
                if module is None:
                    warnings.filterwarnings(
                        "always", category=category, message=message)
                else:
                    module_regex = module.__name__.replace('.', r'\.') + '$'
                    warnings.filterwarnings(
                        "always", category=category, message=message,
                        module=module_regex)
                    self._tmp_modules.add(module)
                    self._clear_registries()

                self._tmp_suppressions.append(
                    (category, message, re.compile(message, re.I), module, record))
            else:
                self._suppressions.append(
                    (category, message, re.compile(message, re.I), module, record))

            return record

        def filter(self, category=Warning, message="", module=None):
            """
            Add a new suppressing filter or apply it if the state is entered.

            Parameters
            ----------
            category : class, optional
                Warning class to filter
            message : string, optional
                Regular expression matching the warning message.
            module : module, optional
                Module to filter for. Note that the module (and its file)
                must match exactly and cannot be a submodule. This may make
                it unreliable for external modules.

            Notes
            -----
            When added within a context, filters are only added inside
            the context and will be forgotten when the context is exited.
            """
            self._filter(category=category, message=message, module=module,
                         record=False)

        def record(self, category=Warning, message="", module=None):
            """
            Append a new recording filter or apply it if the state is entered.

            All warnings matching will be appended to the ``log`` attribute.

            Parameters
            ----------
            category : class, optional
                Warning class to filter
            message : string, optional
                Regular expression matching the warning message.
            module : module, optional
                Module to filter for. Note that the module (and its file)
                must match exactly and cannot be a submodule. This may make
                it unreliable for external modules.

            Returns
            -------
            log : list
                A list which will be filled with all matched warnings.

            Notes
            -----
            When added within a context, filters are only added inside
            the context and will be forgotten when the context is exited.
            """
            return self._filter(category=category, message=message, module=module,
                                record=True)

        def __enter__(self):
            if self._entered:
                raise RuntimeError("cannot enter suppress_warnings twice.")

            self._orig_show = warnings.showwarning
            self._filters = warnings.filters
            warnings.filters = self._filters[:]

            self._entered = True
            self._tmp_suppressions = []
            self._tmp_modules = set()
            self._forwarded = set()

            self.log = []  # reset global log (no need to keep same list)

            for cat, mess, _, mod, log in self._suppressions:
                if log is not None:
                    del log[:]  # clear the log
                if mod is None:
                    warnings.filterwarnings(
                        "always", category=cat, message=mess)
                else:
                    module_regex = mod.__name__.replace('.', r'\.') + '$'
                    warnings.filterwarnings(
                        "always", category=cat, message=mess,
                        module=module_regex)
                    self._tmp_modules.add(mod)
            warnings.showwarning = self._showwarning
            self._clear_registries()

            return self

        def __exit__(self, *exc_info):
            warnings.showwarning = self._orig_show
            warnings.filters = self._filters
            self._clear_registries()
            self._entered = False
            del self._orig_show
            del self._filters

        def _showwarning(self, message, category, filename, lineno,
                         *args, **kwargs):
            use_warnmsg = kwargs.pop("use_warnmsg", None)
            for cat, _, pattern, mod, rec in (
                    self._suppressions + self._tmp_suppressions)[::-1]:
                if (issubclass(category, cat) and
                        pattern.match(message.args[0]) is not None):
                    if mod is None:
                        # Message and category match, either recorded or ignored
                        if rec is not None:
                            msg = WarningMessage(message, category, filename,
                                                 lineno, **kwargs)
                            self.log.append(msg)
                            rec.append(msg)
                        return
                    # Use startswith, because warnings strips the c or o from
                    # .pyc/.pyo files.
                    elif mod.__file__.startswith(filename):
                        # The message and module (filename) match
                        if rec is not None:
                            msg = WarningMessage(message, category, filename,
                                                 lineno, **kwargs)
                            self.log.append(msg)
                            rec.append(msg)
                        return

            # There is no filter in place, so pass to the outside handler
            # unless we should only pass it once
            if self._forwarding_rule == "always":
                if use_warnmsg is None:
                    self._orig_show(message, category, filename, lineno,
                                    *args, **kwargs)
                else:
                    self._orig_showmsg(use_warnmsg)
                return

            if self._forwarding_rule == "once":
                signature = (message.args, category)
            elif self._forwarding_rule == "module":
                signature = (message.args, category, filename)
            elif self._forwarding_rule == "location":
                signature = (message.args, category, filename, lineno)

            if signature in self._forwarded:
                return
            self._forwarded.add(signature)
            if use_warnmsg is None:
                self._orig_show(message, category, filename, lineno, *args,
                                **kwargs)
            else:
                self._orig_showmsg(use_warnmsg)

        def __call__(self, func):
            """
            Function decorator to apply certain suppressions to a whole
            function.
            """
            @wraps(func)
            def new_func(*args, **kwargs):
                with self:
                    return func(*args, **kwargs)

            return new_func

if NumpyVersion(np.__version__) >= '1.10.0':
    from numpy import cov
else:
    from numpy import array, average, dot

    def cov(m, y=None, rowvar=True, bias=False, ddof=None, fweights=None,
            aweights=None):
        """
        Estimate a covariance matrix, given data and weights.

        Covariance indicates the level to which two variables vary together.
        If we examine N-dimensional samples, :math:`X = [x_1, x_2, ... x_N]^T`,
        then the covariance matrix element :math:`C_{ij}` is the covariance of
        :math:`x_i` and :math:`x_j`. The element :math:`C_{ii}` is the variance
        of :math:`x_i`.

        See the notes for an outline of the algorithm.

        Parameters
        ----------
        m : array_like
            A 1-D or 2-D array containing multiple variables and observations.
            Each row of `m` represents a variable, and each column a single
            observation of all those variables. Also see `rowvar` below.
        y : array_like, optional
            An additional set of variables and observations. `y` has the same form
            as that of `m`.
        rowvar : bool, optional
            If `rowvar` is True (default), then each row represents a
            variable, with observations in the columns. Otherwise, the relationship
            is transposed: each column represents a variable, while the rows
            contain observations.
        bias : bool, optional
            Default normalization (False) is by ``(N - 1)``, where ``N`` is the
            number of observations given (unbiased estimate). If `bias` is True,
            then normalization is by ``N``. These values can be overridden by using
            the keyword ``ddof`` in numpy versions >= 1.5.
        ddof : int, optional
            If not ``None`` the default value implied by `bias` is overridden.
            Note that ``ddof=1`` will return the unbiased estimate, even if both
            `fweights` and `aweights` are specified, and ``ddof=0`` will return
            the simple average. See the notes for the details. The default value
            is ``None``.

            .. versionadded:: 1.5
        fweights : array_like, int, optional
            1-D array of integer freguency weights; the number of times each
            observation vector should be repeated.

            .. versionadded:: 1.10
        aweights : array_like, optional
            1-D array of observation vector weights. These relative weights are
            typically large for observations considered "important" and smaller for
            observations considered less "important". If ``ddof=0`` the array of
            weights can be used to assign probabilities to observation vectors.

            .. versionadded:: 1.10

        Returns
        -------
        out : ndarray
            The covariance matrix of the variables.

        See Also
        --------
        corrcoef : Normalized covariance matrix

        Notes
        -----
        Assume that the observations are in the columns of the observation
        array `m` and let ``f = fweights`` and ``a = aweights`` for brevity. The
        steps to compute the weighted covariance are as follows::

            >>> w = f * a
            >>> v1 = np.sum(w)
            >>> v2 = np.sum(w * a)
            >>> m -= np.sum(m * w, axis=1, keepdims=True) / v1
            >>> cov = np.dot(m * w, m.T) * v1 / (v1**2 - ddof * v2)

        Note that when ``a == 1``, the normalization factor
        ``v1 / (v1**2 - ddof * v2)`` goes over to ``1 / (np.sum(f) - ddof)``
        as it should.

        Examples
        --------
        Consider two variables, :math:`x_0` and :math:`x_1`, which
        correlate perfectly, but in opposite directions:

        >>> x = np.array([[0, 2], [1, 1], [2, 0]]).T
        >>> x
        array([[0, 1, 2],
               [2, 1, 0]])

        Note how :math:`x_0` increases while :math:`x_1` decreases. The covariance
        matrix shows this clearly:

        >>> np.cov(x)
        array([[ 1., -1.],
               [-1.,  1.]])

        Note that element :math:`C_{0,1}`, which shows the correlation between
        :math:`x_0` and :math:`x_1`, is negative.

        Further, note how `x` and `y` are combined:

        >>> x = [-2.1, -1,  4.3]
        >>> y = [3,  1.1,  0.12]
        >>> X = np.stack((x, y), axis=0)
        >>> print(np.cov(X))
        [[ 11.71        -4.286     ]
         [ -4.286        2.14413333]]
        >>> print(np.cov(x, y))
        [[ 11.71        -4.286     ]
         [ -4.286        2.14413333]]
        >>> print(np.cov(x))
        11.71

        """
        # Check inputs
        if ddof is not None and ddof != int(ddof):
            raise ValueError(
                "ddof must be integer")

        # Handles complex arrays too
        m = np.asarray(m)
        if m.ndim > 2:
            raise ValueError("m has more than 2 dimensions")

        if y is None:
            dtype = np.result_type(m, np.float64)
        else:
            y = np.asarray(y)
            if y.ndim > 2:
                raise ValueError("y has more than 2 dimensions")
            dtype = np.result_type(m, y, np.float64)

        X = array(m, ndmin=2, dtype=dtype)
        if not rowvar and X.shape[0] != 1:
            X = X.T
        if X.shape[0] == 0:
            return np.array([]).reshape(0, 0)
        if y is not None:
            y = array(y, copy=False, ndmin=2, dtype=dtype)
            if not rowvar and y.shape[0] != 1:
                y = y.T
            X = np.concatenate((X, y), axis=0)

        if ddof is None:
            if bias == 0:
                ddof = 1
            else:
                ddof = 0

        # Get the product of frequencies and weights
        w = None
        if fweights is not None:
            fweights = np.asarray(fweights, dtype=float)
            if not np.all(fweights == np.around(fweights)):
                raise TypeError(
                    "fweights must be integer")
            if fweights.ndim > 1:
                raise RuntimeError(
                    "cannot handle multidimensional fweights")
            if fweights.shape[0] != X.shape[1]:
                raise RuntimeError(
                    "incompatible numbers of samples and fweights")
            if any(fweights < 0):
                raise ValueError(
                    "fweights cannot be negative")
            w = fweights
        if aweights is not None:
            aweights = np.asarray(aweights, dtype=float)
            if aweights.ndim > 1:
                raise RuntimeError(
                    "cannot handle multidimensional aweights")
            if aweights.shape[0] != X.shape[1]:
                raise RuntimeError(
                    "incompatible numbers of samples and aweights")
            if any(aweights < 0):
                raise ValueError(
                    "aweights cannot be negative")
            if w is None:
                w = aweights
            else:
                w *= aweights

        avg, w_sum = average(X, axis=1, weights=w, returned=True)
        w_sum = w_sum[0]

        # Determine the normalization
        if w is None:
            fact = X.shape[1] - ddof
        elif ddof == 0:
            fact = w_sum
        elif aweights is None:
            fact = w_sum - ddof
        else:
            fact = w_sum - ddof*sum(w*aweights)/w_sum

        if fact <= 0:
            warnings.warn("Degrees of freedom <= 0 for slice",
                          RuntimeWarning, stacklevel=2)
            fact = 0.0

        X -= avg[:, None]
        if w is None:
            X_T = X.T
        else:
            X_T = (X*w).T
        c = dot(X, X_T.conj())
        c *= 1. / np.float64(fact)
        return c.squeeze()
