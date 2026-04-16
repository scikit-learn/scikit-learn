import functools
import operator
from collections.abc import Callable
from dataclasses import dataclass
from types import ModuleType

import numpy as np
from scipy._lib._array_api import (
    array_namespace, scipy_namespace_for, is_numpy, is_dask, is_marray,
    xp_promote, xp_capabilities, SCIPY_ARRAY_API, get_native_namespace_name
)
import scipy._lib.array_api_extra as xpx
from . import _basic
from . import _spfun_stats
from . import _ufuncs


@dataclass
class _FuncInfo:
    # NumPy-only function. IT MUST BE ELEMENTWISE.
    func: Callable
    # Number of arguments, not counting out=
    # This is for testing purposes only, due to the fact that
    # inspect.signature() just returns *args for ufuncs.
    n_args: int
    # @xp_capabilities decorator, for the purpose of
    # documentation and unit testing. Omit to indicate
    # full support for all backends.
    xp_capabilities: Callable[[Callable], Callable] | None = None
    # Generic implementation to fall back on if there is no native dispatch
    # available. This is a function that accepts (main namespace, scipy namespace)
    # and returns the final callable, or None if not available.
    generic_impl: Callable[
        [ModuleType, ModuleType | None], Callable | None
    ] | None = None
    # Handle case where a backend uses an alternative name for a function.
    # Should map backend names to alternative function names.
    alt_names_map: dict[str, str] | None = None
    # Some functions only take integer arrays for some arguments.
    int_only: tuple[bool] | None = None
    # For testing purposes, whether tests should only use positive values
    # for some arguments. If bool and equal to True, restrict to positive
    # values for all arguments. To restrict only some arguments to positive
    # values, pass a tuple of bool of the same length as the number of
    # arguments, the ith entry in the tuple controls positive_only for
    # the ith argument. To make backend specific choices for positive_only,
    # pass in a dict mapping backend names to bool or tuple[bool].
    positive_only: bool | tuple[bool] | dict[str, tuple[bool]] = False
    # Some special functions are not ufuncs and ufunc-specific tests
    # should not be applied to these.
    is_ufunc: bool = True
    # Some non-ufunc special functions take only Python ints for some arguments.
    # If so, python_int_only should be a tuple of the same length as the number
    # of arguments,with value True if the corresponding argument needs to be a
    # Python int.
    # Can also take a dict mapping backends to such tuples if an argument being
    # Python int only is backend specific.
    python_int_only: dict[str, tuple[bool]] | tuple[bool] | None = None
    # Some functions which seem to be scalar also accept 0d arrays.
    scalar_or_0d_only: dict[str, tuple[bool]] | tuple[bool] | None = None
    # Some functions may not work well with very large integer valued arguments.
    test_large_ints: bool = True
    # Some non-ufunc special functions don't decay 0d arrays to scalar.
    produces_0d: bool = False
    # Whether or not uses native PyTorch or falls back to NumPy/SciPy. This
    # is needed because in PyTorch, the default dtype affects promotion
    # rules when mixing integer and floating dtypes, so relying on a
    # NumPy/SciPy fallback when the default dtype is other than float64 can lead
    # to float64 output when native PyTorch would have e.g. float32 output. This
    # must be accounted for in tests. Not putting this in xp_capabilities for now,
    # but in the future I think it's likely we may want to add a warning to
    # xp_capabilities when not using native PyTorch on CPU.
    torch_native: bool = True

    @property
    def name(self):
        return self.func.__name__

    # These are needed by @lru_cache below
    def __hash__(self):
        return hash(self.func)

    def __eq__(self, other):
        return isinstance(other, _FuncInfo) and self.func == other.func

    @property
    def wrapper(self):
        if self.name in globals():
            # Already initialised. We are likely in a unit test.
            # Return function potentially overridden by xpx.testing.lazy_xp_function.
            import scipy.special
            return getattr(scipy.special, self.name)

        if SCIPY_ARRAY_API:
            @functools.wraps(self.func)
            def wrapped(*args, **kwargs):
                xp = array_namespace(*args)
                return self._wrapper_for(xp)(*args, **kwargs)

            # Allow pickling the function. Normally this is done by @wraps,
            # but in this case it doesn't work because self.func is a ufunc.
            wrapped.__module__ = "scipy.special"
            wrapped.__qualname__ = self.name
            func = wrapped
        else:
            func = self.func

        capabilities = self.xp_capabilities or xp_capabilities()
        # In order to retain a naked ufunc when SCIPY_ARRAY_API is
        # disabled, xp_capabilities must apply its changes in place.
        cap_func = capabilities(func)
        assert cap_func is func
        return func

    @functools.lru_cache(1000)
    def _wrapper_for(self, xp):
        if is_numpy(xp):
            return self.func

        # If a native implementation is available, use that
        spx = scipy_namespace_for(xp)
        f = _get_native_func(xp, spx, self.name, alt_names_map=self.alt_names_map)
        if f is not None:
            return f

        # If generic Array API implementation is available, use that
        if self.generic_impl is not None:
            f = self.generic_impl(xp, spx)
            if f is not None:
                return f

        if is_marray(xp):
            # Unwrap the array, apply the function on the wrapped namespace,
            # and then re-wrap it.
            # IMPORTANT: this only works because all functions in this module
            # are elementwise. Otherwise, we would not be able to define a
            # general rule for mask propagation.

            _f = globals()[self.name]  # Allow nested wrapping
            def f(*args, _f=_f, xp=xp, **kwargs):
                data_args = [getattr(arg, 'data', arg) for arg in args]
                out = _f(*data_args, **kwargs)
                mask = functools.reduce(operator.or_,
                                        (getattr(arg, 'mask', False) for arg in args))
                return xp.asarray(out, mask=mask)

            return f

        if is_dask(xp):
            # Apply the function to each block of the Dask array.
            # IMPORTANT: map_blocks works only because all functions in this module
            # are elementwise. It would be a grave mistake to apply this to gufuncs
            # or any other function with reductions, as they would change their
            # output depending on chunking!

            _f = globals()[self.name]  # Allow nested wrapping
            def f(*args, _f=_f, xp=xp, **kwargs):
                # Hide dtype kwarg from map_blocks
                return xp.map_blocks(functools.partial(_f, **kwargs), *args)

            return f

        # As a final resort, use the NumPy/SciPy implementation
        _f = self.func
        def f(*args, _f=_f, xp=xp, **kwargs):
            # TODO use xpx.lazy_apply to add jax.jit support
            # (but dtype propagation can be non-trivial)
            args = [np.asarray(arg) for arg in args]
            out = _f(*args, **kwargs)
            return xp.asarray(out)

        return f


def _get_native_func(xp, spx, f_name, *, alt_names_map=None):
    if alt_names_map is None:
        alt_names_map = {}
    f_name = alt_names_map.get(get_native_namespace_name(xp), f_name)
    f = getattr(spx.special, f_name, None) if spx else None
    if f is None and hasattr(xp, 'special'):
        # Currently dead branch, in anticipation of 'special' Array API extension
        # https://github.com/data-apis/array-api/issues/725
        f = getattr(xp.special, f_name, None)
    return f


def _rel_entr(xp, spx):
    def __rel_entr(x, y, *, xp=xp):
        # https://github.com/data-apis/array-api-extra/issues/160
        mxp = array_namespace(x._meta, y._meta) if is_dask(xp) else xp
        x, y = xp_promote(x, y, broadcast=True, force_floating=True, xp=xp)
        xy_pos = (x > 0) & (y > 0)
        xy_inf = xp.isinf(x) & xp.isinf(y)
        res = xpx.apply_where(
            xy_pos & ~xy_inf,
            (x, y),
            # Note: for very large x, this can overflow.
            lambda x, y: x * (mxp.log(x) - mxp.log(y)),
            fill_value=xp.inf
        )
        res = xpx.at(res)[(x == 0) & (y >= 0)].set(0)
        res = xpx.at(res)[xp.isnan(x) | xp.isnan(y) | (xy_pos & xy_inf)].set(xp.nan)
        return res

    return __rel_entr


def _xlogy(xp, spx):
    def __xlogy(x, y, *, xp=xp):
        x, y = xp_promote(x, y, force_floating=True, xp=xp)
        with np.errstate(divide='ignore', invalid='ignore'):
            temp = x * xp.log(y)
        return xp.where(x == 0., 0., temp)
    return __xlogy



def _chdtr(xp, spx):
    # The difference between this and just using `gammainc`
    # defined by `get_array_special_func` is that if `gammainc`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `chdtr`.
    gammainc = _get_native_func(xp, spx, 'gammainc')
    if gammainc is None:
        return None

    def __chdtr(v, x):
        res = gammainc(v / 2, x / 2)  # this is almost all we need
        # The rest can be removed when google/jax#20507 is resolved
        mask = (v == 0) & (x > 0)  # JAX returns NaN
        res = xp.where(mask, 1., res)
        mask = xp.isinf(v) & xp.isinf(x)  # JAX returns 1.0
        return xp.where(mask, xp.nan, res)
    return __chdtr


def _chdtrc(xp, spx):
    # The difference between this and just using `gammaincc`
    # defined by `get_array_special_func` is that if `gammaincc`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `chdtrc`.
    gammaincc = _get_native_func(xp, spx, 'gammaincc')
    if gammaincc is None:
        return None

    def __chdtrc(v, x):
        res = xp.where(x >= 0, gammaincc(v/2, x/2), 1)
        i_nan = ((x == 0) & (v == 0)) | xp.isnan(x) | xp.isnan(v) | (v <= 0)
        res = xp.where(i_nan, xp.nan, res)
        return res
    return __chdtrc


def _betaincc(xp, spx):
    betainc = _get_native_func(xp, spx, 'betainc')
    if betainc is None:
        return None

    def __betaincc(a, b, x):
        # not perfect; might want to just rely on SciPy
        return betainc(b, a, 1-x)
    return __betaincc


def _stdtr(xp, spx):
    betainc = _get_native_func(xp, spx, 'betainc')
    if betainc is None:
        return None

    def __stdtr(df, t):
        x = df / (t ** 2 + df)
        tail = betainc(df / 2, 0.5, x) / 2
        return xp.where(t < 0, tail, 1 - tail)

    return __stdtr


def _stdtrit(xp, spx):
    # Need either native stdtr or native betainc
    stdtr = _get_native_func(xp, spx, 'stdtr') or _stdtr(xp, spx)
    # If betainc is not defined, the root-finding would be done with `xp`
    # despite `stdtr` being evaluated with SciPy/NumPy `stdtr`. Save the
    # conversions: in this case, just evaluate `stdtrit` with SciPy/NumPy.
    if stdtr is None:
        return None

    from scipy.optimize.elementwise import bracket_root, find_root

    def __stdtrit(df, p):
        def fun(t, df, p):  return stdtr(df, t) - p
        res_bracket = bracket_root(fun, xp.zeros_like(p), args=(df, p))
        res_root = find_root(fun, res_bracket.bracket, args=(df, p))
        return res_root.x

    return __stdtrit


# Inventory of automatically dispatched functions
# IMPORTANT: these must all be **elementwise** functions!

# PyTorch doesn't implement `betainc`.
# On torch CPU we can fall back to NumPy, but on GPU it won't work.
_needs_betainc = xp_capabilities(cpu_only=True, exceptions=["jax.numpy", "cupy"])

_special_funcs = (
    _FuncInfo(
        _ufuncs.bdtr, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(False, True, False), torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.bdtrc, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(False, True, False), torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.bdtri, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(False, True, False), torch_native=False,
    ),
    _FuncInfo(_ufuncs.betainc, 3, _needs_betainc, torch_native=False),
    _FuncInfo(_ufuncs.betaincc, 3, _needs_betainc, generic_impl=_betaincc,
              torch_native=False),
    _FuncInfo(
        _ufuncs.betaincinv, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        test_large_ints=False, positive_only=True, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.betaln, 2,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        # For betaln, nan mismatches can occur at negative integer a or b of
        # sufficiently large magnitude.
        positive_only={"jax.numpy": True}, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.binom, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.boxcox, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.boxcox1p, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cbrt, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.chdtr, 2, generic_impl=_chdtr),
    _FuncInfo(_ufuncs.chdtrc, 2, generic_impl=_chdtrc,
              # scipy/scipy#20972
              positive_only={"cupy": True, "jax.numpy": True, "torch": True}),
    _FuncInfo(
        _ufuncs.chdtri, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cosdg, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        test_large_ints=False, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cosm1, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cotdg, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipk, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipkm1, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.entr, 1),
    _FuncInfo(_ufuncs.erf, 1),
    _FuncInfo(_ufuncs.erfc, 1),
    _FuncInfo(
        _ufuncs.erfcx, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.erfinv, 1),
    _FuncInfo(
        _ufuncs.exp1, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.exp10, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.exp2, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.exprel, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.expi, 1,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.expit, 1),
    _FuncInfo(
        _ufuncs.expn, 2,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        # Inconsistent behavior for negative n. expn is not defined here without
        # taking analytic continuation.
        positive_only=True,
        int_only=(True, False), test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fdtr, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fdtrc, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fdtri, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gamma, 1,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.gammainc, 2),
    _FuncInfo(
        _ufuncs.gammaincc, 2,
        # google/jax#20699
        positive_only={"jax.numpy": True},
    ),
    _FuncInfo(
        _ufuncs.gammainccinv, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gammaincinv, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.gammaln, 1),
    _FuncInfo(
        _ufuncs.gammasgn, 1,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gdtr, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gdtrc, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.huber, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.hyp1f1, 3,
        xp_capabilities(cpu_only=True, exceptions=["jax.numpy"]),
        positive_only={"jax.numpy": True}, test_large_ints=False,
        torch_native=False,
    ),
    # Comment out when jax>=0.6.1 is available in Conda for CI.
    # (or add version requirements to xp_capabilities).
    # _FuncInfo(
    #     _ufuncs.hyp2f1, 4,
    #     xp_capabilities(cpu_only=True, exceptions=["jax.numpy"]),
    #     positive_only={"jax.numpy": True}, test_large_ints=False,
    #     torch_native=False,
    # ),
    _FuncInfo(
        _ufuncs.inv_boxcox, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.inv_boxcox1p, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.i0, 1),
    _FuncInfo(_ufuncs.i0e, 1),
    _FuncInfo(_ufuncs.i1, 1),
    _FuncInfo(_ufuncs.i1e, 1),
    _FuncInfo(
        _ufuncs.j0, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "bessel_j0"}, test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.j1, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "bessel_j1"}, test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.k0, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "modified_bessel_k0"},
    ),
    _FuncInfo(
        _ufuncs.k0e, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "scaled_modified_bessel_k0"},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.k1, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "modified_bessel_k1"},
    ),
    _FuncInfo(
        _ufuncs.k1e, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "scaled_modified_bessel_k1"},
        test_large_ints=False),
    _FuncInfo(
        _ufuncs.kl_div, 2,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.log_ndtr, 1),
    _FuncInfo(
        _ufuncs.loggamma, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.logit, 1),
    _FuncInfo(
        _ufuncs.lpmv, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
        test_large_ints=False,
    ),
    _FuncInfo(
        _spfun_stats.multigammaln, 2,
        is_ufunc=False,
        python_int_only={
            "cupy": [False, True],
            "jax.numpy": [False, True],
            "torch": [False, True],
        },
        scalar_or_0d_only={
            "array_api_strict": [False, True],
            "numpy": [False, True],
            "dask.array": [False, True],
            "marray": [False, True],
        },
        int_only=(False, True), test_large_ints=False,
        positive_only=True, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtr, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(True, True, False), positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtrc, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(True, True, False), positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtri, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(True, True, False), positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.ndtr, 1),
    _FuncInfo(_ufuncs.ndtri, 1),
    _FuncInfo(
        _ufuncs.pdtr, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        positive_only=True, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pdtrc, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        positive_only=True, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pdtri, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        int_only=(True, False), positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.poch, 2,
        xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
        test_large_ints=False, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pseudo_huber, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _basic.polygamma, 2, int_only=(True, False), is_ufunc=False,
              scalar_or_0d_only={"torch": (True, False)}, produces_0d=True,
              positive_only={"torch": (True, False), "jax.numpy": True},
              test_large_ints=False,
    ),
    _FuncInfo(_ufuncs.psi, 1, alt_names_map={"jax.numpy": "digamma"}),
    _FuncInfo(
        _ufuncs.radian, 3,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.rel_entr, 2, generic_impl=_rel_entr),
    _FuncInfo(
        _ufuncs.rgamma, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _basic.sinc, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        is_ufunc=False,
    ),
    _FuncInfo(
        _ufuncs.sindg, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        test_large_ints=False, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.spence, 1,
        xp_capabilities(cpu_only=True, exceptions=["jax.numpy"]),
        torch_native=False,
    ),
    _FuncInfo(_ufuncs.stdtr,  2, _needs_betainc, generic_impl=_stdtr,
              torch_native=False),
    _FuncInfo(
        _ufuncs.stdtrit, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],  # needs betainc
            skip_backends=[("jax.numpy", "no scipy.optimize support")],
        ),
        generic_impl=_stdtrit, torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.tandg, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        test_large_ints=False, torch_native=False,
    ),
    _FuncInfo(_ufuncs.xlog1py, 2),
    _FuncInfo(_ufuncs.xlogy, 2, generic_impl=_xlogy),
    _FuncInfo(
        _ufuncs.y0, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "bessel_y0"}, test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.y1, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy", "torch"],
            jax_jit=False,
        ),
        alt_names_map={"torch": "bessel_y1"}, test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.yn, 2,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        positive_only={"cupy": (True, False)}, int_only=(True, False),
        test_large_ints=False, torch_native=False,
    ),
    _FuncInfo(
        _basic.zeta, 2, is_ufunc=False,
        positive_only={"jax.numpy": True, "torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.zetac, 1,
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            jax_jit=False,
        ),
        torch_native=False,
    ),
)

# Override ufuncs.
# When SCIPY_ARRAY_API is disabled, this exclusively updates the docstrings in place
# and populates the xp_capabilities table, while retaining the original ufuncs.
globals().update({nfo.func.__name__: nfo.wrapper for nfo in _special_funcs})
# digamma is an alias for psi. Define here so it also has alternative backend
# support. Add noqa because the linter gets confused by the sneaky way psi
# is inserted into globals above.
digamma = psi  # noqa: F821
__all__ = [nfo.func.__name__ for nfo in _special_funcs] + ["digamma"]
