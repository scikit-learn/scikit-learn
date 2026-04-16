import functools
import types
from scipy._lib._array_api import (
    is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API, xp_capabilities
)

from ._signal_api import *   # noqa: F403
from . import _signal_api
from . import _delegators
__all__ = _signal_api.__all__


MODULE_NAME = 'signal'

# jax.scipy.signal has only partial coverage of scipy.signal, so we keep the list
# of functions we can delegate to JAX
# https://jax.readthedocs.io/en/latest/jax.scipy.html
JAX_SIGNAL_FUNCS = [
    'fftconvolve', 'convolve', 'convolve2d', 'correlate', 'correlate2d',
    'csd', 'detrend', 'istft', 'welch'
]

# some cupyx.scipy.signal functions are incompatible with their scipy counterparts
CUPY_BLACKLIST = [
    'lfilter_zi', 'sosfilt_zi', 'get_window', 'besselap', 'envelope', 'remez', 'bessel'
]

# freqz_sos is a sosfreqz rename, and cupy does not have the new name yet (in v13.x)
CUPY_RENAMES = {'freqz_sos': 'sosfreqz'}


def delegate_xp(delegator, module_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwds):
            try:
                xp = delegator(*args, **kwds)
            except TypeError:
                # object arrays
                if func.__name__ == "tf2ss":
                    import numpy as np
                    xp = np
                else:
                    raise

            # try delegating to a cupyx/jax namesake
            if is_cupy(xp) and func.__name__ not in CUPY_BLACKLIST:
                func_name = CUPY_RENAMES.get(func.__name__, func.__name__)

                # https://github.com/cupy/cupy/issues/8336
                import importlib
                cupyx_module = importlib.import_module(f"cupyx.scipy.{module_name}")
                cupyx_func = getattr(cupyx_module, func_name)
                kwds.pop('xp', None)
                return cupyx_func(*args, **kwds)
            elif is_jax(xp) and func.__name__ in JAX_SIGNAL_FUNCS:
                spx = scipy_namespace_for(xp)
                jax_module = getattr(spx, module_name)
                jax_func = getattr(jax_module, func.__name__)
                kwds.pop('xp', None)
                return jax_func(*args, **kwds)
            else:
                # the original function
                return func(*args, **kwds)
        return wrapper
    return inner


# Although most of these functions currently exist in CuPy and some in JAX,
# there are no alternative backend tests for any of them in the current
# test suite. Each will be documented as np_only until tests are added.
untested = {
    "argrelextrema",
    "argrelmax",
    "argrelmin",
    "band_stop_obj",
    "check_NOLA",
    "chirp",
    "coherence",
    "csd",
    "czt_points",
    "dbode",
    "dfreqresp",
    "dlsim",
    "dstep",
    "find_peaks",
    "find_peaks_cwt",
    "findfreqs",
    "freqresp",
    "gausspulse",
    "lombscargle",
    "lsim",
    "max_len_seq",
    "peak_prominences",
    "peak_widths",
    "periodogram",
    "place_pols",
    "sawtooth",
    "sepfir2d",
    "square",
    "ss2tf",
    "ss2zpk",
    "step",
    "sweep_poly",
    "symiirorder1",
    "symiirorder2",
    "tf2ss",
    "unit_impulse",
    "welch",
    "zoom_fft",
    "zpk2ss",
}


def get_default_capabilities(func_name, delegator):
    if delegator is None or func_name in untested:
        return xp_capabilities(np_only=True)
    return xp_capabilities()

bilinear_extra_note = \
    """CuPy does not accept complex inputs.

    """

uses_choose_conv_extra_note = \
    """CuPy does not support inputs with ``ndim>1`` when ``method="auto"``
    but does support higher dimensional arrays for ``method="direct"``
    and ``method="fft"``.

    """

resample_poly_extra_note = \
    """CuPy only supports ``padtype="constant"``.

    """

upfirdn_extra_note = \
    """CuPy only supports ``mode="constant"`` and ``cval=0.0``.

    """

xord_extra_note = \
    """The ``torch`` backend on GPU does not support the case where
    `wp` and `ws` specify a Bandstop filter.

    """

convolve2d_extra_note = \
    """The JAX backend only supports ``boundary="fill"`` and ``fillvalue=0``.

    """

zpk2tf_extra_note = \
    """The CuPy and JAX backends both support only 1d input.

    """

capabilities_overrides = {
    "bessel": xp_capabilities(cpu_only=True, jax_jit=False, allow_dask_compute=True),
    "bilinear": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                jax_jit=False, allow_dask_compute=True,
                                reason="Uses np.polynomial.Polynomial",
                                extra_note=bilinear_extra_note),
    "bilinear_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                    jax_jit=False, allow_dask_compute=True),
    "butter": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False, 
                              allow_dask_compute=True),
    "buttord": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                               jax_jit=False, allow_dask_compute=True,
                               extra_note=xord_extra_note),
    "cheb1ord": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                jax_jit=False, allow_dask_compute=True,
                                extra_note=xord_extra_note),
    "cheb2ord": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                jax_jit=False, allow_dask_compute=True,
                                extra_note=xord_extra_note),
    "cheby1": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),

    "cheby2": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "cont2discrete": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "convolve": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                allow_dask_compute=True,
                                extra_note=uses_choose_conv_extra_note),
    "convolve2d": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                  allow_dask_compute=True,
                                  extra_note=convolve2d_extra_note),
    "correlate": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                 allow_dask_compute=True,
                                 extra_note=uses_choose_conv_extra_note),
    "correlate2d": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                   allow_dask_compute=True,
                                   extra_note=convolve2d_extra_note),
    "correlation_lags": xp_capabilities(out_of_scope=True),
    "cspline1d": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                 jax_jit=False, allow_dask_compute=True),
    "cspline1d_eval": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                      jax_jit=False, allow_dask_compute=True),
    "cspline2d": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                 jax_jit=False, allow_dask_compute=True),
    "czt": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "deconvolve": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                  allow_dask_compute=True,
                                  skip_backends=[("jax.numpy", "item assignment")]),
    "decimate": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "detrend": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                               allow_dask_compute=True),
    "dimpulse": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "dlti": xp_capabilities(np_only=True,
                            reason="works in CuPy but delegation isn't set up yet"),
    "ellip": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False, 
                             allow_dask_compute=True,
                             reason="scipy.special.ellipk"),
    "ellipord": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                jax_jit=False, allow_dask_compute=True,
                                reason="scipy.special.ellipk"),
    "firls": xp_capabilities(cpu_only=True, allow_dask_compute=True, jax_jit=False,
                             reason="lstsq"),
    "firwin": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                              jax_jit=False, allow_dask_compute=True),
    "firwin2": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               jax_jit=False, allow_dask_compute=True,
                               reason="firwin uses np.interp"),
    "fftconvolve": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
    "freqs": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             jax_jit=False, allow_dask_compute=True),
    "freqs_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "freqz": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             jax_jit=False, allow_dask_compute=True),
    "freqz_sos": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "group_delay": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                   jax_jit=False, allow_dask_compute=True),
    "hilbert": xp_capabilities(
        cpu_only=True, exceptions=["cupy", "torch"],
        skip_backends=[("jax.numpy", "item assignment")],
    ),
    "hilbert2": xp_capabilities(
        cpu_only=True, exceptions=["cupy", "torch"],
        skip_backends=[("jax.numpy", "item assignment")],
    ),
    "invres": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "invresz": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "iircomb": xp_capabilities(xfail_backends=[("jax.numpy", "inaccurate")]),
    "iirfilter": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "kaiser_atten": xp_capabilities(
        out_of_scope=True, reason="scalars in, scalars out"
    ),
    "kaiser_beta": xp_capabilities(out_of_scope=True, reason="scalars in, scalars out"),
    "kaiserord": xp_capabilities(out_of_scope=True, reason="scalars in, scalars out"),
    "lfilter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True, jax_jit=False),
    "lfilter_zi": xp_capabilities(cpu_only=True, allow_dask_compute=True,
                                  jax_jit=False),
    "lfiltic": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True),
    "lp2bp": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True,
                             skip_backends=[("jax.numpy", "in-place item assignment")]),
    "lp2bp_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lp2bs": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True,
                             skip_backends=[("jax.numpy", "in-place item assignment")]),
    "lp2bs_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lp2lp": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True, jax_jit=False),
    "lp2lp_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lp2hp": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True,
                             skip_backends=[("jax.numpy", "in-place item assignment")]),
    "lp2hp_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lti": xp_capabilities(np_only=True,
                            reason="works in CuPy but delegation isn't set up yet"),
    "medfilt": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True, jax_jit=False,
                               reason="uses scipy.ndimage.rank_filter"),
    "medfilt2d": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                 allow_dask_compute=True, jax_jit=False,
                                 reason="c extension module"),
    "minimum_phase": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                     allow_dask_compute=True, jax_jit=False),
    "normalize": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "oaconvolve": xp_capabilities(
        cpu_only=True, exceptions=["cupy", "torch"],
        skip_backends=[("jax.numpy", "fails all around")],
        xfail_backends=[("dask.array", "wrong answer")],
    ),
    "order_filter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                    allow_dask_compute=True, jax_jit=False,
                                    reason="uses scipy.ndimage.rank_filter"),
    "qspline1d": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                 jax_jit=False, allow_dask_compute=True),
    "qspline1d_eval": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                      jax_jit=False, allow_dask_compute=True),
    "qspline2d": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "remez": xp_capabilities(cpu_only=True, allow_dask_compute=True, jax_jit=False),
    "resample": xp_capabilities(
        cpu_only=True, exceptions=["cupy"],
        skip_backends=[
            ("dask.array", "XXX something in dask"),
            ("jax.numpy", "XXX: immutable arrays"),
        ]
    ),
    "resample_poly": xp_capabilities(
        cpu_only=True, exceptions=["cupy"],
        jax_jit=False, skip_backends=[("dask.array", "XXX something in dask")],
        extra_note=resample_poly_extra_note,
    ),
    "residue": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "residuez": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "savgol_filter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                     jax_jit=False,
                                     reason="convolve1d is cpu-only"),
    "sepfir2d": xp_capabilities(np_only=True),
    "sos2zpk": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                               allow_dask_compute=True),
    "sos2tf": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "sosfilt": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True),
    "sosfiltfilt": xp_capabilities(
        cpu_only=True, exceptions=["cupy"],
        skip_backends=[
            (
                "dask.array",
                "sosfiltfilt directly sets shape attributes on arrays"
                " which dask doesn't like"
            ),
            ("torch", "negative strides"),
            ("jax.numpy", "sosfilt works in-place"),
        ],
    ),
    "sosfreqz": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                jax_jit=False, allow_dask_compute=True),
    "spline_filter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                     jax_jit=False, allow_dask_compute=True),
    "tf2sos": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "tf2zpk": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "unique_roots": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "upfirdn": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                               allow_dask_compute=True,
                               reason="Cython implementation",
                               extra_note=upfirdn_extra_note),
    "vectorstrength": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                      allow_dask_compute=True, jax_jit=False),
    "wiener": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                              allow_dask_compute=True, jax_jit=False,
                              reason="uses scipy.signal.correlate"),
    "zpk2sos": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                               allow_dask_compute=True),
    "zpk2tf": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True,
                              extra_note=zpk2tf_extra_note),
    "spectrogram": xp_capabilities(out_of_scope=True),  # legacy
    "stft": xp_capabilities(out_of_scope=True),  # legacy
    "istft": xp_capabilities(out_of_scope=True),  # legacy
    "check_COLA": xp_capabilities(out_of_scope=True),  # legacy
}


# ### decorate ###
for obj_name in _signal_api.__all__:
    bare_obj = getattr(_signal_api, obj_name)
    delegator = getattr(_delegators, obj_name + "_signature", None)

    if SCIPY_ARRAY_API and delegator is not None:
        f = delegate_xp(delegator, MODULE_NAME)(bare_obj)
    else:
        f = bare_obj

    if not isinstance(f, types.ModuleType):
        capabilities = capabilities_overrides.get(
            obj_name, get_default_capabilities(obj_name, delegator)
        )
        f = capabilities(f)

    # add the decorated function to the namespace, to be imported in __init__.py
    vars()[obj_name] = f
