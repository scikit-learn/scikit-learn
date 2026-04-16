# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.odr` namespace for importing the functions
# included below.


__all__ = [  # noqa: F822
    'odr', 'OdrWarning', 'OdrError', 'OdrStop',
    'Data', 'RealData', 'Model', 'Output', 'ODR',
    'odr_error', 'odr_stop'
]


def __dir__():
    return __all__


def __getattr__(name):
    msg = ("`scipy.odr` is deprecated as of version 1.17.0 and will be removed in "
           "SciPy 1.19.0. Please use `https://pypi.org/project/odrpack/` instead.")
    if name not in __all__:
        raise AttributeError(
            f"`scipy.odr.odrpack` has no attribute {name}. In addition, {msg}")

    import warnings
    from . import _odrpack
    warnings.warn(msg, category=DeprecationWarning, stacklevel=2)

    return getattr(_odrpack, name)
