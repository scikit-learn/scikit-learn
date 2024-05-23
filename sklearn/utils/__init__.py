"""Various utilities to help with development."""

import warnings
from collections.abc import Sequence

import numpy as np

from ..exceptions import DataConversionWarning
from . import _joblib, metadata_routing
from ._bunch import Bunch
from ._chunking import gen_batches, gen_even_slices
from ._estimator_html_repr import estimator_html_repr

# Make _safe_indexing importable from here for backward compat as this particular
# helper is considered semi-private and typically very useful for third-party
# libraries that want to comply with scikit-learn's estimator API. In particular,
# _safe_indexing was included in our public API documentation despite the leading
# `_` in its name.
from ._indexing import (
    _safe_indexing,  # noqa
    resample,
    shuffle,
)
from ._mask import safe_mask
from .class_weight import compute_class_weight, compute_sample_weight
from .deprecation import deprecated
from .discovery import all_estimators
from .extmath import safe_sqr
from .murmurhash import murmurhash3_32
from .validation import (
    as_float_array,
    assert_all_finite,
    check_array,
    check_consistent_length,
    check_random_state,
    check_scalar,
    check_symmetric,
    check_X_y,
    column_or_1d,
    indexable,
)

# TODO(1.7): remove parallel_backend and register_parallel_backend
msg = "deprecated in 1.5 to be removed in 1.7. Use joblib.{} instead."
register_parallel_backend = deprecated(msg)(_joblib.register_parallel_backend)


# if a class, deprecated will change the object in _joblib module so we need to subclass
@deprecated(msg)
class parallel_backend(_joblib.parallel_backend):
    pass


__all__ = [
    "murmurhash3_32",
    "as_float_array",
    "assert_all_finite",
    "check_array",
    "check_random_state",
    "compute_class_weight",
    "compute_sample_weight",
    "column_or_1d",
    "check_consistent_length",
    "check_X_y",
    "check_scalar",
    "indexable",
    "check_symmetric",
    "deprecated",
    "parallel_backend",
    "register_parallel_backend",
    "resample",
    "shuffle",
    "all_estimators",
    "DataConversionWarning",
    "estimator_html_repr",
    "Bunch",
    "metadata_routing",
    "safe_sqr",
    "safe_mask",
    "gen_batches",
    "gen_even_slices",
]


# TODO(1.7): remove
def __getattr__(name):
    if name == "IS_PYPY":
        warnings.warn(
            "IS_PYPY is deprecated and will be removed in 1.7.",
            FutureWarning,
        )
        from .fixes import _IS_PYPY

        return _IS_PYPY
    raise AttributeError(f"module {__name__} has no attribute {name}")


# TODO(1.7): remove tosequence
@deprecated("tosequence was deprecated in 1.5 and will be removed in 1.7")
def tosequence(x):
    """Cast iterable x to a Sequence, avoiding a copy if possible.

    Parameters
    ----------
    x : iterable
        The iterable to be converted.

    Returns
    -------
    x : Sequence
        If `x` is a NumPy array, it returns it as a `ndarray`. If `x`
        is a `Sequence`, `x` is returned as-is. If `x` is from any other
        type, `x` is returned casted as a list.
    """
    if isinstance(x, np.ndarray):
        return np.asarray(x)
    elif isinstance(x, Sequence):
        return x
    else:
        return list(x)
