"""
The :mod:`sklearn.utils` module includes various utilities.
"""
import platform
import struct

from ..exceptions import DataConversionWarning
from . import _joblib, metadata_routing
from ._bunch import Bunch
from ._chunking import gen_batches, gen_even_slices
from ._estimator_html_repr import estimator_html_repr
from ._indexing import _safe_indexing, resample, shuffle
from ._mask import safe_mask
from ._misc import _message_with_time, _print_elapsed_time, _to_object_array, tosequence
from ._optional_dependencies import check_matplotlib_support, check_pandas_support
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

# Do not deprecate parallel_backend and register_parallel_backend as they are
# needed to tune `scikit-learn` behavior and have different effect if called
# from the vendored version or or the site-package version. The other are
# utilities that are independent of scikit-learn so they are not part of
# scikit-learn public API.
parallel_backend = _joblib.parallel_backend
register_parallel_backend = _joblib.register_parallel_backend

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
    "check_matplotlib_support",
    "check_pandas_support",
    "all_estimators",
    "DataConversionWarning",
    "estimator_html_repr",
    "Bunch",
    "metadata_routing",
    "_safe_indexing",
    "gen_batches",
    "gen_even_slices",
    "safe_mask",
    "safe_sqr",
    "tosequence",
    "_to_object_array",
    "_print_elapsed_time",
    "_message_with_time",
]

IS_PYPY = platform.python_implementation() == "PyPy"
_IS_32BIT = 8 * struct.calcsize("P") == 32
