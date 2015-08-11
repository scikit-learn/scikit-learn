"""
The :mod:`sklearn.utils` module includes various utilities.
"""
from .validation import as_float_array
from .validation import assert_all_finite
from .validation import check_array
from .validation import check_consistent_length
from .validation import check_random_state
from .validation import check_symmetric
from .validation import check_X_y
from .validation import column_or_1d
from .validation import indexable
from .validation import DataConversionWarning

from .base import deprecated
from .base import gen_batches
from .base import gen_even_slices
from .base import resample
from .base import safe_indexing
from .base import safe_mask
from .base import safe_sqr
from .base import shuffle
from .base import tosequence
from .base import ConvergenceWarning
from .base import DataDimensionalityWarning

from .class_weight import compute_class_weight, compute_sample_weight

from .murmurhash import murmurhash3_32


__all__ = ["as_float_array", "assert_all_finite", "check_array",
           "check_consistent_length", "check_random_state",
           "check_symmetric", "check_X_y", "column_or_1d",
           "compute_class_weight", "compute_sample_weight",
           "deprecated", "gen_batches", "gen_even_slices", "indexable",
           "murmurhash3_32", "resample", "safe_indexing", "safe_mask",
           "safe_sqr", "shuffle", "tosequence",
           "ConvergenceWarning", "DataConversionWarning",
           "DataDimensionalityWarning"]
