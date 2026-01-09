"""Various utilities to help with development."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.exceptions import DataConversionWarning
from sklearn.utils import metadata_routing
from sklearn.utils._bunch import Bunch
from sklearn.utils._chunking import gen_batches, gen_even_slices

# Make _safe_indexing importable from here for backward compat as this particular
# helper is considered semi-private and typically very useful for third-party
# libraries that want to comply with scikit-learn's estimator API. In particular,
# _safe_indexing was included in our public API documentation despite the leading
# `_` in its name.
from sklearn.utils._indexing import _safe_indexing, resample, shuffle
from sklearn.utils._mask import safe_mask
from sklearn.utils._repr_html.base import _HTMLDocumentationLinkMixin  # noqa: F401
from sklearn.utils._repr_html.estimator import estimator_html_repr
from sklearn.utils._tags import (
    ClassifierTags,
    InputTags,
    RegressorTags,
    Tags,
    TargetTags,
    TransformerTags,
    get_tags,
)
from sklearn.utils.class_weight import compute_class_weight, compute_sample_weight
from sklearn.utils.deprecation import deprecated
from sklearn.utils.discovery import all_estimators
from sklearn.utils.extmath import safe_sqr
from sklearn.utils.murmurhash import murmurhash3_32
from sklearn.utils.validation import (
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

__all__ = [
    "Bunch",
    "ClassifierTags",
    "DataConversionWarning",
    "InputTags",
    "RegressorTags",
    "Tags",
    "TargetTags",
    "TransformerTags",
    "_safe_indexing",
    "all_estimators",
    "as_float_array",
    "assert_all_finite",
    "check_X_y",
    "check_array",
    "check_consistent_length",
    "check_random_state",
    "check_scalar",
    "check_symmetric",
    "column_or_1d",
    "compute_class_weight",
    "compute_sample_weight",
    "deprecated",
    "estimator_html_repr",
    "gen_batches",
    "gen_even_slices",
    "get_tags",
    "indexable",
    "metadata_routing",
    "murmurhash3_32",
    "resample",
    "safe_mask",
    "safe_sqr",
    "shuffle",
]
