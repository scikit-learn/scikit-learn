# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from dataclasses import dataclass, field

import numpy as np

_DEFAULT_TAGS = {
    "array_api_support": False,
    "non_deterministic": False,
    "requires_positive_X": False,
    "requires_positive_y": False,
    "X_types": ["2darray"],
    "poor_score": False,
    "no_validation": False,
    "multioutput": False,
    "allow_nan": False,
    "stateless": False,
    "multilabel": False,
    "_skip_test": False,
    "_xfail_checks": False,
    "multioutput_only": False,
    "binary_only": False,
    "requires_fit": True,
    "preserves_dtype": [np.float64],
    "requires_y": False,
    "pairwise": False,
}


@dataclass
class InputTags:
    two_d_array: bool = True
    sparse: bool = False
    categorical: bool = False
    string: bool = False
    positive_only: bool = False
    nan_allowed: bool = False
    pairwise: bool = False


@dataclass
class TargetTags:
    required: bool
    positive_only: bool = False
    multi_output: bool = False
    single_output: bool = True


@dataclass
class TransformerTags:
    preserve_dtype: list[str] = field(default_factory=lambda: ["float64"])


@dataclass
class ClassifierTags:
    poor_score: bool = False
    binary: bool = True
    multiclass: bool = True
    multilabel: bool = False


@dataclass
class RegressorTags:
    poor_score: bool = False


@dataclass
class Tags:
    target_flags: TargetTags
    transformer_flags: TransformerTags
    classifier_flags: ClassifierTags
    regressor_flags: RegressorTags
    array_api_support: bool = False
    no_validation: bool = False
    stateless: bool = False
    non_deterministic: bool = False
    requires_fit: bool = True
    _skip_test: bool = False
    _xfail_checks: dict[str, str] = field(default_factory=dict)
    input_flags: InputTags = field(default_factory=InputTags)


def _safe_tags(estimator, key=None):
    """Safely get estimator tags.

    :class:`~sklearn.BaseEstimator` provides the estimator tags machinery.
    However, if an estimator does not inherit from this base class, we should
    fall-back to the default tags.

    For scikit-learn built-in estimators, we should still rely on
    `self.__sklearn_tags__()`. `_safe_tags(est)` should be used when we are not sure
    where `est` comes from: typically `_safe_tags(self.base_estimator)` where
    `self` is a meta-estimator, or in the common checks.

    Parameters
    ----------
    estimator : estimator object
        The estimator from which to get the tag.

    key : str, default=None
        Tag name to get. By default (`None`), all tags are returned.

    Returns
    -------
    tags : dict or tag value
        The estimator tags. A single value is returned if `key` is not None.
    """
    if hasattr(estimator, "__sklearn_tags__"):
        tags_provider = "__sklearn_tags__()"
        tags = estimator.__sklearn_tags__()
    else:
        tags_provider = "_DEFAULT_TAGS"
        tags = _DEFAULT_TAGS

    if key is not None:
        if key not in tags:
            raise ValueError(
                f"The key {key} is not defined in {tags_provider} for the "
                f"class {estimator.__class__.__name__}."
            )
        return tags[key]
    return tags
