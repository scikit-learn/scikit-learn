# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from dataclasses import dataclass, field

import numpy as np


@dataclass
class InputTags:
    two_d_array: bool = True
    sparse: bool = False
    categorical: bool = False
    string: bool = False
    dict: bool = False
    positive_only: bool = False
    allow_nan: bool = False
    pairwise: bool = False


@dataclass
class TargetTags:
    required: bool
    positive_only: bool = False
    multi_output: bool = False
    single_output: bool = True


@dataclass
class TransformerTags:
    preserves_dtype: list[object] = field(default_factory=lambda: [np.float64])


@dataclass
class ClassifierTags:
    poor_score: bool = False
    binary: bool = True
    multi_class: bool = True
    multi_label: bool = False


@dataclass
class RegressorTags:
    poor_score: bool = False
    multi_label: bool = False


@dataclass
class Tags:
    target_tags: TargetTags
    transformer_tags: TransformerTags
    classifier_tags: ClassifierTags
    regressor_tags: RegressorTags
    array_api_support: bool = False
    no_validation: bool = False
    stateless: bool = False
    non_deterministic: bool = False
    requires_fit: bool = True
    _skip_test: bool = False
    _xfail_checks: dict[str, str] = field(default_factory=dict)
    input_tags: InputTags = field(default_factory=InputTags)


def default_tags(estimator):
    """Get the default tags for an estimator.

    Parameters
    ----------
    estimator : estimator object
        The estimator for which to get the default tags.

    Returns
    -------
    tags : Tags
        The default tags for the estimator.
    """
    from ..base import is_classifier, is_regressor

    target_required = is_classifier(estimator) or is_regressor(estimator)

    return Tags(
        target_tags=TargetTags(required=target_required),
        transformer_tags=TransformerTags() if hasattr(estimator, "transform") else None,
        classifier_tags=ClassifierTags() if is_classifier(estimator) else None,
        regressor_tags=RegressorTags() if is_regressor(estimator) else None,
    )


def _safe_tags(estimator):
    """Safely get estimator tags.

    :class:`~sklearn.BaseEstimator` provides the estimator tags machinery.
    However, if an estimator does not inherit from this base class, we should
    fall-back to the default tags.

    For scikit-learn built-in estimators, we should still rely on
    `self.__sklearn_tags__()`. `_safe_tags(est)` should be used when we are not sure
    where `est` comes from: typically `_safe_tags(self.estimator)` where
    `self` is a meta-estimator, or in the common checks.

    Parameters
    ----------
    estimator : estimator object
        The estimator from which to get the tag.

    Returns
    -------
    tags : Tags
        The estimator tags.
    """
    if hasattr(estimator, "__sklearn_tags__"):
        tags = estimator.__sklearn_tags__()
    else:
        tags = default_tags(estimator)

    return tags
