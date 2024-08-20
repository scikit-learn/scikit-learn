# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from dataclasses import dataclass, field

import numpy as np


@dataclass
class InputTags:
    """Tags for the input data.

    Parameters
    ----------
    one_d_array : bool
        Whether the input can be a 1D array.

    two_d_array : bool
        Whether the input can be a 2D array.

    three_d_array : bool
        Whether the input can be a 3D array.

    one_d_labels : bool
        Whether the input is a 1D labels(y).

    two_d_labels : bool
        Whether the input is a 2D labels(y).

    sparse : bool
        Whether the input can be a sparse matrix.

    categorical : bool
        Whether the input can be categorical.

    string : bool
        Whether the input can be an array-like of strings.

    dict : bool
        Whether the input can be a dictionary.

    positive_only : bool
        Whether the input has to be positive.

    allow_nan : bool
        Whether the input can contain NaNs.

    pairwise : bool
        Whether the input is in the form of a calculated pairwise distances or computed
        kernel values.
    """

    one_d_array = False
    two_d_array: bool = True
    three_d_array: bool = False
    one_d_labels = False
    two_d_labels: bool = False
    sparse: bool = False
    categorical: bool = False
    string: bool = False
    dict: bool = False
    positive_only: bool = False
    allow_nan: bool = False
    pairwise: bool = False


@dataclass
class TargetTags:
    """Tags for the target data.

    Parameters
    ----------
    required : bool
        Whether the target is required.

    positive_only : bool
        Whether the target has to be positive.

    multi_output : bool
        Whether the target can be multi-output.

    single_output : bool
        Whether the target can be single-output.
    """

    required: bool
    positive_only: bool = False
    multi_output: bool = False
    single_output: bool = True


@dataclass
class TransformerTags:
    """Tags for the transformer.

    Parameters
    ----------
    preserves_dtype : list[object]
        The data types that the transformer preserves.
    """

    preserves_dtype: list[object] = field(default_factory=lambda: [np.float64])


@dataclass
class ClassifierTags:
    """Tags for the classifier.

    Parameters
    ----------
    poor_score : bool
        Whether the classifier can have a poor score in tests.

    binary : bool
        Whether the classifier can handle binary classification.

    multi_class : bool
        Whether the classifier can handle multi-class classification.

    multi_label : bool
        Whether the classifier can handle multi-label classification.
    """

    poor_score: bool = False
    binary: bool = True
    multi_class: bool = True
    multi_label: bool = False


@dataclass
class RegressorTags:
    """Tags for the regressor.

    Parameters
    ----------
    poor_score : bool
        Whether the regressor can have a poor score in tests.

    multi_label : bool
        Whether the regressor can handle multi-label regression.
    """

    poor_score: bool = False
    multi_label: bool = False


@dataclass
class Tags:
    """Tags for the estimator.

    Parameters
    ----------
    target_tags : TargetTags
        The target(y) tags.

    transformer_tags : TransformerTags
        The transformer tags.

    classifier_tags : ClassifierTags
        The classifier tags.

    regressor_tags : RegressorTags
        The regressor tags.

    array_api_support : bool
        Whether the estimator supports array API supporting input.

    no_validation : bool
        Whether the estimator does not validate input.

    non_deterministic : bool
        Whether the estimator is non-deterministic.

    requires_fit : bool
        Whether the estimator requires fitting before other methods can be called.

    _skip_test : bool
        Whether the estimator should be skipped in tests.

    _xfail_checks : dict[str, str]
        Checks that should be xfailed.

    input_tags : InputTags
        The input data(X) tags.
    """

    target_tags: TargetTags
    transformer_tags: TransformerTags
    classifier_tags: ClassifierTags
    regressor_tags: RegressorTags
    array_api_support: bool = False
    no_validation: bool = False
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
        transformer_tags=(
            TransformerTags()
            if hasattr(estimator, "transform") or hasattr(estimator, "fit_transform")
            else None
        ),
        classifier_tags=ClassifierTags() if is_classifier(estimator) else None,
        regressor_tags=RegressorTags() if is_regressor(estimator) else None,
    )


def _safe_tags(estimator) -> Tags:
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
