from __future__ import annotations

from dataclasses import dataclass, field

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


@dataclass
class InputTags:
    """Tags for the input data.

    Parameters
    ----------
    one_d_array : bool, default=False
        Whether the input can be a 1D array.

    two_d_array : bool, default=True
        Whether the input can be a 2D array. Note that most common
        tests currently run only if this flag is set to ``True``.

    three_d_array : bool, default=False
        Whether the input can be a 3D array.

    sparse : bool, default=False
        Whether the input can be a sparse matrix.

    categorical : bool, default=False
        Whether the input can be categorical.

    string : bool, default=False
        Whether the input can be an array-like of strings.

    dict : bool, default=False
        Whether the input can be a dictionary.

    positive_only : bool, default=False
        Whether the estimator requires positive X.

    allow_nan : bool, default=False
        Whether the estimator supports data with missing values encoded as `np.nan`.

    pairwise : bool, default=False
        This boolean attribute indicates whether the data (`X`),
        :term:`fit` and similar methods consists of pairwise measures
        over samples rather than a feature representation for each
        sample.  It is usually `True` where an estimator has a
        `metric` or `affinity` or `kernel` parameter with value
        'precomputed'. Its primary purpose is to support a
        :term:`meta-estimator` or a cross validation procedure that
        extracts a sub-sample of data intended for a pairwise
        estimator, where the data needs to be indexed on both axes.
        Specifically, this tag is used by
        `sklearn.utils.metaestimators._safe_split` to slice rows and
        columns.
    """

    one_d_array: bool = False
    two_d_array: bool = True
    three_d_array: bool = False
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
        Whether the estimator requires y to be passed to `fit`,
        `fit_predict` or `fit_transform` methods. The tag is ``True``
        for estimators inheriting from `~sklearn.base.RegressorMixin`
        and `~sklearn.base.ClassifierMixin`.

    one_d_labels : bool, default=False
        Whether the input is a 1D labels (y).

    two_d_labels : bool, default=False
        Whether the input is a 2D labels (y).

    positive_only : bool, default=False
        Whether the estimator requires a positive y (only applicable
        for regression).

    multi_output : bool, default=False
        Whether a regressor supports multi-target outputs or a classifier supports
        multi-class multi-output.

    single_output : bool, default=True
        Whether the target can be single-output. This can be ``False`` if the
        estimator supports only multi-output cases.
    """

    required: bool
    one_d_labels: bool = False
    two_d_labels: bool = False
    positive_only: bool = False
    multi_output: bool = False
    single_output: bool = True


@dataclass
class TransformerTags:
    """Tags for the transformer.

    Parameters
    ----------
    preserves_dtype : list[str], default=["float64"]
        Applies only on transformers. It corresponds to the data types
        which will be preserved such that `X_trans.dtype` is the same
        as `X.dtype` after calling `transformer.transform(X)`. If this
        list is empty, then the transformer is not expected to
        preserve the data type. The first value in the list is
        considered as the default data type, corresponding to the data
        type of the output when the input data type is not going to be
        preserved.
    """

    preserves_dtype: list[str] = field(default_factory=lambda: ["float64"])


@dataclass
class ClassifierTags:
    """Tags for the classifier.

    Parameters
    ----------
    poor_score : bool, default=False
        Whether the estimator fails to provide a "reasonable" test-set
        score, which currently for classification is an accuracy of
        0.83 on ``make_blobs(n_samples=300, random_state=0)``. The
        datasets and values are based on current estimators in scikit-learn
        and might be replaced by something more systematic.

    multi_class : bool, default=True
        Whether the classifier can handle multi-class
        classification. Note that all classifiers support binary
        classification. Therefore this flag indicates whether the
        classifier is a binary-classifier-only or not.

    multi_label : bool, default=False
        Whether the classifier supports multi-label output.
    """

    poor_score: bool = False
    multi_class: bool = True
    multi_label: bool = False


@dataclass
class RegressorTags:
    """Tags for the regressor.

    Parameters
    ----------
    poor_score : bool, default=False
        Whether the estimator fails to provide a "reasonable" test-set
        score, which currently for regression is an R2 of 0.5 on
        ``make_regression(n_samples=200, n_features=10,
        n_informative=1, bias=5.0, noise=20, random_state=42)``. The
        dataset and values are based on current estimators in scikit-learn
        and might be replaced by something more systematic.

    multi_label : bool, default=False
        Whether the regressor supports multilabel output.
    """

    poor_score: bool = False
    multi_label: bool = False


@dataclass
class Tags:
    """Tags for the estimator.

    Parameters
    ----------
    target_tags : :class:`TargetTags`
        The target(y) tags.

    transformer_tags : :class:`TransformerTags` or None
        The transformer tags.

    classifier_tags : :class:`ClassifierTags` or None
        The classifier tags.

    regressor_tags : :class:`RegressorTags` or None
        The regressor tags.

    array_api_support : bool, default=False
        Whether the estimator supports Array API compatible inputs.

    no_validation : bool, default=False
        Whether the estimator skips input-validation. This is only meant for
        stateless and dummy transformers!

    non_deterministic : bool, default=False
        Whether the estimator is not deterministic given a fixed ``random_state``.

    requires_fit : bool, default=True
        Whether the estimator requires to be fitted before calling one of
        `transform`, `predict`, `predict_proba`, or `decision_function`.

    _skip_test : bool, default=False
        Whether to skip common tests entirely. Don't use this unless
        you have a *very good* reason.

    _xfail_checks : dict[str, str], default={}
        Dictionary ``{check_name: reason}`` of common checks that will
        be marked as `XFAIL` for pytest, when using
        :func:`~sklearn.utils.estimator_checks.parametrize_with_checks`. These
        checks will be simply ignored and not run by
        :func:`~sklearn.utils.estimator_checks.check_estimator`, but a
        `SkipTestWarning` will be raised.  Don't use this unless there
        is a *very good* reason for your estimator not to pass the
        check.  Also note that the usage of this tag is highly subject
        to change because we are trying to make it more flexible: be
        prepared for breaking changes in the future.

    input_tags : :class:`InputTags`
        The input data(X) tags.
    """

    target_tags: TargetTags
    transformer_tags: TransformerTags | None
    classifier_tags: ClassifierTags | None
    regressor_tags: RegressorTags | None
    array_api_support: bool = False
    no_validation: bool = False
    non_deterministic: bool = False
    requires_fit: bool = True
    _skip_test: bool = False
    _xfail_checks: dict[str, str] = field(default_factory=dict)
    input_tags: InputTags = field(default_factory=InputTags)


def default_tags(estimator) -> Tags:
    """Get the default tags for an estimator.

    This ignores any ``__sklearn_tags__`` method that the estimator may have.

    If the estimator is a classifier or a regressor, ``target_tags.required``
    will be set to ``True``, otherwise it will be set to ``False``.

    ``transformer_tags`` will be set to :class:`~.sklearn.utils. TransformerTags` if the
    estimator has a ``transform`` or ``fit_transform`` method, otherwise it will be set
    to ``None``.

    ``classifier_tags`` will be set to :class:`~.sklearn.utils.ClassifierTags` if the
    estimator is a classifier, otherwise it will be set to ``None``.
    a classifier, otherwise it will be set to ``None``.

    ``regressor_tags`` will be set to :class:`~.sklearn.utils.RegressorTags` if the
    estimator is a regressor, otherwise it will be set to ``None``.

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


def get_tags(estimator) -> Tags:
    """Get estimator tags.

    :class:`~sklearn.BaseEstimator` provides the estimator tags machinery.
    However, if an estimator does not inherit from this base class, we should
    fall-back to the default tags.

    For scikit-learn built-in estimators, we should still rely on
    `self.__sklearn_tags__()`. `get_tags(est)` should be used when we
    are not sure where `est` comes from: typically
    `get_tags(self.estimator)` where `self` is a meta-estimator, or in
    the common checks.

    Parameters
    ----------
    estimator : estimator object
        The estimator from which to get the tag.

    Returns
    -------
    tags : :class:`~.sklearn.utils.Tags`
        The estimator tags.
    """
    if hasattr(estimator, "__sklearn_tags__"):
        tags = estimator.__sklearn_tags__()
    else:
        tags = default_tags(estimator)

    return tags
