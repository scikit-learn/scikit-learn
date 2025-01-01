from __future__ import annotations

import warnings
from collections import OrderedDict
from dataclasses import dataclass, field
from itertools import chain

from .fixes import _dataclass_args

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


@dataclass(**_dataclass_args())
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

        Note that if setting this tag to ``True`` means the estimator can take only
        positive values, the `positive_only` tag must reflect it and also be set to
        ``True``.
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


@dataclass(**_dataclass_args())
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

        See :term:`multi-output` in the glossary.

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


@dataclass(**_dataclass_args())
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


@dataclass(**_dataclass_args())
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

        See :term:`multi-class` in the glossary.

    multi_label : bool, default=False
        Whether the classifier supports multi-label output: a data point can
        be predicted to belong to a variable number of classes.

        See :term:`multi-label` in the glossary.
    """

    poor_score: bool = False
    multi_class: bool = True
    multi_label: bool = False


@dataclass(**_dataclass_args())
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
    """

    poor_score: bool = False


@dataclass(**_dataclass_args())
class Tags:
    """Tags for the estimator.

    See :ref:`estimator_tags` for more information.

    Parameters
    ----------
    estimator_type : str or None
        The type of the estimator. Can be one of:
        - "classifier"
        - "regressor"
        - "transformer"
        - "clusterer"
        - "outlier_detector"
        - "density_estimator"

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

    input_tags : :class:`InputTags`
        The input data(X) tags.
    """

    estimator_type: str | None
    target_tags: TargetTags
    transformer_tags: TransformerTags | None = None
    classifier_tags: ClassifierTags | None = None
    regressor_tags: RegressorTags | None = None
    array_api_support: bool = False
    no_validation: bool = False
    non_deterministic: bool = False
    requires_fit: bool = True
    _skip_test: bool = False
    input_tags: InputTags = field(default_factory=InputTags)


# TODO(1.8): Remove this function
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
    est_is_classifier = getattr(estimator, "_estimator_type", None) == "classifier"
    est_is_regressor = getattr(estimator, "_estimator_type", None) == "regressor"
    target_required = est_is_classifier or est_is_regressor

    return Tags(
        estimator_type=getattr(estimator, "_estimator_type", None),
        target_tags=TargetTags(required=target_required),
        transformer_tags=(
            TransformerTags()
            if hasattr(estimator, "transform") or hasattr(estimator, "fit_transform")
            else None
        ),
        classifier_tags=ClassifierTags() if est_is_classifier else None,
        regressor_tags=RegressorTags() if est_is_regressor else None,
    )


# TODO(1.7): Remove this function
def _find_tags_provider(estimator, warn=True):
    """Find the tags provider for an estimator.

    Parameters
    ----------
    estimator : estimator object
        The estimator to find the tags provider for.

    warn : bool, default=True
        Whether to warn if the tags provider is not found.

    Returns
    -------
    tag_provider : str
        The tags provider for the estimator. Can be one of:
        - "_get_tags": to use the old tags infrastructure
        - "__sklearn_tags__": to use the new tags infrastructure
    """
    mro_model = type(estimator).mro()
    tags_mro = OrderedDict()
    for klass in mro_model:
        tags_provider = []
        if "_more_tags" in vars(klass):
            tags_provider.append("_more_tags")
        if "_get_tags" in vars(klass):
            tags_provider.append("_get_tags")
        if "__sklearn_tags__" in vars(klass):
            tags_provider.append("__sklearn_tags__")
        tags_mro[klass.__name__] = tags_provider

    all_providers = set(chain.from_iterable(tags_mro.values()))
    if "__sklearn_tags__" not in all_providers:
        # default on the old tags infrastructure
        return "_get_tags"

    tag_provider = "__sklearn_tags__"
    for klass in tags_mro:
        has_get_or_more_tags = any(
            provider in tags_mro[klass] for provider in ("_get_tags", "_more_tags")
        )
        has_sklearn_tags = "__sklearn_tags__" in tags_mro[klass]

        if tags_mro[klass] and tag_provider == "__sklearn_tags__":  # is it empty
            if has_get_or_more_tags and not has_sklearn_tags:
                # Case where a class does not implement __sklearn_tags__ and we fallback
                # to _get_tags. We should therefore warn for implementing
                # __sklearn_tags__.
                tag_provider = "_get_tags"
                break

    if warn and tag_provider == "_get_tags":
        warnings.warn(
            f"The {estimator.__class__.__name__} or classes from which it inherits "
            "use `_get_tags` and `_more_tags`. Please define the "
            "`__sklearn_tags__` method, or inherit from `sklearn.base.BaseEstimator` "
            "and/or other appropriate mixins such as `sklearn.base.TransformerMixin`, "
            "`sklearn.base.ClassifierMixin`, `sklearn.base.RegressorMixin`, and "
            "`sklearn.base.OutlierMixin`. From scikit-learn 1.7, not defining "
            "`__sklearn_tags__` will raise an error.",
            category=FutureWarning,
        )
    return tag_provider


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

    .. versionadded:: 1.6

    Parameters
    ----------
    estimator : estimator object
        The estimator from which to get the tag.

    Returns
    -------
    tags : :class:`~.sklearn.utils.Tags`
        The estimator tags.
    """

    tag_provider = _find_tags_provider(estimator)

    if tag_provider == "__sklearn_tags__":
        tags = estimator.__sklearn_tags__()
    else:
        # TODO(1.7): Remove this branch of the code
        # Let's go through the MRO and patch each class implementing _more_tags
        sklearn_tags_provider = {}
        more_tags_provider = {}
        class_order = []
        for klass in reversed(type(estimator).mro()):
            if "__sklearn_tags__" in vars(klass):
                sklearn_tags_provider[klass] = klass.__sklearn_tags__(estimator)  # type: ignore[attr-defined]
                class_order.append(klass)
            elif "_more_tags" in vars(klass):
                more_tags_provider[klass] = klass._more_tags(estimator)  # type: ignore[attr-defined]
                class_order.append(klass)

        # Find differences between consecutive in the case of __sklearn_tags__
        # inheritance
        sklearn_tags_diff = {}
        items = list(sklearn_tags_provider.items())
        for current_item, next_item in zip(items[:-1], items[1:]):
            current_name, current_tags = current_item
            next_name, next_tags = next_item
            current_tags = _to_old_tags(current_tags)
            next_tags = _to_old_tags(next_tags)

            # Compare tags and store differences
            diff = {}
            for key in current_tags:
                if current_tags[key] != next_tags[key]:
                    diff[key] = next_tags[key]

            sklearn_tags_diff[next_name] = diff

        tags = {}
        for klass in class_order:
            if klass in sklearn_tags_diff:
                tags.update(sklearn_tags_diff[klass])
            elif klass in more_tags_provider:
                tags.update(more_tags_provider[klass])

        tags = _to_new_tags(
            {**_to_old_tags(default_tags(estimator)), **tags}, estimator
        )

    return tags


# TODO(1.7): Remove this function
def _safe_tags(estimator, key=None):
    warnings.warn(
        "The `_safe_tags` function is deprecated in 1.6 and will be removed in "
        "1.7. Use the public `get_tags` function instead and make sure to implement "
        "the `__sklearn_tags__` method.",
        category=FutureWarning,
    )
    tags = _to_old_tags(get_tags(estimator))

    if key is not None:
        if key not in tags:
            raise ValueError(
                f"The key {key} is not defined for the class "
                f"{estimator.__class__.__name__}."
            )
        return tags[key]
    return tags


# TODO(1.7): Remove this function
def _to_new_tags(old_tags, estimator=None):
    """Utility function convert old tags (dictionary) to new tags (dataclass)."""
    input_tags = InputTags(
        one_d_array="1darray" in old_tags["X_types"],
        two_d_array="2darray" in old_tags["X_types"],
        three_d_array="3darray" in old_tags["X_types"],
        sparse="sparse" in old_tags["X_types"],
        categorical="categorical" in old_tags["X_types"],
        string="string" in old_tags["X_types"],
        dict="dict" in old_tags["X_types"],
        positive_only=old_tags["requires_positive_X"],
        allow_nan=old_tags["allow_nan"],
        pairwise=old_tags["pairwise"],
    )
    target_tags = TargetTags(
        required=old_tags["requires_y"],
        one_d_labels="1dlabels" in old_tags["X_types"],
        two_d_labels="2dlabels" in old_tags["X_types"],
        positive_only=old_tags["requires_positive_y"],
        multi_output=old_tags["multioutput"] or old_tags["multioutput_only"],
        single_output=not old_tags["multioutput_only"],
    )
    if estimator is not None and (
        hasattr(estimator, "transform") or hasattr(estimator, "fit_transform")
    ):
        transformer_tags = TransformerTags(
            preserves_dtype=old_tags["preserves_dtype"],
        )
    else:
        transformer_tags = None
    estimator_type = getattr(estimator, "_estimator_type", None)
    if estimator_type == "classifier":
        classifier_tags = ClassifierTags(
            poor_score=old_tags["poor_score"],
            multi_class=not old_tags["binary_only"],
            multi_label=old_tags["multilabel"],
        )
    else:
        classifier_tags = None
    if estimator_type == "regressor":
        regressor_tags = RegressorTags(
            poor_score=old_tags["poor_score"],
        )
    else:
        regressor_tags = None
    return Tags(
        estimator_type=estimator_type,
        target_tags=target_tags,
        transformer_tags=transformer_tags,
        classifier_tags=classifier_tags,
        regressor_tags=regressor_tags,
        input_tags=input_tags,
        array_api_support=old_tags["array_api_support"],
        no_validation=old_tags["no_validation"],
        non_deterministic=old_tags["non_deterministic"],
        requires_fit=old_tags["requires_fit"],
        _skip_test=old_tags["_skip_test"],
    )


# TODO(1.7): Remove this function
def _to_old_tags(new_tags):
    """Utility function convert old tags (dictionary) to new tags (dataclass)."""
    if new_tags.classifier_tags:
        binary_only = not new_tags.classifier_tags.multi_class
        multilabel = new_tags.classifier_tags.multi_label
        poor_score_clf = new_tags.classifier_tags.poor_score
    else:
        binary_only = False
        multilabel = False
        poor_score_clf = False

    if new_tags.regressor_tags:
        poor_score_reg = new_tags.regressor_tags.poor_score
    else:
        poor_score_reg = False

    if new_tags.transformer_tags:
        preserves_dtype = new_tags.transformer_tags.preserves_dtype
    else:
        preserves_dtype = ["float64"]

    tags = {
        "allow_nan": new_tags.input_tags.allow_nan,
        "array_api_support": new_tags.array_api_support,
        "binary_only": binary_only,
        "multilabel": multilabel,
        "multioutput": new_tags.target_tags.multi_output,
        "multioutput_only": (
            not new_tags.target_tags.single_output and new_tags.target_tags.multi_output
        ),
        "no_validation": new_tags.no_validation,
        "non_deterministic": new_tags.non_deterministic,
        "pairwise": new_tags.input_tags.pairwise,
        "preserves_dtype": preserves_dtype,
        "poor_score": poor_score_clf or poor_score_reg,
        "requires_fit": new_tags.requires_fit,
        "requires_positive_X": new_tags.input_tags.positive_only,
        "requires_y": new_tags.target_tags.required,
        "requires_positive_y": new_tags.target_tags.positive_only,
        "_skip_test": new_tags._skip_test,
        "stateless": new_tags.requires_fit,
    }
    X_types = []
    if new_tags.input_tags.one_d_array:
        X_types.append("1darray")
    if new_tags.input_tags.two_d_array:
        X_types.append("2darray")
    if new_tags.input_tags.three_d_array:
        X_types.append("3darray")
    if new_tags.input_tags.sparse:
        X_types.append("sparse")
    if new_tags.input_tags.categorical:
        X_types.append("categorical")
    if new_tags.input_tags.string:
        X_types.append("string")
    if new_tags.input_tags.dict:
        X_types.append("dict")
    if new_tags.target_tags.one_d_labels:
        X_types.append("1dlabels")
    if new_tags.target_tags.two_d_labels:
        X_types.append("2dlabels")
    tags["X_types"] = X_types
    return tags
