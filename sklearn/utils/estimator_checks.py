"""Various utilities to check the compatibility of estimators with scikit-learn API."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
from __future__ import annotations

import pickle
import re
import textwrap
import warnings
from contextlib import nullcontext
from copy import deepcopy
from functools import partial, wraps
from inspect import signature
from numbers import Integral, Real
from typing import Callable, Literal

import joblib
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

from sklearn.base import (
    BaseEstimator,
    BiclusterMixin,
    ClassifierMixin,
    ClassNamePrefixFeaturesOutMixin,
    DensityMixin,
    MetaEstimatorMixin,
    MultiOutputMixin,
    OneToOneFeatureMixin,
    OutlierMixin,
    RegressorMixin,
    TransformerMixin,
)

from .. import config_context
from ..base import (
    ClusterMixin,
    clone,
    is_classifier,
    is_outlier_detector,
    is_regressor,
)
from ..datasets import (
    load_iris,
    make_blobs,
    make_classification,
    make_multilabel_classification,
    make_regression,
)
from ..exceptions import (
    DataConversionWarning,
    EstimatorCheckFailedWarning,
    NotFittedError,
    SkipTestWarning,
)
from ..linear_model._base import LinearClassifierMixin
from ..metrics import accuracy_score, adjusted_rand_score, f1_score
from ..metrics.pairwise import linear_kernel, pairwise_distances, rbf_kernel
from ..model_selection import LeaveOneGroupOut, ShuffleSplit, train_test_split
from ..model_selection._validation import _safe_split
from ..pipeline import make_pipeline
from ..preprocessing import StandardScaler, scale
from ..utils import _safe_indexing
from ..utils._array_api import (
    _atol_for_type,
    _convert_to_numpy,
    get_namespace,
    yield_namespace_device_dtype_combinations,
)
from ..utils._array_api import device as array_device
from ..utils._param_validation import (
    InvalidParameterError,
    generate_invalid_param_val,
    make_constraint,
)
from . import shuffle
from ._missing import is_scalar_nan
from ._param_validation import Interval, StrOptions, validate_params
from ._tags import (
    ClassifierTags,
    InputTags,
    RegressorTags,
    TargetTags,
    TransformerTags,
    get_tags,
)
from ._test_common.instance_generator import (
    CROSS_DECOMPOSITION,
    _get_check_estimator_ids,
    _yield_instances_for_check,
)
from ._testing import (
    SkipTest,
    _array_api_for_tests,
    _get_args,
    assert_allclose,
    assert_allclose_dense_sparse,
    assert_array_almost_equal,
    assert_array_equal,
    assert_array_less,
    create_memmap_backed_data,
    ignore_warnings,
    raises,
    set_random_state,
)
from .validation import _num_samples, check_is_fitted, has_fit_parameter

REGRESSION_DATASET = None


def _raise_for_missing_tags(estimator, tag_name, Mixin):
    tags = get_tags(estimator)
    estimator_type = Mixin.__name__.replace("Mixin", "")
    if getattr(tags, tag_name) is None:
        raise RuntimeError(
            f"Estimator {estimator.__class__.__name__} seems to be a {estimator_type},"
            f" but the `{tag_name}` tag is not set. Either set the tag manually"
            f" or inherit from the {Mixin.__name__}. Note that the order of inheritance"
            f" matters, the {Mixin.__name__} should come before BaseEstimator."
        )


def _yield_api_checks(estimator):
    if not isinstance(estimator, BaseEstimator):
        warnings.warn(
            f"Estimator {estimator.__class__.__name__} does not inherit from"
            " `sklearn.base.BaseEstimator`. This might lead to unexpected behavior, or"
            " even errors when collecting tests.",
            category=UserWarning,
        )

    tags = get_tags(estimator)
    yield check_estimator_cloneable
    yield check_estimator_tags_renamed
    yield check_valid_tag_types
    yield check_estimator_repr
    yield check_no_attributes_set_in_init
    yield check_fit_score_takes_y
    yield check_estimators_overwrite_params
    yield check_dont_overwrite_parameters
    yield check_estimators_fit_returns_self
    yield check_readonly_memmap_input
    if tags.requires_fit:
        yield check_estimators_unfitted
    yield check_do_not_raise_errors_in_init_or_set_params
    yield check_n_features_in_after_fitting
    yield check_mixin_order
    yield check_positive_only_tag_during_fit


def _yield_checks(estimator):
    name = estimator.__class__.__name__
    tags = get_tags(estimator)

    yield check_estimators_dtypes
    if has_fit_parameter(estimator, "sample_weight"):
        yield check_sample_weights_pandas_series
        yield check_sample_weights_not_an_array
        yield check_sample_weights_list
        if not tags.input_tags.pairwise:
            # We skip pairwise because the data is not pairwise
            yield check_sample_weights_shape
            yield check_sample_weights_not_overwritten
            yield check_sample_weight_equivalence_on_dense_data
            if tags.input_tags.sparse:
                yield check_sample_weight_equivalence_on_sparse_data

    # Check that all estimator yield informative messages when
    # trained on empty datasets
    if not tags.no_validation:
        yield check_complex_data
        yield check_dtype_object
        yield check_estimators_empty_data_messages

    if name not in CROSS_DECOMPOSITION:
        # cross-decomposition's "transform" returns X and Y
        yield check_pipeline_consistency

    if not tags.input_tags.allow_nan and not tags.no_validation:
        # Test that all estimators check their input for NaN's and infs
        yield check_estimators_nan_inf

    if tags.input_tags.pairwise:
        # Check that pairwise estimator throws error on non-square input
        yield check_nonsquare_error

    if hasattr(estimator, "sparsify"):
        yield check_sparsify_coefficients

    yield check_estimator_sparse_tag
    yield check_estimator_sparse_array
    yield check_estimator_sparse_matrix

    # Test that estimators can be pickled, and once pickled
    # give the same answer as before.
    yield check_estimators_pickle
    yield partial(check_estimators_pickle, readonly_memmap=True)

    if tags.array_api_support:
        for check in _yield_array_api_checks(estimator):
            yield check

    yield check_f_contiguous_array_estimator


def _yield_classifier_checks(classifier):
    _raise_for_missing_tags(classifier, "classifier_tags", ClassifierMixin)
    tags = get_tags(classifier)

    # test classifiers can handle non-array data and pandas objects
    yield check_classifier_data_not_an_array
    # test classifiers trained on a single label always return this label
    yield check_classifiers_one_label
    yield check_classifiers_one_label_sample_weights
    yield check_classifiers_classes
    yield check_estimators_partial_fit_n_features
    if tags.target_tags.multi_output:
        yield check_classifier_multioutput
    # basic consistency testing
    yield check_classifiers_train
    yield partial(check_classifiers_train, readonly_memmap=True)
    yield partial(check_classifiers_train, readonly_memmap=True, X_dtype="float32")
    yield check_classifiers_regression_target
    if tags.classifier_tags.multi_label:
        yield check_classifiers_multilabel_representation_invariance
        yield check_classifiers_multilabel_output_format_predict
        yield check_classifiers_multilabel_output_format_predict_proba
        yield check_classifiers_multilabel_output_format_decision_function
    if not tags.no_validation:
        yield check_supervised_y_no_nan
        if tags.target_tags.single_output:
            yield check_supervised_y_2d
    if "class_weight" in classifier.get_params().keys():
        yield check_class_weight_classifiers

    yield check_non_transformer_estimators_n_iter
    # test if predict_proba is a monotonic transformation of decision_function
    yield check_decision_proba_consistency

    if isinstance(classifier, LinearClassifierMixin):
        if "class_weight" in classifier.get_params().keys():
            yield check_class_weight_balanced_linear_classifier
    if (
        isinstance(classifier, LinearClassifierMixin)
        and "class_weight" in classifier.get_params().keys()
    ):
        yield check_class_weight_balanced_linear_classifier

    if not tags.classifier_tags.multi_class:
        yield check_classifier_not_supporting_multiclass


def _yield_regressor_checks(regressor):
    _raise_for_missing_tags(regressor, "regressor_tags", RegressorMixin)
    tags = get_tags(regressor)
    # TODO: test with intercept
    # TODO: test with multiple responses
    # basic testing
    yield check_regressors_train
    yield partial(check_regressors_train, readonly_memmap=True)
    yield partial(check_regressors_train, readonly_memmap=True, X_dtype="float32")
    yield check_regressor_data_not_an_array
    yield check_estimators_partial_fit_n_features
    if tags.target_tags.multi_output:
        yield check_regressor_multioutput
    yield check_regressors_no_decision_function
    if not tags.no_validation and tags.target_tags.single_output:
        yield check_supervised_y_2d
    yield check_supervised_y_no_nan
    name = regressor.__class__.__name__
    if name != "CCA":
        # check that the regressor handles int input
        yield check_regressors_int
    yield check_non_transformer_estimators_n_iter


def _yield_transformer_checks(transformer):
    _raise_for_missing_tags(transformer, "transformer_tags", TransformerMixin)
    tags = get_tags(transformer)
    # All transformers should either deal with sparse data or raise an
    # exception with type TypeError and an intelligible error message
    if not tags.no_validation:
        yield check_transformer_data_not_an_array
    # these don't actually fit the data, so don't raise errors
    yield check_transformer_general
    if tags.transformer_tags.preserves_dtype:
        yield check_transformer_preserve_dtypes
    yield partial(check_transformer_general, readonly_memmap=True)
    if get_tags(transformer).requires_fit:
        yield check_transformers_unfitted
    else:
        yield check_transformers_unfitted_stateless
    # Dependent on external solvers and hence accessing the iter
    # param is non-trivial.
    external_solver = [
        "Isomap",
        "KernelPCA",
        "LocallyLinearEmbedding",
        "LogisticRegressionCV",
        "BisectingKMeans",
    ]

    name = transformer.__class__.__name__
    if name not in external_solver:
        yield check_transformer_n_iter


def _yield_clustering_checks(clusterer):
    yield check_clusterer_compute_labels_predict
    name = clusterer.__class__.__name__
    if name not in ("WardAgglomeration", "FeatureAgglomeration"):
        # this is clustering on the features
        # let's not test that here.
        yield check_clustering
        yield partial(check_clustering, readonly_memmap=True)
        yield check_estimators_partial_fit_n_features
    if not hasattr(clusterer, "transform"):
        yield check_non_transformer_estimators_n_iter


def _yield_outliers_checks(estimator):
    # checks for the contamination parameter
    if hasattr(estimator, "contamination"):
        yield check_outlier_contamination

    # checks for outlier detectors that have a fit_predict method
    if hasattr(estimator, "fit_predict"):
        yield check_outliers_fit_predict

    # checks for estimators that can be used on a test set
    if hasattr(estimator, "predict"):
        yield check_outliers_train
        yield partial(check_outliers_train, readonly_memmap=True)
        # test outlier detectors can handle non-array data
        yield check_classifier_data_not_an_array
    yield check_non_transformer_estimators_n_iter


def _yield_array_api_checks(estimator):
    for (
        array_namespace,
        device,
        dtype_name,
    ) in yield_namespace_device_dtype_combinations():
        yield partial(
            check_array_api_input,
            array_namespace=array_namespace,
            dtype_name=dtype_name,
            device=device,
        )


def _yield_all_checks(estimator, legacy: bool):
    name = estimator.__class__.__name__
    tags = get_tags(estimator)
    if not tags.input_tags.two_d_array:
        warnings.warn(
            "Can't test estimator {} which requires input  of type {}".format(
                name, tags.input_tags
            ),
            SkipTestWarning,
        )
        return
    if tags._skip_test:
        warnings.warn(
            "Explicit SKIP via _skip_test tag for estimator {}.".format(name),
            SkipTestWarning,
        )
        return

    for check in _yield_api_checks(estimator):
        yield check

    if not legacy:
        return  # pragma: no cover

    for check in _yield_checks(estimator):
        yield check
    if is_classifier(estimator):
        for check in _yield_classifier_checks(estimator):
            yield check
    if is_regressor(estimator):
        for check in _yield_regressor_checks(estimator):
            yield check
    if hasattr(estimator, "transform"):
        for check in _yield_transformer_checks(estimator):
            yield check
    if isinstance(estimator, ClusterMixin):
        for check in _yield_clustering_checks(estimator):
            yield check
    if is_outlier_detector(estimator):
        for check in _yield_outliers_checks(estimator):
            yield check
    yield check_parameters_default_constructible
    if not tags.non_deterministic:
        yield check_methods_sample_order_invariance
        yield check_methods_subset_invariance
    yield check_fit2d_1sample
    yield check_fit2d_1feature
    yield check_get_params_invariance
    yield check_set_params
    yield check_dict_unchanged
    yield check_fit_idempotent
    yield check_fit_check_is_fitted
    if not tags.no_validation:
        yield check_n_features_in
        yield check_fit1d
        yield check_fit2d_predict1d
        if tags.target_tags.required:
            yield check_requires_y_none
    if tags.input_tags.positive_only:
        yield check_fit_non_negative


def _check_name(check):
    if hasattr(check, "__wrapped__"):
        return _check_name(check.__wrapped__)
    return check.func.__name__ if isinstance(check, partial) else check.__name__


def _maybe_mark(
    estimator,
    check,
    expected_failed_checks: dict[str, str] | None = None,
    mark: Literal["xfail", "skip", None] = None,
    pytest=None,
):
    """Mark the test as xfail or skip if needed.

    Parameters
    ----------
    estimator : estimator object
        Estimator instance for which to generate checks.
    check : partial or callable
        Check to be marked.
    expected_failed_checks : dict[str, str], default=None
        Dictionary of the form {check_name: reason} for checks that are expected to
        fail.
    mark : "xfail" or "skip" or None
        Whether to mark the check as xfail or skip.
    pytest : pytest module, default=None
        Pytest module to use to mark the check. This is only needed if ``mark`` is
        `"xfail"`. Note that one can run `check_estimator` without having `pytest`
        installed. This is used in combination with `parametrize_with_checks` only.
    """
    should_be_marked, reason = _should_be_skipped_or_marked(
        estimator, check, expected_failed_checks
    )
    if not should_be_marked or mark is None:
        return estimator, check

    estimator_name = estimator.__class__.__name__
    if mark == "xfail":
        return pytest.param(estimator, check, marks=pytest.mark.xfail(reason=reason))
    else:

        @wraps(check)
        def wrapped(*args, **kwargs):
            raise SkipTest(
                f"Skipping {_check_name(check)} for {estimator_name}: {reason}"
            )

        return estimator, wrapped


def _should_be_skipped_or_marked(
    estimator, check, expected_failed_checks: dict[str, str] | None = None
) -> tuple[bool, str]:
    """Check whether a check should be skipped or marked as xfail.

    Parameters
    ----------
    estimator : estimator object
        Estimator instance for which to generate checks.
    check : partial or callable
        Check to be marked.
    expected_failed_checks : dict[str, str], default=None
        Dictionary of the form {check_name: reason} for checks that are expected to
        fail.

    Returns
    -------
    should_be_marked : bool
        Whether the check should be marked as xfail or skipped.
    reason : str
        Reason for skipping the check.
    """

    expected_failed_checks = expected_failed_checks or {}

    check_name = _check_name(check)
    if check_name in expected_failed_checks:
        return True, expected_failed_checks[check_name]

    return False, "Check is not expected to fail"


def estimator_checks_generator(
    estimator,
    *,
    legacy: bool = True,
    expected_failed_checks: dict[str, str] | None = None,
    mark: Literal["xfail", "skip", None] = None,
):
    """Iteratively yield all check callables for an estimator.

    This function is used by
    :func:`~sklearn.utils.estimator_checks.parametrize_with_checks` and
    :func:`~sklearn.utils.estimator_checks.check_estimator` to yield all check callables
    for an estimator. In most cases, these functions should be used instead. When
    implementing a custom equivalent, please refer to their source code to
    understand how `estimator_checks_generator` is intended to be used.

    .. versionadded:: 1.6

    Parameters
    ----------
    estimator : estimator object
        Estimator instance for which to generate checks.
    legacy : bool, default=True
        Whether to include legacy checks. Over time we remove checks from this category
        and move them into their specific category.
    expected_failed_checks : dict[str, str], default=None
        Dictionary of the form {check_name: reason} for checks that are expected to
        fail.
    mark : {"xfail", "skip"} or None, default=None
        Whether to mark the checks that are expected to fail as
        xfail(`pytest.mark.xfail`) or skip. Marking a test as "skip" is done via
        wrapping the check in a function that raises a
        :class:`~sklearn.exceptions.SkipTest` exception.

    Returns
    -------
    estimator_checks_generator : generator
        Generator that yields (estimator, check) tuples.
    """
    if mark == "xfail":
        import pytest
    else:
        pytest = None  # type: ignore[assignment]

    name = type(estimator).__name__
    # First check that the estimator is cloneable which is needed for the rest
    # of the checks to run
    yield estimator, partial(check_estimator_cloneable, name)
    for check in _yield_all_checks(estimator, legacy=legacy):
        check_with_name = partial(check, name)
        for check_instance in _yield_instances_for_check(check, estimator):
            yield _maybe_mark(
                check_instance,
                check_with_name,
                expected_failed_checks=expected_failed_checks,
                mark=mark,
                pytest=pytest,
            )


def parametrize_with_checks(
    estimators,
    *,
    legacy: bool = True,
    expected_failed_checks: Callable | None = None,
):
    """Pytest specific decorator for parametrizing estimator checks.

    Checks are categorised into the following groups:

    - API checks: a set of checks to ensure API compatibility with scikit-learn.
      Refer to https://scikit-learn.org/dev/developers/develop.html a requirement of
      scikit-learn estimators.
    - legacy: a set of checks which gradually will be grouped into other categories.

    The `id` of each check is set to be a pprint version of the estimator
    and the name of the check with its keyword arguments.
    This allows to use `pytest -k` to specify which tests to run::

        pytest test_check_estimators.py -k check_estimators_fit_returns_self

    Parameters
    ----------
    estimators : list of estimators instances
        Estimators to generated checks for.

        .. versionchanged:: 0.24
           Passing a class was deprecated in version 0.23, and support for
           classes was removed in 0.24. Pass an instance instead.

        .. versionadded:: 0.24


    legacy : bool, default=True
        Whether to include legacy checks. Over time we remove checks from this category
        and move them into their specific category.

        .. versionadded:: 1.6

    expected_failed_checks : callable, default=None
        A callable that takes an estimator as input and returns a dictionary of the
        form::

            {
                "check_name": "my reason",
            }

        Where `"check_name"` is the name of the check, and `"my reason"` is why
        the check fails. These tests will be marked as xfail if the check fails.


        .. versionadded:: 1.6

    Returns
    -------
    decorator : `pytest.mark.parametrize`

    See Also
    --------
    check_estimator : Check if estimator adheres to scikit-learn conventions.

    Examples
    --------
    >>> from sklearn.utils.estimator_checks import parametrize_with_checks
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.tree import DecisionTreeRegressor

    >>> @parametrize_with_checks([LogisticRegression(),
    ...                           DecisionTreeRegressor()])
    ... def test_sklearn_compatible_estimator(estimator, check):
    ...     check(estimator)

    """
    import pytest

    if any(isinstance(est, type) for est in estimators):
        msg = (
            "Passing a class was deprecated in version 0.23 "
            "and isn't supported anymore from 0.24."
            "Please pass an instance instead."
        )
        raise TypeError(msg)

    def _checks_generator(estimators, legacy, expected_failed_checks):
        for estimator in estimators:
            args = {"estimator": estimator, "legacy": legacy, "mark": "xfail"}
            if callable(expected_failed_checks):
                args["expected_failed_checks"] = expected_failed_checks(estimator)
            yield from estimator_checks_generator(**args)

    return pytest.mark.parametrize(
        "estimator, check",
        _checks_generator(estimators, legacy, expected_failed_checks),
        ids=_get_check_estimator_ids,
    )


@validate_params(
    {
        "generate_only": ["boolean"],
        "legacy": ["boolean"],
        "expected_failed_checks": [dict, None],
        "on_skip": [StrOptions({"warn"}), None],
        "on_fail": [StrOptions({"raise", "warn"}), None],
        "callback": [callable, None],
    },
    prefer_skip_nested_validation=False,
)
def check_estimator(
    estimator=None,
    generate_only=False,
    *,
    legacy: bool = True,
    expected_failed_checks: dict[str, str] | None = None,
    on_skip: Literal["warn"] | None = "warn",
    on_fail: Literal["raise", "warn"] | None = "raise",
    callback: Callable | None = None,
):
    """Check if estimator adheres to scikit-learn conventions.

    This function will run an extensive test-suite for input validation,
    shapes, etc, making sure that the estimator complies with `scikit-learn`
    conventions as detailed in :ref:`rolling_your_own_estimator`.
    Additional tests for classifiers, regressors, clustering or transformers
    will be run if the Estimator class inherits from the corresponding mixin
    from sklearn.base.

    scikit-learn also provides a pytest specific decorator,
    :func:`~sklearn.utils.estimator_checks.parametrize_with_checks`, making it
    easier to test multiple estimators.

    Checks are categorised into the following groups:

    - API checks: a set of checks to ensure API compatibility with scikit-learn.
      Refer to https://scikit-learn.org/dev/developers/develop.html a requirement of
      scikit-learn estimators.
    - legacy: a set of checks which gradually will be grouped into other categories.

    Parameters
    ----------
    estimator : estimator object
        Estimator instance to check.

    generate_only : bool, default=False
        When `False`, checks are evaluated when `check_estimator` is called.
        When `True`, `check_estimator` returns a generator that yields
        (estimator, check) tuples. The check is run by calling
        `check(estimator)`.

        .. versionadded:: 0.22

        .. deprecated:: 1.6
            `generate_only` will be removed in 1.8. Use
            :func:`~sklearn.utils.estimator_checks.estimator_checks_generator` instead.

    legacy : bool, default=True
        Whether to include legacy checks. Over time we remove checks from this category
        and move them into their specific category.

        .. versionadded:: 1.6

    expected_failed_checks : dict, default=None
        A dictionary of the form::

            {
                "check_name": "this check is expected to fail because ...",
            }

        Where `"check_name"` is the name of the check, and `"my reason"` is why
        the check fails.

        .. versionadded:: 1.6

    on_skip : "warn", None, default="warn"
        This parameter controls what happens when a check is skipped.

        - "warn": A :class:`~sklearn.exceptions.SkipTestWarning` is logged
          and running tests continue.
        - None: No warning is logged and running tests continue.

        .. versionadded:: 1.6

    on_fail : {"raise", "warn"}, None, default="raise"
        This parameter controls what happens when a check fails.

        - "raise": The exception raised by the first failing check is raised and
          running tests are aborted. This does not included tests that are expected
          to fail.
        - "warn": A :class:`~sklearn.exceptions.EstimatorCheckFailedWarning` is logged
          and running tests continue.
        - None: No exception is raised and no warning is logged.

        Note that if ``on_fail != "raise"``, no exception is raised, even if the checks
        fail. You'd need to inspect the return result of ``check_estimator`` to check
        if any checks failed.

        .. versionadded:: 1.6

    callback : callable, or None, default=None
        This callback will be called with the estimator and the check name,
        the exception (if any), the status of the check (xfail, failed, skipped,
        passed), and the reason for the expected failure if the check is
        expected to fail. The callable's signature needs to be::

            def callback(
                estimator,
                check_name: str,
                exception: Exception,
                status: Literal["xfail", "failed", "skipped", "passed"],
                expected_to_fail: bool,
                expected_to_fail_reason: str,
            )

        ``callback`` cannot be provided together with ``on_fail="raise"``.

        .. versionadded:: 1.6

    Returns
    -------
    test_results : list
        List of dictionaries with the results of the failing tests, of the form::

            {
                "estimator": estimator,
                "check_name": check_name,
                "exception": exception,
                "status": status (one of "xfail", "failed", "skipped", "passed"),
                "expected_to_fail": expected_to_fail,
                "expected_to_fail_reason": expected_to_fail_reason,
            }

    estimator_checks_generator : generator
        Generator that yields (estimator, check) tuples. Returned when
        `generate_only=True`.

        ..
            TODO(1.8): remove return value

        .. deprecated:: 1.6
            ``generate_only`` will be removed in 1.8. Use
            :func:`~sklearn.utils.estimator_checks.estimator_checks_generator` instead.

    Raises
    ------
    Exception
        If ``on_fail="raise"``, the exception raised by the first failing check is
        raised and running tests are aborted.

        Note that if ``on_fail != "raise"``, no exception is raised, even if the checks
        fail. You'd need to inspect the return result of ``check_estimator`` to check
        if any checks failed.

    See Also
    --------
    parametrize_with_checks : Pytest specific decorator for parametrizing estimator
        checks.
    estimator_checks_generator : Generator that yields (estimator, check) tuples.

    Examples
    --------
    >>> from sklearn.utils.estimator_checks import check_estimator
    >>> from sklearn.linear_model import LogisticRegression
    >>> check_estimator(LogisticRegression())
    [...]
    """
    if isinstance(estimator, type):
        msg = (
            "Passing a class was deprecated in version 0.23 "
            "and isn't supported anymore from 0.24."
            "Please pass an instance instead."
        )
        raise TypeError(msg)

    if on_fail == "raise" and callback is not None:
        raise ValueError("callback cannot be provided together with on_fail='raise'")

    name = type(estimator).__name__

    # TODO(1.8): remove generate_only
    if generate_only:
        warnings.warn(
            "`generate_only` is deprecated in 1.6 and will be removed in 1.8. "
            "Use :func:`~sklearn.utils.estimator_checks.estimator_checks_generator` "
            "instead.",
            FutureWarning,
        )
        return estimator_checks_generator(
            estimator, legacy=legacy, expected_failed_checks=None, mark="skip"
        )

    test_results = []

    for estimator, check in estimator_checks_generator(
        estimator,
        legacy=legacy,
        expected_failed_checks=expected_failed_checks,
        # Not marking tests to be skipped here, we run and simulate an xfail behavior
        mark=None,
    ):
        test_can_fail, reason = _should_be_skipped_or_marked(
            estimator, check, expected_failed_checks
        )
        try:
            check(estimator)
        except SkipTest as e:
            # We get here if the test raises SkipTest, which is expected in cases where
            # the check cannot run for instance if a required dependency is not
            # installed.
            check_result = {
                "estimator": estimator,
                "check_name": _check_name(check),
                "exception": e,
                "status": "skipped",
                "expected_to_fail": test_can_fail,
                "expected_to_fail_reason": reason,
            }
            if on_skip == "warn":
                warnings.warn(
                    f"Skipping check {_check_name(check)} for {name} because it raised "
                    f"{type(e).__name__}: {e}",
                    SkipTestWarning,
                )
        except Exception as e:
            if on_fail == "raise" and not test_can_fail:
                raise

            check_result = {
                "estimator": estimator,
                "check_name": _check_name(check),
                "exception": e,
                "expected_to_fail": test_can_fail,
                "expected_to_fail_reason": reason,
            }

            if test_can_fail:
                # This check failed, but could be expected to fail, therefore we mark it
                # as xfail.
                check_result["status"] = "xfail"
            else:
                check_result["status"] = "failed"

            if on_fail == "warn":
                warning = EstimatorCheckFailedWarning(**check_result)
                warnings.warn(warning)
        else:
            check_result = {
                "estimator": estimator,
                "check_name": _check_name(check),
                "exception": None,
                "status": "passed",
                "expected_to_fail": test_can_fail,
                "expected_to_fail_reason": reason,
            }

        test_results.append(check_result)

        if callback:
            callback(**check_result)

    return test_results


def _regression_dataset():
    global REGRESSION_DATASET
    if REGRESSION_DATASET is None:
        X, y = make_regression(
            n_samples=200,
            n_features=10,
            n_informative=1,
            bias=5.0,
            noise=20,
            random_state=42,
        )
        X = StandardScaler().fit_transform(X)
        REGRESSION_DATASET = X, y
    return REGRESSION_DATASET


class _NotAnArray:
    """An object that is convertible to an array.

    Parameters
    ----------
    data : array-like
        The data.
    """

    def __init__(self, data):
        self.data = np.asarray(data)

    def __array__(self, dtype=None, copy=None):
        return self.data

    def __array_function__(self, func, types, args, kwargs):
        if func.__name__ == "may_share_memory":
            return True
        raise TypeError("Don't want to call array_function {}!".format(func.__name__))


def _is_pairwise_metric(estimator):
    """Returns True if estimator accepts pairwise metric.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    Returns
    -------
    out : bool
        True if _pairwise is set to True and False otherwise.
    """
    metric = getattr(estimator, "metric", None)

    return bool(metric == "precomputed")


def _generate_sparse_data(X_csr):
    """Generate sparse matrices or arrays with {32,64}bit indices of diverse format.

    Parameters
    ----------
    X_csr: scipy.sparse.csr_matrix or scipy.sparse.csr_array
        Input in CSR format.

    Returns
    -------
    out: iter(Matrices) or iter(Arrays)
        In format['dok', 'lil', 'dia', 'bsr', 'csr', 'csc', 'coo',
        'coo_64', 'csc_64', 'csr_64']
    """

    assert X_csr.format == "csr"
    yield "csr", X_csr.copy()
    for sparse_format in ["dok", "lil", "dia", "bsr", "csc", "coo"]:
        yield sparse_format, X_csr.asformat(sparse_format)

    # Generate large indices matrix only if its supported by scipy
    X_coo = X_csr.asformat("coo")
    X_coo.row = X_coo.row.astype("int64")
    X_coo.col = X_coo.col.astype("int64")
    yield "coo_64", X_coo

    for sparse_format in ["csc", "csr"]:
        X = X_csr.asformat(sparse_format)
        X.indices = X.indices.astype("int64")
        X.indptr = X.indptr.astype("int64")
        yield sparse_format + "_64", X


@ignore_warnings(category=FutureWarning)
def check_supervised_y_no_nan(name, estimator_orig):
    # Checks that the Estimator targets are not NaN.
    estimator = clone(estimator_orig)
    rng = np.random.RandomState(888)
    X = rng.standard_normal(size=(10, 5))

    for value in [np.nan, np.inf]:
        y = np.full(10, value)
        y = _enforce_estimator_tags_y(estimator, y)

        module_name = estimator.__module__
        if module_name.startswith("sklearn.") and not (
            "test_" in module_name or module_name.endswith("_testing")
        ):
            # In scikit-learn we want the error message to mention the input
            # name and be specific about the kind of unexpected value.
            if np.isinf(value):
                match = (
                    r"Input (y|Y) contains infinity or a value too large for"
                    r" dtype\('float64'\)."
                )
            else:
                match = r"Input (y|Y) contains NaN."
        else:
            # Do not impose a particular error message to third-party libraries.
            match = None
        err_msg = (
            f"Estimator {name} should have raised error on fitting array y with inf"
            " value."
        )
        with raises(ValueError, match=match, err_msg=err_msg):
            estimator.fit(X, y)


def check_array_api_input(
    name,
    estimator_orig,
    array_namespace,
    device=None,
    dtype_name="float64",
    check_values=False,
):
    """Check that the estimator can work consistently with the Array API

    By default, this just checks that the types and shapes of the arrays are
    consistent with calling the same estimator with numpy arrays.

    When check_values is True, it also checks that calling the estimator on the
    array_api Array gives the same results as ndarrays.
    """
    xp = _array_api_for_tests(array_namespace, device)

    X, y = make_classification(random_state=42)
    X = X.astype(dtype_name, copy=False)

    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = _enforce_estimator_tags_y(estimator_orig, y)

    est = clone(estimator_orig)

    X_xp = xp.asarray(X, device=device)
    y_xp = xp.asarray(y, device=device)

    est.fit(X, y)

    array_attributes = {
        key: value for key, value in vars(est).items() if isinstance(value, np.ndarray)
    }

    est_xp = clone(est)
    with config_context(array_api_dispatch=True):
        est_xp.fit(X_xp, y_xp)
        input_ns = get_namespace(X_xp)[0].__name__

    # Fitted attributes which are arrays must have the same
    # namespace as the one of the training data.
    for key, attribute in array_attributes.items():
        est_xp_param = getattr(est_xp, key)
        with config_context(array_api_dispatch=True):
            attribute_ns = get_namespace(est_xp_param)[0].__name__
        assert attribute_ns == input_ns, (
            f"'{key}' attribute is in wrong namespace, expected {input_ns} "
            f"got {attribute_ns}"
        )

        assert array_device(est_xp_param) == array_device(X_xp)

        est_xp_param_np = _convert_to_numpy(est_xp_param, xp=xp)
        if check_values:
            assert_allclose(
                attribute,
                est_xp_param_np,
                err_msg=f"{key} not the same",
                atol=_atol_for_type(X.dtype),
            )
        else:
            assert attribute.shape == est_xp_param_np.shape
            assert attribute.dtype == est_xp_param_np.dtype

    # Check estimator methods, if supported, give the same results
    methods = (
        "score",
        "score_samples",
        "decision_function",
        "predict",
        "predict_log_proba",
        "predict_proba",
        "transform",
    )

    try:
        np.asarray(X_xp)
        np.asarray(y_xp)
        # TODO There are a few errors in SearchCV with array-api-strict because
        # we end up doing X[train_indices] where X is an array-api-strict array
        # and train_indices is a numpy array. array-api-strict insists
        # train_indices should be an array-api-strict array. On the other hand,
        # all the array API libraries (PyTorch, jax, CuPy) accept indexing with a
        # numpy array. This is probably not worth doing anything about for
        # now since array-api-strict seems a bit too strict ...
        numpy_asarray_works = xp.__name__ != "array_api_strict"

    except (TypeError, RuntimeError, ValueError):
        # PyTorch with CUDA device and CuPy raise TypeError consistently.
        # array-api-strict chose to raise RuntimeError instead. NumPy emits
        # a ValueError if `__array__` dunder does not return an array.
        # Exception type may need to be updated in the future for other libraries.
        numpy_asarray_works = False

    if numpy_asarray_works:
        # In this case, array_api_dispatch is disabled and we rely on np.asarray
        # being called to convert the non-NumPy inputs to NumPy arrays when needed.
        est_fitted_with_as_array = clone(est).fit(X_xp, y_xp)
        # We only do a smoke test for now, in order to avoid complicating the
        # test function even further.
        for method_name in methods:
            method = getattr(est_fitted_with_as_array, method_name, None)
            if method is None:
                continue

            if method_name == "score":
                method(X_xp, y_xp)
            else:
                method(X_xp)

    for method_name in methods:
        method = getattr(est, method_name, None)
        if method is None:
            continue

        if method_name == "score":
            result = method(X, y)
            with config_context(array_api_dispatch=True):
                result_xp = getattr(est_xp, method_name)(X_xp, y_xp)
            # score typically returns a Python float
            assert isinstance(result, float)
            assert isinstance(result_xp, float)
            if check_values:
                assert abs(result - result_xp) < _atol_for_type(X.dtype)
            continue
        else:
            result = method(X)
            with config_context(array_api_dispatch=True):
                result_xp = getattr(est_xp, method_name)(X_xp)

        with config_context(array_api_dispatch=True):
            result_ns = get_namespace(result_xp)[0].__name__
        assert result_ns == input_ns, (
            f"'{method}' output is in wrong namespace, expected {input_ns}, "
            f"got {result_ns}."
        )

        assert array_device(result_xp) == array_device(X_xp)
        result_xp_np = _convert_to_numpy(result_xp, xp=xp)

        if check_values:
            assert_allclose(
                result,
                result_xp_np,
                err_msg=f"{method} did not the return the same result",
                atol=_atol_for_type(X.dtype),
            )
        else:
            if hasattr(result, "shape"):
                assert result.shape == result_xp_np.shape
                assert result.dtype == result_xp_np.dtype

        if method_name == "transform" and hasattr(est, "inverse_transform"):
            inverse_result = est.inverse_transform(result)
            with config_context(array_api_dispatch=True):
                invese_result_xp = est_xp.inverse_transform(result_xp)
                inverse_result_ns = get_namespace(invese_result_xp)[0].__name__
            assert inverse_result_ns == input_ns, (
                "'inverse_transform' output is in wrong namespace, expected"
                f" {input_ns}, got {inverse_result_ns}."
            )

            assert array_device(invese_result_xp) == array_device(X_xp)

            invese_result_xp_np = _convert_to_numpy(invese_result_xp, xp=xp)
            if check_values:
                assert_allclose(
                    inverse_result,
                    invese_result_xp_np,
                    err_msg="inverse_transform did not the return the same result",
                    atol=_atol_for_type(X.dtype),
                )
            else:
                assert inverse_result.shape == invese_result_xp_np.shape
                assert inverse_result.dtype == invese_result_xp_np.dtype


def check_array_api_input_and_values(
    name,
    estimator_orig,
    array_namespace,
    device=None,
    dtype_name="float64",
):
    return check_array_api_input(
        name,
        estimator_orig,
        array_namespace=array_namespace,
        device=device,
        dtype_name=dtype_name,
        check_values=True,
    )


def check_estimator_sparse_tag(name, estimator_orig):
    """Check that estimator tag related with accepting sparse data is properly set."""
    estimator = clone(estimator_orig)

    rng = np.random.RandomState(0)
    n_samples = 15 if name == "SpectralCoclustering" else 40
    X = rng.uniform(size=(n_samples, 3))
    X[X < 0.6] = 0
    y = rng.randint(0, 3, size=n_samples)
    X = _enforce_estimator_tags_X(estimator, X)
    y = _enforce_estimator_tags_y(estimator, y)
    X = sparse.csr_array(X)

    tags = get_tags(estimator)
    if tags.input_tags.sparse:
        try:
            estimator.fit(X, y)  # should pass
        except Exception as e:
            err_msg = (
                f"Estimator {name} raised an exception. "
                f"The tag self.input_tags.sparse={tags.input_tags.sparse} "
                "might not be consistent with the estimator's ability to "
                "handle sparse data (i.e. controlled by the parameter `accept_sparse`"
                " in `validate_data` or `check_array` functions)."
            )
            raise AssertionError(err_msg) from e
    else:
        err_msg = (
            f"Estimator {name} raised an exception. "
            "The estimator failed when fitted on sparse data in accordance "
            f"with its tag self.input_tags.sparse={tags.input_tags.sparse} "
            "but didn't raise the appropriate error: error message should "
            "state explicitly that sparse input is not supported if this is "
            "not the case, e.g. by using check_array(X, accept_sparse=False)."
        )
        try:
            estimator.fit(X, y)  # should fail with appropriate error
        except (ValueError, TypeError) as e:
            if re.search("[Ss]parse", str(e)):
                # Got the right error type and mentioning sparse issue
                return
            raise AssertionError(err_msg) from e
        except Exception as e:
            raise AssertionError(err_msg) from e
        raise AssertionError(
            f"Estimator {name} didn't fail when fitted on sparse data "
            "but should have according to its tag "
            f"self.input_tags.sparse={tags.input_tags.sparse}. "
            f"The tag is inconsistent and must be fixed."
        )


def _check_estimator_sparse_container(name, estimator_orig, sparse_type):
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(40, 3))
    X[X < 0.6] = 0
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = (4 * rng.uniform(size=X.shape[0])).astype(np.int32)
    # catch deprecation warnings
    with ignore_warnings(category=FutureWarning):
        estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)
    tags = get_tags(estimator_orig)
    for matrix_format, X in _generate_sparse_data(sparse_type(X)):
        # catch deprecation warnings
        with ignore_warnings(category=FutureWarning):
            estimator = clone(estimator_orig)
            if name in ["Scaler", "StandardScaler"]:
                estimator.set_params(with_mean=False)
        # fit and predict
        if "64" in matrix_format:
            err_msg = (
                f"Estimator {name} doesn't seem to support {matrix_format} "
                "matrix, and is not failing gracefully, e.g. by using "
                "check_array(X, accept_large_sparse=False)."
            )
        else:
            err_msg = (
                f"Estimator {name} doesn't seem to fail gracefully on sparse "
                "data: error message should state explicitly that sparse "
                "input is not supported if this is not the case, e.g. by using "
                "check_array(X, accept_sparse=False)."
            )
        with raises(
            (TypeError, ValueError),
            match=["sparse", "Sparse"],
            may_pass=True,
            err_msg=err_msg,
        ):
            with ignore_warnings(category=FutureWarning):
                estimator.fit(X, y)
            if hasattr(estimator, "predict"):
                pred = estimator.predict(X)
                if tags.target_tags.multi_output and not tags.target_tags.single_output:
                    assert pred.shape == (X.shape[0], 1)
                else:
                    assert pred.shape == (X.shape[0],)
            if hasattr(estimator, "predict_proba"):
                probs = estimator.predict_proba(X)
                if not tags.classifier_tags.multi_class:
                    expected_probs_shape = (X.shape[0], 2)
                else:
                    expected_probs_shape = (X.shape[0], 4)
                assert probs.shape == expected_probs_shape


def check_estimator_sparse_matrix(name, estimator_orig):
    _check_estimator_sparse_container(name, estimator_orig, sparse.csr_matrix)


def check_estimator_sparse_array(name, estimator_orig):
    _check_estimator_sparse_container(name, estimator_orig, sparse.csr_array)


def check_f_contiguous_array_estimator(name, estimator_orig):
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/23988
    # https://github.com/scikit-learn/scikit-learn/issues/24013
    estimator = clone(estimator_orig)

    rng = np.random.RandomState(0)
    X = 3 * rng.uniform(size=(20, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    X = np.asfortranarray(X)
    y = X[:, 0].astype(int)
    y = _enforce_estimator_tags_y(estimator_orig, y)

    estimator.fit(X, y)

    if hasattr(estimator, "transform"):
        estimator.transform(X)

    if hasattr(estimator, "predict"):
        estimator.predict(X)


@ignore_warnings(category=FutureWarning)
def check_sample_weights_pandas_series(name, estimator_orig):
    # check that estimators will accept a 'sample_weight' parameter of
    # type pandas.Series in the 'fit' function.
    estimator = clone(estimator_orig)
    try:
        import pandas as pd

        X = np.array(
            [
                [1, 1],
                [1, 2],
                [1, 3],
                [1, 4],
                [2, 1],
                [2, 2],
                [2, 3],
                [2, 4],
                [3, 1],
                [3, 2],
                [3, 3],
                [3, 4],
            ]
        )
        X = pd.DataFrame(_enforce_estimator_tags_X(estimator_orig, X), copy=False)
        y = pd.Series([1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2])
        weights = pd.Series([1] * 12)
        if (
            not get_tags(estimator).target_tags.single_output
            and get_tags(estimator).target_tags.multi_output
        ):
            y = pd.DataFrame(y, copy=False)
        try:
            estimator.fit(X, y, sample_weight=weights)
        except ValueError:
            raise ValueError(
                "Estimator {0} raises error if "
                "'sample_weight' parameter is of "
                "type pandas.Series".format(name)
            )
    except ImportError:
        raise SkipTest(
            "pandas is not installed: not testing for "
            "input of type pandas.Series to class weight."
        )


@ignore_warnings(category=(FutureWarning))
def check_sample_weights_not_an_array(name, estimator_orig):
    # check that estimators will accept a 'sample_weight' parameter of
    # type _NotAnArray in the 'fit' function.
    estimator = clone(estimator_orig)
    X = np.array(
        [
            [1, 1],
            [1, 2],
            [1, 3],
            [1, 4],
            [2, 1],
            [2, 2],
            [2, 3],
            [2, 4],
            [3, 1],
            [3, 2],
            [3, 3],
            [3, 4],
        ]
    )
    X = _NotAnArray(_enforce_estimator_tags_X(estimator_orig, X))
    y = _NotAnArray([1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2])
    weights = _NotAnArray([1] * 12)
    tags = get_tags(estimator)
    if not tags.target_tags.single_output and tags.target_tags.multi_output:
        y = _NotAnArray(y.data.reshape(-1, 1))
    estimator.fit(X, y, sample_weight=weights)


@ignore_warnings(category=(FutureWarning))
def check_sample_weights_list(name, estimator_orig):
    # check that estimators will accept a 'sample_weight' parameter of
    # type list in the 'fit' function.
    estimator = clone(estimator_orig)
    rnd = np.random.RandomState(0)
    n_samples = 30
    X = _enforce_estimator_tags_X(estimator_orig, rnd.uniform(size=(n_samples, 3)))
    y = np.arange(n_samples) % 3
    y = _enforce_estimator_tags_y(estimator, y)
    sample_weight = [3] * n_samples
    # Test that estimators don't raise any exception
    estimator.fit(X, y, sample_weight=sample_weight)


@ignore_warnings(category=FutureWarning)
def check_sample_weights_shape(name, estimator_orig):
    # check that estimators raise an error if sample_weight
    # shape mismatches the input
    estimator = clone(estimator_orig)
    X = np.array(
        [
            [1, 3],
            [1, 3],
            [1, 3],
            [1, 3],
            [2, 1],
            [2, 1],
            [2, 1],
            [2, 1],
            [3, 3],
            [3, 3],
            [3, 3],
            [3, 3],
            [4, 1],
            [4, 1],
            [4, 1],
            [4, 1],
        ]
    )
    y = np.array([1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2])
    y = _enforce_estimator_tags_y(estimator, y)

    estimator.fit(X, y, sample_weight=np.ones(len(y)))

    with raises(ValueError):
        estimator.fit(X, y, sample_weight=np.ones(2 * len(y)))

    with raises(ValueError):
        estimator.fit(X, y, sample_weight=np.ones((len(y), 2)))


@ignore_warnings(category=FutureWarning)
def _check_sample_weight_equivalence(name, estimator_orig, sparse_container):
    # check that setting sample_weight to zero / integer is equivalent
    # to removing / repeating corresponding samples.
    estimator_weighted = clone(estimator_orig)
    estimator_repeated = clone(estimator_orig)
    set_random_state(estimator_weighted, random_state=0)
    set_random_state(estimator_repeated, random_state=0)

    rng = np.random.RandomState(42)
    n_samples = 15
    X = rng.rand(n_samples, n_samples * 2)
    y = rng.randint(0, 3, size=n_samples)
    # Use random integers (including zero) as weights.
    sw = rng.randint(0, 5, size=n_samples)

    X_weighted = X
    y_weighted = y
    # repeat samples according to weights
    X_repeated = X_weighted.repeat(repeats=sw, axis=0)
    y_repeated = y_weighted.repeat(repeats=sw)

    X_weighted, y_weighted, sw = shuffle(X_weighted, y_weighted, sw, random_state=0)

    # when the estimator has an internal CV scheme
    # we only use weights / repetitions in a specific CV group (here group=0)
    if "cv" in estimator_orig.get_params():
        groups_weighted = np.hstack(
            [np.full_like(y_weighted, 0), np.full_like(y, 1), np.full_like(y, 2)]
        )
        sw = np.hstack([sw, np.ones_like(y), np.ones_like(y)])
        X_weighted = np.vstack([X_weighted, X, X])
        y_weighted = np.hstack([y_weighted, y, y])
        splits_weighted = list(
            LeaveOneGroupOut().split(X_weighted, groups=groups_weighted)
        )
        estimator_weighted.set_params(cv=splits_weighted)

        groups_repeated = np.hstack(
            [np.full_like(y_repeated, 0), np.full_like(y, 1), np.full_like(y, 2)]
        )
        X_repeated = np.vstack([X_repeated, X, X])
        y_repeated = np.hstack([y_repeated, y, y])
        splits_repeated = list(
            LeaveOneGroupOut().split(X_repeated, groups=groups_repeated)
        )
        estimator_repeated.set_params(cv=splits_repeated)

    y_weighted = _enforce_estimator_tags_y(estimator_weighted, y_weighted)
    y_repeated = _enforce_estimator_tags_y(estimator_repeated, y_repeated)

    # convert to sparse X if needed
    if sparse_container is not None:
        X_weighted = sparse_container(X_weighted)
        X_repeated = sparse_container(X_repeated)

    estimator_repeated.fit(X_repeated, y=y_repeated, sample_weight=None)
    estimator_weighted.fit(X_weighted, y=y_weighted, sample_weight=sw)

    for method in ["predict_proba", "decision_function", "predict", "transform"]:
        if hasattr(estimator_orig, method):
            X_pred1 = getattr(estimator_repeated, method)(X)
            X_pred2 = getattr(estimator_weighted, method)(X)
            err_msg = (
                f"Comparing the output of {name}.{method} revealed that fitting "
                "with `sample_weight` is not equivalent to fitting with removed "
                "or repeated data points."
            )
            assert_allclose_dense_sparse(X_pred1, X_pred2, err_msg=err_msg)


def check_sample_weight_equivalence_on_dense_data(name, estimator_orig):
    _check_sample_weight_equivalence(name, estimator_orig, sparse_container=None)


def check_sample_weight_equivalence_on_sparse_data(name, estimator_orig):
    _check_sample_weight_equivalence(name, estimator_orig, sparse.csr_array)


def check_sample_weights_not_overwritten(name, estimator_orig):
    # check that estimators don't override the passed sample_weight parameter
    estimator = clone(estimator_orig)
    set_random_state(estimator, random_state=0)

    X = np.array(
        [
            [1, 3],
            [1, 3],
            [1, 3],
            [1, 3],
            [2, 1],
            [2, 1],
            [2, 1],
            [2, 1],
            [3, 3],
            [3, 3],
            [3, 3],
            [3, 3],
            [4, 1],
            [4, 1],
            [4, 1],
            [4, 1],
        ],
        dtype=np.float64,
    )
    y = np.array([1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2], dtype=int)
    y = _enforce_estimator_tags_y(estimator, y)

    sample_weight_original = np.ones(y.shape[0])
    sample_weight_original[0] = 10.0

    sample_weight_fit = sample_weight_original.copy()

    estimator.fit(X, y, sample_weight=sample_weight_fit)

    err_msg = f"{name} overwrote the original `sample_weight` given during fit"
    assert_allclose(sample_weight_fit, sample_weight_original, err_msg=err_msg)


@ignore_warnings(category=(FutureWarning, UserWarning))
def check_dtype_object(name, estimator_orig):
    # check that estimators treat dtype object as numeric if possible
    rng = np.random.RandomState(0)
    X = _enforce_estimator_tags_X(estimator_orig, rng.uniform(size=(40, 10)))
    X = X.astype(object)
    tags = get_tags(estimator_orig)
    y = (X[:, 0] * 4).astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    estimator.fit(X, y)
    if hasattr(estimator, "predict"):
        estimator.predict(X)

    if hasattr(estimator, "transform"):
        estimator.transform(X)

    err_msg = (
        "y with unknown label type is passed, but an error with no proper message "
        "is raised. You can use `type_of_target(..., raise_unknown=True)` to check "
        "and raise the right error, or include 'Unknown label type' in the error "
        "message."
    )
    with raises(Exception, match="Unknown label type", may_pass=True, err_msg=err_msg):
        estimator.fit(X, y.astype(object))

    if not tags.input_tags.string:
        X[0, 0] = {"foo": "bar"}
        # This error is raised by:
        # - `np.asarray` in `check_array`
        # - `_unique_python` for encoders
        msg = "argument must be .* string.* number"
        with raises(TypeError, match=msg):
            estimator.fit(X, y)
    else:
        # Estimators supporting string will not call np.asarray to convert the
        # data to numeric and therefore, the error will not be raised.
        # Checking for each element dtype in the input array will be costly.
        # Refer to #11401 for full discussion.
        estimator.fit(X, y)


def check_complex_data(name, estimator_orig):
    rng = np.random.RandomState(42)
    # check that estimators raise an exception on providing complex data
    X = rng.uniform(size=10) + 1j * rng.uniform(size=10)
    X = X.reshape(-1, 1)

    # Something both valid for classification and regression
    y = rng.randint(low=0, high=2, size=10) + 1j
    estimator = clone(estimator_orig)
    set_random_state(estimator, random_state=0)
    with raises(ValueError, match="Complex data not supported"):
        estimator.fit(X, y)


@ignore_warnings
def check_dict_unchanged(name, estimator_orig):
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)

    y = X[:, 0].astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)
    set_random_state(estimator, 1)

    estimator.fit(X, y)
    for method in ["predict", "transform", "decision_function", "predict_proba"]:
        if hasattr(estimator, method):
            dict_before = estimator.__dict__.copy()
            getattr(estimator, method)(X)
            assert estimator.__dict__ == dict_before, (
                "Estimator changes __dict__ during %s" % method
            )


def _is_public_parameter(attr):
    return not (attr.startswith("_") or attr.endswith("_"))


@ignore_warnings(category=FutureWarning)
def check_dont_overwrite_parameters(name, estimator_orig):
    # check that fit method only changes or sets private attributes
    if hasattr(estimator_orig.__init__, "deprecated_original"):
        # to not check deprecated classes
        return
    estimator = clone(estimator_orig)
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = X[:, 0].astype(int)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    dict_before_fit = estimator.__dict__.copy()
    estimator.fit(X, y)

    dict_after_fit = estimator.__dict__

    public_keys_after_fit = [
        key for key in dict_after_fit.keys() if _is_public_parameter(key)
    ]

    attrs_added_by_fit = [
        key for key in public_keys_after_fit if key not in dict_before_fit.keys()
    ]

    # check that fit doesn't add any public attribute
    assert not attrs_added_by_fit, (
        "Estimator adds public attribute(s) during"
        " the fit method."
        " Estimators are only allowed to add private attributes"
        " either started with _ or ended"
        " with _ but %s added" % ", ".join(attrs_added_by_fit)
    )

    # check that fit doesn't change any public attribute
    attrs_changed_by_fit = [
        key
        for key in public_keys_after_fit
        if (dict_before_fit[key] is not dict_after_fit[key])
    ]

    assert not attrs_changed_by_fit, (
        "Estimator changes public attribute(s) during"
        " the fit method. Estimators are only allowed"
        " to change attributes started"
        " or ended with _, but"
        " %s changed" % ", ".join(attrs_changed_by_fit)
    )


@ignore_warnings(category=FutureWarning)
def check_fit2d_predict1d(name, estimator_orig):
    # check by fitting a 2d array and predicting with a 1d array
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = X[:, 0].astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    estimator.fit(X, y)

    for method in ["predict", "transform", "decision_function", "predict_proba"]:
        if hasattr(estimator, method):
            with raises(ValueError, match="Reshape your data"):
                getattr(estimator, method)(X[0])


def _apply_on_subsets(func, X):
    # apply function on the whole set and on mini batches
    result_full = func(X)
    n_features = X.shape[1]
    result_by_batch = [func(batch.reshape(1, n_features)) for batch in X]

    # func can output tuple (e.g. score_samples)
    if isinstance(result_full, tuple):
        result_full = result_full[0]
        result_by_batch = list(map(lambda x: x[0], result_by_batch))

    if sparse.issparse(result_full):
        result_full = result_full.toarray()
        result_by_batch = [x.toarray() for x in result_by_batch]

    return np.ravel(result_full), np.ravel(result_by_batch)


@ignore_warnings(category=FutureWarning)
def check_methods_subset_invariance(name, estimator_orig):
    # check that method gives invariant results if applied
    # on mini batches or the whole set
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = X[:, 0].astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    estimator.fit(X, y)

    for method in [
        "predict",
        "transform",
        "decision_function",
        "score_samples",
        "predict_proba",
    ]:
        msg = ("{method} of {name} is not invariant when applied to a subset.").format(
            method=method, name=name
        )

        if hasattr(estimator, method):
            result_full, result_by_batch = _apply_on_subsets(
                getattr(estimator, method), X
            )
            assert_allclose(result_full, result_by_batch, atol=1e-7, err_msg=msg)


@ignore_warnings(category=FutureWarning)
def check_methods_sample_order_invariance(name, estimator_orig):
    # check that method gives invariant results if applied
    # on a subset with different sample order
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = X[:, 0].astype(np.int64)
    tags = get_tags(estimator_orig)
    if tags.classifier_tags is not None and not tags.classifier_tags.multi_class:
        y[y == 2] = 1
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 2

    set_random_state(estimator, 1)
    estimator.fit(X, y)

    idx = np.random.permutation(X.shape[0])

    for method in [
        "predict",
        "transform",
        "decision_function",
        "score_samples",
        "predict_proba",
    ]:
        msg = (
            "{method} of {name} is not invariant when applied to a dataset"
            "with different sample order."
        ).format(method=method, name=name)

        if hasattr(estimator, method):
            assert_allclose_dense_sparse(
                _safe_indexing(getattr(estimator, method)(X), idx),
                getattr(estimator, method)(_safe_indexing(X, idx)),
                atol=1e-9,
                err_msg=msg,
            )


@ignore_warnings
def check_fit2d_1sample(name, estimator_orig):
    # Check that fitting a 2d array with only one sample either works or
    # returns an informative message. The error message should either mention
    # the number of samples or the number of classes.
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(1, 10))
    X = _enforce_estimator_tags_X(estimator_orig, X)

    y = X[:, 0].astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)

    # min_cluster_size cannot be less than the data size for OPTICS.
    if name == "OPTICS":
        estimator.set_params(min_samples=1.0)

    # perplexity cannot be more than the number of samples for TSNE.
    if name == "TSNE":
        estimator.set_params(perplexity=0.5)

    msgs = [
        "1 sample",
        "n_samples = 1",
        "n_samples=1",
        "one sample",
        "1 class",
        "one class",
    ]

    with raises(ValueError, match=msgs, may_pass=True):
        estimator.fit(X, y)


@ignore_warnings
def check_fit2d_1feature(name, estimator_orig):
    # check fitting a 2d array with only 1 feature either works or returns
    # informative message
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(10, 1))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = X[:, 0].astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1
    # ensure two labels in subsample for RandomizedLogisticRegression
    if name == "RandomizedLogisticRegression":
        estimator.sample_fraction = 1
    # ensure non skipped trials for RANSACRegressor
    if name == "RANSACRegressor":
        estimator.residual_threshold = 0.5

    y = _enforce_estimator_tags_y(estimator, y)
    set_random_state(estimator, 1)

    msgs = [r"1 feature\(s\)", "n_features = 1", "n_features=1"]

    with raises(ValueError, match=msgs, may_pass=True):
        estimator.fit(X, y)


@ignore_warnings
def check_fit1d(name, estimator_orig):
    # check fitting 1d X array raises a ValueError
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20))
    y = X.astype(int)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    with raises(ValueError):
        estimator.fit(X, y)


@ignore_warnings(category=FutureWarning)
def check_transformer_general(name, transformer, readonly_memmap=False):
    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        n_features=2,
        cluster_std=0.1,
    )
    X = StandardScaler().fit_transform(X)
    X = _enforce_estimator_tags_X(transformer, X)

    if readonly_memmap:
        X, y = create_memmap_backed_data([X, y])

    _check_transformer(name, transformer, X, y)


@ignore_warnings(category=FutureWarning)
def check_transformer_data_not_an_array(name, transformer):
    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        n_features=2,
        cluster_std=0.1,
    )
    X = StandardScaler().fit_transform(X)
    X = _enforce_estimator_tags_X(transformer, X)
    this_X = _NotAnArray(X)
    this_y = _NotAnArray(np.asarray(y))
    _check_transformer(name, transformer, this_X, this_y)
    # try the same with some list
    _check_transformer(name, transformer, X.tolist(), y.tolist())


@ignore_warnings(category=FutureWarning)
def check_transformers_unfitted(name, transformer):
    X, y = _regression_dataset()

    transformer = clone(transformer)
    with raises(
        (AttributeError, ValueError),
        err_msg=(
            "The unfitted "
            f"transformer {name} does not raise an error when "
            "transform is called. Perhaps use "
            "check_is_fitted in transform."
        ),
    ):
        transformer.transform(X)


@ignore_warnings(category=FutureWarning)
def check_transformers_unfitted_stateless(name, transformer):
    """Check that using transform without prior fitting
    doesn't raise a NotFittedError for stateless transformers.
    """
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(20, 5))
    X = _enforce_estimator_tags_X(transformer, X)

    transformer = clone(transformer)
    X_trans = transformer.transform(X)

    assert X_trans.shape[0] == X.shape[0]


def _check_transformer(name, transformer_orig, X, y):
    n_samples, n_features = np.asarray(X).shape
    transformer = clone(transformer_orig)
    set_random_state(transformer)

    # fit

    if name in CROSS_DECOMPOSITION:
        y_ = np.c_[np.asarray(y), np.asarray(y)]
        y_[::2, 1] *= 2
        if isinstance(X, _NotAnArray):
            y_ = _NotAnArray(y_)
    else:
        y_ = y

    transformer.fit(X, y_)
    # fit_transform method should work on non fitted estimator
    transformer_clone = clone(transformer)
    X_pred = transformer_clone.fit_transform(X, y=y_)

    if isinstance(X_pred, tuple):
        for x_pred in X_pred:
            assert x_pred.shape[0] == n_samples
    else:
        # check for consistent n_samples
        assert X_pred.shape[0] == n_samples

    if hasattr(transformer, "transform"):
        if name in CROSS_DECOMPOSITION:
            X_pred2 = transformer.transform(X, y_)
            X_pred3 = transformer.fit_transform(X, y=y_)
        else:
            X_pred2 = transformer.transform(X)
            X_pred3 = transformer.fit_transform(X, y=y_)

        if get_tags(transformer_orig).non_deterministic:
            msg = name + " is non deterministic"
            raise SkipTest(msg)
        if isinstance(X_pred, tuple) and isinstance(X_pred2, tuple):
            for x_pred, x_pred2, x_pred3 in zip(X_pred, X_pred2, X_pred3):
                assert_allclose_dense_sparse(
                    x_pred,
                    x_pred2,
                    atol=1e-2,
                    err_msg="fit_transform and transform outcomes not consistent in %s"
                    % transformer,
                )
                assert_allclose_dense_sparse(
                    x_pred,
                    x_pred3,
                    atol=1e-2,
                    err_msg="consecutive fit_transform outcomes not consistent in %s"
                    % transformer,
                )
        else:
            assert_allclose_dense_sparse(
                X_pred,
                X_pred2,
                err_msg="fit_transform and transform outcomes not consistent in %s"
                % transformer,
                atol=1e-2,
            )
            assert_allclose_dense_sparse(
                X_pred,
                X_pred3,
                atol=1e-2,
                err_msg="consecutive fit_transform outcomes not consistent in %s"
                % transformer,
            )
            assert _num_samples(X_pred2) == n_samples
            assert _num_samples(X_pred3) == n_samples

        # raises error on malformed input for transform
        if (
            hasattr(X, "shape")
            and get_tags(transformer).requires_fit
            and X.ndim == 2
            and X.shape[1] > 1
        ):
            # If it's not an array, it does not have a 'T' property
            with raises(
                ValueError,
                err_msg=(
                    f"The transformer {name} does not raise an error "
                    "when the number of features in transform is different from "
                    "the number of features in fit."
                ),
            ):
                transformer.transform(X[:, :-1])


@ignore_warnings
def check_pipeline_consistency(name, estimator_orig):
    if get_tags(estimator_orig).non_deterministic:
        msg = name + " is non deterministic"
        raise SkipTest(msg)

    # check that make_pipeline(est) gives same score as est
    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        n_features=2,
        cluster_std=0.1,
    )
    X = _enforce_estimator_tags_X(estimator_orig, X, kernel=rbf_kernel)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)
    set_random_state(estimator)
    pipeline = make_pipeline(estimator)
    estimator.fit(X, y)
    pipeline.fit(X, y)

    funcs = ["score", "fit_transform"]

    for func_name in funcs:
        func = getattr(estimator, func_name, None)
        if func is not None:
            func_pipeline = getattr(pipeline, func_name)
            result = func(X, y)
            result_pipe = func_pipeline(X, y)
            assert_allclose_dense_sparse(result, result_pipe)


@ignore_warnings
def check_mixin_order(name, estimator_orig):
    """Check that mixins are inherited in the correct order."""
    # We define a list of edges, which in effect define a DAG of mixins and their
    # required order of inheritance.
    # This is of the form (mixin_a_should_be_before, mixin_b_should_be_after)
    dag = [
        (ClassifierMixin, BaseEstimator),
        (RegressorMixin, BaseEstimator),
        (ClusterMixin, BaseEstimator),
        (TransformerMixin, BaseEstimator),
        (BiclusterMixin, BaseEstimator),
        (OneToOneFeatureMixin, BaseEstimator),
        (ClassNamePrefixFeaturesOutMixin, BaseEstimator),
        (DensityMixin, BaseEstimator),
        (OutlierMixin, BaseEstimator),
        (MetaEstimatorMixin, BaseEstimator),
        (MultiOutputMixin, BaseEstimator),
    ]
    violations = []
    mro = type(estimator_orig).mro()
    for mixin_a, mixin_b in dag:
        if (
            mixin_a in mro
            and mixin_b in mro
            and mro.index(mixin_a) > mro.index(mixin_b)
        ):
            violations.append((mixin_a, mixin_b))
    violation_str = "\n".join(
        f"{mixin_a.__name__} comes before/left side of {mixin_b.__name__}"
        for mixin_a, mixin_b in violations
    )
    assert not violations, (
        f"{name} is inheriting from mixins in the wrong order. In general, in mixin "
        "inheritance, more specialized mixins must come before more general ones. "
        "This means, for instance, `BaseEstimator` should be on the right side of most "
        "other mixins. You need to change the order so that:\n"
        f"{violation_str}"
    )


@ignore_warnings
def check_fit_score_takes_y(name, estimator_orig):
    # check that all estimators accept an optional y
    # in fit and score so they can be used in pipelines
    rnd = np.random.RandomState(0)
    n_samples = 30
    X = rnd.uniform(size=(n_samples, 3))
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = np.arange(n_samples) % 3
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)
    set_random_state(estimator)

    funcs = ["fit", "score", "partial_fit", "fit_predict", "fit_transform"]
    for func_name in funcs:
        func = getattr(estimator, func_name, None)
        if func is not None:
            func(X, y)
            args = [p.name for p in signature(func).parameters.values()]
            if args[0] == "self":
                # available_if makes methods into functions
                # with an explicit "self", so need to shift arguments
                args = args[1:]
            assert args[1] in ["y", "Y"], (
                "Expected y or Y as second argument for method "
                "%s of %s. Got arguments: %r."
                % (func_name, type(estimator).__name__, args)
            )


@ignore_warnings
def check_estimators_dtypes(name, estimator_orig):
    rnd = np.random.RandomState(0)
    X_train_32 = 3 * rnd.uniform(size=(20, 5)).astype(np.float32)
    X_train_32 = _enforce_estimator_tags_X(estimator_orig, X_train_32)
    X_train_64 = X_train_32.astype(np.float64)
    X_train_int_64 = X_train_32.astype(np.int64)
    X_train_int_32 = X_train_32.astype(np.int32)
    y = np.array([1, 2] * 10, dtype=np.int64)
    y = _enforce_estimator_tags_y(estimator_orig, y)

    methods = ["predict", "transform", "decision_function", "predict_proba"]

    for X_train in [X_train_32, X_train_64, X_train_int_64, X_train_int_32]:
        estimator = clone(estimator_orig)
        set_random_state(estimator, 1)
        estimator.fit(X_train, y)

        for method in methods:
            if hasattr(estimator, method):
                getattr(estimator, method)(X_train)


def check_transformer_preserve_dtypes(name, transformer_orig):
    # check that dtype are preserved meaning if input X is of some dtype
    # X_transformed should be from the same dtype.
    transformer = clone(transformer_orig)
    if hasattr(transformer, "set_output"):
        transformer.set_output(transform="default")
    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        cluster_std=0.1,
    )
    X = StandardScaler().fit_transform(X)
    X = _enforce_estimator_tags_X(transformer_orig, X)

    for dtype in get_tags(transformer_orig).transformer_tags.preserves_dtype:
        X_cast = X.astype(dtype)
        set_random_state(transformer)
        X_trans1 = transformer.fit_transform(X_cast, y)
        X_trans2 = transformer.fit(X_cast, y).transform(X_cast)

        for Xt, method in zip([X_trans1, X_trans2], ["fit_transform", "transform"]):
            if isinstance(Xt, tuple):
                # cross-decompostion returns a tuple of (x_scores, y_scores)
                # when given y with fit_transform; only check the first element
                Xt = Xt[0]

            # check that the output dtype is preserved
            assert Xt.dtype == dtype, (
                f"{name} (method={method}) does not preserve dtype. "
                f"Original/Expected dtype={dtype}, got dtype={Xt.dtype}."
            )


@ignore_warnings(category=FutureWarning)
def check_estimators_empty_data_messages(name, estimator_orig):
    e = clone(estimator_orig)
    set_random_state(e, 1)

    X_zero_samples = np.empty(0).reshape(0, 3)
    # The precise message can change depending on whether X or y is
    # validated first. Let us test the type of exception only:
    err_msg = (
        f"The estimator {name} does not raise a ValueError when an "
        "empty data is used to train. Perhaps use check_array in train."
    )
    with raises(ValueError, err_msg=err_msg):
        e.fit(X_zero_samples, [])

    X_zero_features = np.empty(0).reshape(12, 0)
    # the following y should be accepted by both classifiers and regressors
    # and ignored by unsupervised models
    y = _enforce_estimator_tags_y(e, np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]))
    msg = r"0 feature\(s\) \(shape=\(\d*, 0\)\) while a minimum of \d* is required."
    with raises(ValueError, match=msg):
        e.fit(X_zero_features, y)


@ignore_warnings(category=FutureWarning)
def check_estimators_nan_inf(name, estimator_orig):
    # Checks that Estimator X's do not contain NaN or inf.
    rnd = np.random.RandomState(0)
    X_train_finite = _enforce_estimator_tags_X(
        estimator_orig, rnd.uniform(size=(10, 3))
    )
    X_train_nan = rnd.uniform(size=(10, 3))
    X_train_nan[0, 0] = np.nan
    X_train_inf = rnd.uniform(size=(10, 3))
    X_train_inf[0, 0] = np.inf
    y = np.ones(10)
    y[:5] = 0
    y = _enforce_estimator_tags_y(estimator_orig, y)
    error_string_fit = f"Estimator {name} doesn't check for NaN and inf in fit."
    error_string_predict = f"Estimator {name} doesn't check for NaN and inf in predict."
    error_string_transform = (
        f"Estimator {name} doesn't check for NaN and inf in transform."
    )
    for X_train in [X_train_nan, X_train_inf]:
        # catch deprecation warnings
        with ignore_warnings(category=FutureWarning):
            estimator = clone(estimator_orig)
            set_random_state(estimator, 1)
            # try to fit
            with raises(ValueError, match=["inf", "NaN"], err_msg=error_string_fit):
                estimator.fit(X_train, y)
            # actually fit
            estimator.fit(X_train_finite, y)

            # predict
            if hasattr(estimator, "predict"):
                with raises(
                    ValueError,
                    match=["inf", "NaN"],
                    err_msg=error_string_predict,
                ):
                    estimator.predict(X_train)

            # transform
            if hasattr(estimator, "transform"):
                with raises(
                    ValueError,
                    match=["inf", "NaN"],
                    err_msg=error_string_transform,
                ):
                    estimator.transform(X_train)


@ignore_warnings
def check_nonsquare_error(name, estimator_orig):
    """Test that error is thrown when non-square data provided."""

    X, y = make_blobs(n_samples=20, n_features=10)
    estimator = clone(estimator_orig)

    with raises(
        ValueError,
        err_msg=(
            f"The pairwise estimator {name} does not raise an error on non-square data"
        ),
    ):
        estimator.fit(X, y)


@ignore_warnings
def check_estimators_pickle(name, estimator_orig, readonly_memmap=False):
    """Test that we can pickle all estimators."""
    check_methods = ["predict", "transform", "decision_function", "predict_proba"]

    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        n_features=2,
        cluster_std=0.1,
    )

    X = _enforce_estimator_tags_X(estimator_orig, X, kernel=rbf_kernel)

    tags = get_tags(estimator_orig)
    # include NaN values when the estimator should deal with them
    if tags.input_tags.allow_nan:
        # set randomly 10 elements to np.nan
        rng = np.random.RandomState(42)
        mask = rng.choice(X.size, 10, replace=False)
        X.reshape(-1)[mask] = np.nan

    estimator = clone(estimator_orig)

    y = _enforce_estimator_tags_y(estimator, y)

    set_random_state(estimator)
    estimator.fit(X, y)

    if readonly_memmap:
        unpickled_estimator = create_memmap_backed_data(estimator)
    else:
        # No need to touch the file system in that case.
        pickled_estimator = pickle.dumps(estimator)
        module_name = estimator.__module__
        if module_name.startswith("sklearn.") and not (
            "test_" in module_name or module_name.endswith("_testing")
        ):
            # strict check for sklearn estimators that are not implemented in test
            # modules.
            assert b"_sklearn_version" in pickled_estimator
        unpickled_estimator = pickle.loads(pickled_estimator)

    result = dict()
    for method in check_methods:
        if hasattr(estimator, method):
            result[method] = getattr(estimator, method)(X)

    for method in result:
        unpickled_result = getattr(unpickled_estimator, method)(X)
        assert_allclose_dense_sparse(result[method], unpickled_result)


@ignore_warnings(category=FutureWarning)
def check_estimators_partial_fit_n_features(name, estimator_orig):
    # check if number of features changes between calls to partial_fit.
    if not hasattr(estimator_orig, "partial_fit"):
        return
    estimator = clone(estimator_orig)
    X, y = make_blobs(n_samples=50, random_state=1)
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = _enforce_estimator_tags_y(estimator_orig, y)

    try:
        if is_classifier(estimator):
            classes = np.unique(y)
            estimator.partial_fit(X, y, classes=classes)
        else:
            estimator.partial_fit(X, y)
    except NotImplementedError:
        return

    with raises(
        ValueError,
        err_msg=(
            f"The estimator {name} does not raise an error when the "
            "number of features changes between calls to partial_fit."
        ),
    ):
        estimator.partial_fit(X[:, :-1], y)


@ignore_warnings(category=FutureWarning)
def check_classifier_multioutput(name, estimator_orig):
    n_samples, n_labels, n_classes = 42, 5, 3
    tags = get_tags(estimator_orig)
    estimator = clone(estimator_orig)
    X, y = make_multilabel_classification(
        random_state=42, n_samples=n_samples, n_labels=n_labels, n_classes=n_classes
    )
    X = _enforce_estimator_tags_X(estimator, X)
    estimator.fit(X, y)
    y_pred = estimator.predict(X)

    assert y_pred.shape == (n_samples, n_classes), (
        "The shape of the prediction for multioutput data is "
        "incorrect. Expected {}, got {}.".format((n_samples, n_labels), y_pred.shape)
    )
    assert y_pred.dtype.kind == "i"

    if hasattr(estimator, "decision_function"):
        decision = estimator.decision_function(X)
        assert isinstance(decision, np.ndarray)
        assert decision.shape == (n_samples, n_classes), (
            "The shape of the decision function output for "
            "multioutput data is incorrect. Expected {}, got {}.".format(
                (n_samples, n_classes), decision.shape
            )
        )

        dec_pred = (decision > 0).astype(int)
        dec_exp = estimator.classes_[dec_pred]
        assert_array_equal(dec_exp, y_pred)

    if hasattr(estimator, "predict_proba"):
        y_prob = estimator.predict_proba(X)

        if isinstance(y_prob, list) and not tags.classifier_tags.poor_score:
            for i in range(n_classes):
                assert y_prob[i].shape == (n_samples, 2), (
                    "The shape of the probability for multioutput data is"
                    " incorrect. Expected {}, got {}.".format(
                        (n_samples, 2), y_prob[i].shape
                    )
                )
                assert_array_equal(
                    np.argmax(y_prob[i], axis=1).astype(int), y_pred[:, i]
                )
        elif not tags.classifier_tags.poor_score:
            assert y_prob.shape == (n_samples, n_classes), (
                "The shape of the probability for multioutput data is"
                " incorrect. Expected {}, got {}.".format(
                    (n_samples, n_classes), y_prob.shape
                )
            )
            assert_array_equal(y_prob.round().astype(int), y_pred)

    if hasattr(estimator, "decision_function") and hasattr(estimator, "predict_proba"):
        for i in range(n_classes):
            y_proba = estimator.predict_proba(X)[:, i]
            y_decision = estimator.decision_function(X)
            assert_array_equal(rankdata(y_proba), rankdata(y_decision[:, i]))


@ignore_warnings(category=FutureWarning)
def check_regressor_multioutput(name, estimator):
    estimator = clone(estimator)
    n_samples = n_features = 10

    if not _is_pairwise_metric(estimator):
        n_samples = n_samples + 1

    X, y = make_regression(
        random_state=42, n_targets=5, n_samples=n_samples, n_features=n_features
    )
    X = _enforce_estimator_tags_X(estimator, X)

    estimator.fit(X, y)
    y_pred = estimator.predict(X)

    assert y_pred.dtype == np.dtype("float64"), (
        "Multioutput predictions by a regressor are expected to be"
        f" floating-point precision. Got {y_pred.dtype} instead"
    )
    assert y_pred.shape == y.shape, (
        "The shape of the prediction for multioutput data is incorrect."
        f" Expected {y_pred.shape}, got {y.shape}."
    )


@ignore_warnings(category=FutureWarning)
def check_clustering(name, clusterer_orig, readonly_memmap=False):
    clusterer = clone(clusterer_orig)
    X, y = make_blobs(n_samples=50, random_state=1)
    X, y = shuffle(X, y, random_state=7)
    X = StandardScaler().fit_transform(X)
    rng = np.random.RandomState(7)
    X_noise = np.concatenate([X, rng.uniform(low=-3, high=3, size=(5, 2))])

    if readonly_memmap:
        X, y, X_noise = create_memmap_backed_data([X, y, X_noise])

    n_samples, n_features = X.shape
    # catch deprecation and neighbors warnings
    if hasattr(clusterer, "n_clusters"):
        clusterer.set_params(n_clusters=3)
    set_random_state(clusterer)
    if name == "AffinityPropagation":
        clusterer.set_params(preference=-100)
        clusterer.set_params(max_iter=100)

    # fit
    clusterer.fit(X)
    # with lists
    clusterer.fit(X.tolist())

    pred = clusterer.labels_
    assert pred.shape == (n_samples,)
    assert adjusted_rand_score(pred, y) > 0.4
    if get_tags(clusterer).non_deterministic:
        return
    set_random_state(clusterer)
    with warnings.catch_warnings(record=True):
        pred2 = clusterer.fit_predict(X)
    assert_array_equal(pred, pred2)

    # fit_predict(X) and labels_ should be of type int
    assert pred.dtype in [np.dtype("int32"), np.dtype("int64")]
    assert pred2.dtype in [np.dtype("int32"), np.dtype("int64")]

    # Add noise to X to test the possible values of the labels
    labels = clusterer.fit_predict(X_noise)

    # There should be at least one sample in every cluster. Equivalently
    # labels_ should contain all the consecutive values between its
    # min and its max.
    labels_sorted = np.unique(labels)
    assert_array_equal(
        labels_sorted, np.arange(labels_sorted[0], labels_sorted[-1] + 1)
    )

    # Labels are expected to start at 0 (no noise) or -1 (if noise)
    assert labels_sorted[0] in [0, -1]
    # Labels should be less than n_clusters - 1
    if hasattr(clusterer, "n_clusters"):
        n_clusters = getattr(clusterer, "n_clusters")
        assert n_clusters - 1 >= labels_sorted[-1]
    # else labels should be less than max(labels_) which is necessarily true


@ignore_warnings(category=FutureWarning)
def check_clusterer_compute_labels_predict(name, clusterer_orig):
    """Check that predict is invariant of compute_labels."""
    X, y = make_blobs(n_samples=20, random_state=0)
    clusterer = clone(clusterer_orig)
    set_random_state(clusterer)

    if hasattr(clusterer, "compute_labels"):
        # MiniBatchKMeans
        X_pred1 = clusterer.fit(X).predict(X)
        clusterer.set_params(compute_labels=False)
        X_pred2 = clusterer.fit(X).predict(X)
        assert_array_equal(X_pred1, X_pred2)


@ignore_warnings(category=FutureWarning)
def check_classifiers_one_label(name, classifier_orig):
    error_string_fit = "Classifier can't train when only one class is present."
    error_string_predict = "Classifier can't predict when only one class is present."
    classifier = clone(classifier_orig)
    rnd = np.random.RandomState(0)
    X_train = rnd.uniform(size=(10, 3))
    X_test = rnd.uniform(size=(10, 3))
    X_train, X_test = _enforce_estimator_tags_X(classifier, X_train, X_test=X_test)
    y = np.ones(10)
    # catch deprecation warnings
    with ignore_warnings(category=FutureWarning):
        with raises(
            ValueError, match="class", may_pass=True, err_msg=error_string_fit
        ) as cm:
            classifier.fit(X_train, y)

        if cm.raised_and_matched:
            # ValueError was raised with proper error message
            return

        assert_array_equal(classifier.predict(X_test), y, err_msg=error_string_predict)


@ignore_warnings(category=FutureWarning)
def check_classifiers_one_label_sample_weights(name, classifier_orig):
    """Check that classifiers accepting sample_weight fit or throws a ValueError with
    an explicit message if the problem is reduced to one class.
    """
    error_fit = (
        f"{name} failed when fitted on one label after sample_weight trimming. Error "
        "message is not explicit, it should have 'class'."
    )
    error_predict = f"{name} prediction results should only output the remaining class."
    rnd = np.random.RandomState(0)
    # X should be square for test on SVC with precomputed kernel
    X_train = rnd.uniform(size=(10, 10))
    X_test = rnd.uniform(size=(10, 10))
    y = np.arange(10) % 2
    sample_weight = y.copy()  # select a single class
    classifier = clone(classifier_orig)

    if has_fit_parameter(classifier, "sample_weight"):
        match = [r"\bclass(es)?\b", error_predict]
        err_type, err_msg = (AssertionError, ValueError), error_fit
    else:
        match = r"\bsample_weight\b"
        err_type, err_msg = (TypeError, ValueError), None

    with raises(err_type, match=match, may_pass=True, err_msg=err_msg) as cm:
        classifier.fit(X_train, y, sample_weight=sample_weight)
        if cm.raised_and_matched:
            # raise the proper error type with the proper error message
            return
        # for estimators that do not fail, they should be able to predict the only
        # class remaining during fit
        assert_array_equal(
            classifier.predict(X_test), np.ones(10), err_msg=error_predict
        )


@ignore_warnings  # Warnings are raised by decision function
def check_classifiers_train(
    name, classifier_orig, readonly_memmap=False, X_dtype="float64"
):
    X_m, y_m = make_blobs(n_samples=300, random_state=0)
    X_m = X_m.astype(X_dtype)
    X_m, y_m = shuffle(X_m, y_m, random_state=7)
    X_m = StandardScaler().fit_transform(X_m)
    # generate binary problem from multi-class one
    y_b = y_m[y_m != 2]
    X_b = X_m[y_m != 2]

    if name in ["BernoulliNB", "MultinomialNB", "ComplementNB", "CategoricalNB"]:
        X_m -= X_m.min()
        X_b -= X_b.min()

    if readonly_memmap:
        X_m, y_m, X_b, y_b = create_memmap_backed_data([X_m, y_m, X_b, y_b])

    problems = [(X_b, y_b)]
    tags = get_tags(classifier_orig)
    if tags.classifier_tags.multi_class:
        problems.append((X_m, y_m))

    for X, y in problems:
        classes = np.unique(y)
        n_classes = len(classes)
        n_samples, n_features = X.shape
        classifier = clone(classifier_orig)
        X = _enforce_estimator_tags_X(classifier, X)
        y = _enforce_estimator_tags_y(classifier, y)

        set_random_state(classifier)
        # raises error on malformed input for fit
        if not tags.no_validation:
            with raises(
                ValueError,
                err_msg=(
                    f"The classifier {name} does not raise an error when "
                    "incorrect/malformed input data for fit is passed. The number "
                    "of training examples is not the same as the number of "
                    "labels. Perhaps use check_X_y in fit."
                ),
            ):
                classifier.fit(X, y[:-1])

        # fit
        classifier.fit(X, y)
        # with lists
        classifier.fit(X.tolist(), y.tolist())
        assert hasattr(classifier, "classes_")
        y_pred = classifier.predict(X)

        assert y_pred.shape == (n_samples,)
        # training set performance
        if not tags.classifier_tags.poor_score:
            assert accuracy_score(y, y_pred) > 0.83

        # raises error on malformed input for predict
        msg_pairwise = (
            "The classifier {} does not raise an error when shape of X in "
            " {} is not equal to (n_test_samples, n_training_samples)"
        )
        msg = (
            "The classifier {} does not raise an error when the number of "
            "features in {} is different from the number of features in "
            "fit."
        )

        if not tags.no_validation:
            if tags.input_tags.pairwise:
                with raises(
                    ValueError,
                    err_msg=msg_pairwise.format(name, "predict"),
                ):
                    classifier.predict(X.reshape(-1, 1))
            else:
                with raises(ValueError, err_msg=msg.format(name, "predict")):
                    classifier.predict(X.T)
        if hasattr(classifier, "decision_function"):
            try:
                # decision_function agrees with predict
                decision = classifier.decision_function(X)
                if n_classes == 2:
                    if tags.target_tags.single_output:
                        assert decision.shape == (n_samples,)
                    else:
                        assert decision.shape == (n_samples, 1)
                    dec_pred = (decision.ravel() > 0).astype(int)
                    assert_array_equal(dec_pred, y_pred)
                else:
                    assert decision.shape == (n_samples, n_classes)
                    assert_array_equal(np.argmax(decision, axis=1), y_pred)

                # raises error on malformed input for decision_function
                if not tags.no_validation:
                    if tags.input_tags.pairwise:
                        with raises(
                            ValueError,
                            err_msg=msg_pairwise.format(name, "decision_function"),
                        ):
                            classifier.decision_function(X.reshape(-1, 1))
                    else:
                        with raises(
                            ValueError,
                            err_msg=msg.format(name, "decision_function"),
                        ):
                            classifier.decision_function(X.T)
            except NotImplementedError:
                pass

        if hasattr(classifier, "predict_proba"):
            # predict_proba agrees with predict
            y_prob = classifier.predict_proba(X)
            assert y_prob.shape == (n_samples, n_classes)
            assert_array_equal(np.argmax(y_prob, axis=1), y_pred)
            # check that probas for all classes sum to one
            assert_array_almost_equal(np.sum(y_prob, axis=1), np.ones(n_samples))
            if not tags.no_validation:
                # raises error on malformed input for predict_proba
                if tags.input_tags.pairwise:
                    with raises(
                        ValueError,
                        err_msg=msg_pairwise.format(name, "predict_proba"),
                    ):
                        classifier.predict_proba(X.reshape(-1, 1))
                else:
                    with raises(
                        ValueError,
                        err_msg=msg.format(name, "predict_proba"),
                    ):
                        classifier.predict_proba(X.T)
            if hasattr(classifier, "predict_log_proba"):
                # predict_log_proba is a transformation of predict_proba
                y_log_prob = classifier.predict_log_proba(X)
                assert_allclose(y_log_prob, np.log(y_prob), 8, atol=1e-9)
                assert_array_equal(np.argsort(y_log_prob), np.argsort(y_prob))


def check_outlier_corruption(num_outliers, expected_outliers, decision):
    # Check for deviation from the precise given contamination level that may
    # be due to ties in the anomaly scores.
    if num_outliers < expected_outliers:
        start = num_outliers
        end = expected_outliers + 1
    else:
        start = expected_outliers
        end = num_outliers + 1

    # ensure that all values in the 'critical area' are tied,
    # leading to the observed discrepancy between provided
    # and actual contamination levels.
    sorted_decision = np.sort(decision)
    msg = (
        "The number of predicted outliers is not equal to the expected "
        "number of outliers and this difference is not explained by the "
        "number of ties in the decision_function values"
    )
    assert len(np.unique(sorted_decision[start:end])) == 1, msg


def check_outliers_train(name, estimator_orig, readonly_memmap=True):
    n_samples = 300
    X, _ = make_blobs(n_samples=n_samples, random_state=0)
    X = shuffle(X, random_state=7)

    if readonly_memmap:
        X = create_memmap_backed_data(X)

    n_samples, n_features = X.shape
    estimator = clone(estimator_orig)
    set_random_state(estimator)

    # fit
    estimator.fit(X)
    # with lists
    estimator.fit(X.tolist())

    y_pred = estimator.predict(X)
    assert y_pred.shape == (n_samples,)
    assert y_pred.dtype.kind == "i"
    assert_array_equal(np.unique(y_pred), np.array([-1, 1]))

    decision = estimator.decision_function(X)
    scores = estimator.score_samples(X)
    for output in [decision, scores]:
        assert output.dtype == np.dtype("float")
        assert output.shape == (n_samples,)

    # raises error on malformed input for predict
    with raises(ValueError):
        estimator.predict(X.T)

    # decision_function agrees with predict
    dec_pred = (decision >= 0).astype(int)
    dec_pred[dec_pred == 0] = -1
    assert_array_equal(dec_pred, y_pred)

    # raises error on malformed input for decision_function
    with raises(ValueError):
        estimator.decision_function(X.T)

    # decision_function is a translation of score_samples
    y_dec = scores - estimator.offset_
    assert_allclose(y_dec, decision)

    # raises error on malformed input for score_samples
    with raises(ValueError):
        estimator.score_samples(X.T)

    # contamination parameter (not for OneClassSVM which has the nu parameter)
    if hasattr(estimator, "contamination") and not hasattr(estimator, "novelty"):
        # proportion of outliers equal to contamination parameter when not
        # set to 'auto'. This is true for the training set and cannot thus be
        # checked as follows for estimators with a novelty parameter such as
        # LocalOutlierFactor (tested in check_outliers_fit_predict)
        expected_outliers = 30
        contamination = expected_outliers / n_samples
        estimator.set_params(contamination=contamination)
        estimator.fit(X)
        y_pred = estimator.predict(X)

        num_outliers = np.sum(y_pred != 1)
        # num_outliers should be equal to expected_outliers unless
        # there are ties in the decision_function values. this can
        # only be tested for estimators with a decision_function
        # method, i.e. all estimators except LOF which is already
        # excluded from this if branch.
        if num_outliers != expected_outliers:
            decision = estimator.decision_function(X)
            check_outlier_corruption(num_outliers, expected_outliers, decision)


def check_outlier_contamination(name, estimator_orig):
    # Check that the contamination parameter is in (0.0, 0.5] when it is an
    # interval constraint.

    if not hasattr(estimator_orig, "_parameter_constraints"):
        # Only estimator implementing parameter constraints will be checked
        return

    if "contamination" not in estimator_orig._parameter_constraints:
        return

    contamination_constraints = estimator_orig._parameter_constraints["contamination"]
    if not any([isinstance(c, Interval) for c in contamination_constraints]):
        raise AssertionError(
            "contamination constraints should contain a Real Interval constraint."
        )

    for constraint in contamination_constraints:
        if isinstance(constraint, Interval):
            assert (
                constraint.type == Real
                and constraint.left >= 0.0
                and constraint.right <= 0.5
                and (constraint.left > 0 or constraint.closed in {"right", "neither"})
            ), "contamination constraint should be an interval in (0, 0.5]"


@ignore_warnings(category=FutureWarning)
def check_classifiers_multilabel_representation_invariance(name, classifier_orig):
    X, y = make_multilabel_classification(
        n_samples=100,
        n_features=2,
        n_classes=5,
        n_labels=3,
        length=50,
        allow_unlabeled=True,
        random_state=0,
    )
    X = scale(X)

    X_train, y_train = X[:80], y[:80]
    X_test = X[80:]
    X_train, X_test = _enforce_estimator_tags_X(classifier_orig, X_train, X_test=X_test)

    y_train_list_of_lists = y_train.tolist()
    y_train_list_of_arrays = list(y_train)

    classifier = clone(classifier_orig)
    set_random_state(classifier)

    y_pred = classifier.fit(X_train, y_train).predict(X_test)

    y_pred_list_of_lists = classifier.fit(X_train, y_train_list_of_lists).predict(
        X_test
    )

    y_pred_list_of_arrays = classifier.fit(X_train, y_train_list_of_arrays).predict(
        X_test
    )

    assert_array_equal(y_pred, y_pred_list_of_arrays)
    assert_array_equal(y_pred, y_pred_list_of_lists)

    assert y_pred.dtype == y_pred_list_of_arrays.dtype
    assert y_pred.dtype == y_pred_list_of_lists.dtype
    assert type(y_pred) == type(y_pred_list_of_arrays)
    assert type(y_pred) == type(y_pred_list_of_lists)


@ignore_warnings(category=FutureWarning)
def check_classifiers_multilabel_output_format_predict(name, classifier_orig):
    """Check the output of the `predict` method for classifiers supporting
    multilabel-indicator targets."""
    classifier = clone(classifier_orig)
    set_random_state(classifier)

    n_samples, test_size, n_outputs = 100, 25, 5
    X, y = make_multilabel_classification(
        n_samples=n_samples,
        n_features=2,
        n_classes=n_outputs,
        n_labels=3,
        length=50,
        allow_unlabeled=True,
        random_state=0,
    )
    X = scale(X)

    X_train, X_test = X[:-test_size], X[-test_size:]
    y_train, y_test = y[:-test_size], y[-test_size:]
    X_train, X_test = _enforce_estimator_tags_X(classifier_orig, X_train, X_test=X_test)
    classifier.fit(X_train, y_train)

    response_method_name = "predict"
    predict_method = getattr(classifier, response_method_name, None)
    if predict_method is None:
        raise SkipTest(f"{name} does not have a {response_method_name} method.")

    y_pred = predict_method(X_test)

    # y_pred.shape -> y_test.shape with the same dtype
    assert isinstance(y_pred, np.ndarray), (
        f"{name}.predict is expected to output a NumPy array. Got "
        f"{type(y_pred)} instead."
    )
    assert y_pred.shape == y_test.shape, (
        f"{name}.predict outputs a NumPy array of shape {y_pred.shape} "
        f"instead of {y_test.shape}."
    )
    assert y_pred.dtype == y_test.dtype, (
        f"{name}.predict does not output the same dtype than the targets. "
        f"Got {y_pred.dtype} instead of {y_test.dtype}."
    )


@ignore_warnings(category=FutureWarning)
def check_classifiers_multilabel_output_format_predict_proba(name, classifier_orig):
    """Check the output of the `predict_proba` method for classifiers supporting
    multilabel-indicator targets."""
    classifier = clone(classifier_orig)
    set_random_state(classifier)

    n_samples, test_size, n_outputs = 100, 25, 5
    X, y = make_multilabel_classification(
        n_samples=n_samples,
        n_features=2,
        n_classes=n_outputs,
        n_labels=3,
        length=50,
        allow_unlabeled=True,
        random_state=0,
    )
    X = scale(X)

    X_train, X_test = X[:-test_size], X[-test_size:]
    y_train = y[:-test_size]
    X_train, X_test = _enforce_estimator_tags_X(classifier_orig, X_train, X_test=X_test)
    classifier.fit(X_train, y_train)

    response_method_name = "predict_proba"
    predict_proba_method = getattr(classifier, response_method_name, None)
    if predict_proba_method is None:
        raise SkipTest(f"{name} does not have a {response_method_name} method.")

    y_pred = predict_proba_method(X_test)

    # y_pred.shape -> 2 possibilities:
    # - list of length n_outputs of shape (n_samples, 2);
    # - ndarray of shape (n_samples, n_outputs).
    # dtype should be floating
    if isinstance(y_pred, list):
        assert len(y_pred) == n_outputs, (
            f"When {name}.predict_proba returns a list, the list should "
            "be of length n_outputs and contain NumPy arrays. Got length "
            f"of {len(y_pred)} instead of {n_outputs}."
        )
        for pred in y_pred:
            assert pred.shape == (test_size, 2), (
                f"When {name}.predict_proba returns a list, this list "
                "should contain NumPy arrays of shape (n_samples, 2). Got "
                f"NumPy arrays of shape {pred.shape} instead of "
                f"{(test_size, 2)}."
            )
            assert pred.dtype.kind == "f", (
                f"When {name}.predict_proba returns a list, it should "
                "contain NumPy arrays with floating dtype. Got "
                f"{pred.dtype} instead."
            )
            # check that we have the correct probabilities
            err_msg = (
                f"When {name}.predict_proba returns a list, each NumPy "
                "array should contain probabilities for each class and "
                "thus each row should sum to 1 (or close to 1 due to "
                "numerical errors)."
            )
            assert_allclose(pred.sum(axis=1), 1, err_msg=err_msg)
    elif isinstance(y_pred, np.ndarray):
        assert y_pred.shape == (test_size, n_outputs), (
            f"When {name}.predict_proba returns a NumPy array, the "
            f"expected shape is (n_samples, n_outputs). Got {y_pred.shape}"
            f" instead of {(test_size, n_outputs)}."
        )
        assert y_pred.dtype.kind == "f", (
            f"When {name}.predict_proba returns a NumPy array, the "
            f"expected data type is floating. Got {y_pred.dtype} instead."
        )
        err_msg = (
            f"When {name}.predict_proba returns a NumPy array, this array "
            "is expected to provide probabilities of the positive class "
            "and should therefore contain values between 0 and 1."
        )
        assert_array_less(0, y_pred, err_msg=err_msg)
        assert_array_less(y_pred, 1, err_msg=err_msg)
    else:
        raise ValueError(
            f"Unknown returned type {type(y_pred)} by {name}."
            "predict_proba. A list or a Numpy array is expected."
        )


@ignore_warnings(category=FutureWarning)
def check_classifiers_multilabel_output_format_decision_function(name, classifier_orig):
    """Check the output of the `decision_function` method for classifiers supporting
    multilabel-indicator targets."""
    classifier = clone(classifier_orig)
    set_random_state(classifier)

    n_samples, test_size, n_outputs = 100, 25, 5
    X, y = make_multilabel_classification(
        n_samples=n_samples,
        n_features=2,
        n_classes=n_outputs,
        n_labels=3,
        length=50,
        allow_unlabeled=True,
        random_state=0,
    )
    X = scale(X)

    X_train, X_test = X[:-test_size], X[-test_size:]
    y_train = y[:-test_size]
    X_train, X_test = _enforce_estimator_tags_X(classifier_orig, X_train, X_test=X_test)
    classifier.fit(X_train, y_train)

    response_method_name = "decision_function"
    decision_function_method = getattr(classifier, response_method_name, None)
    if decision_function_method is None:
        raise SkipTest(f"{name} does not have a {response_method_name} method.")

    y_pred = decision_function_method(X_test)

    # y_pred.shape -> y_test.shape with floating dtype
    assert isinstance(y_pred, np.ndarray), (
        f"{name}.decision_function is expected to output a NumPy array."
        f" Got {type(y_pred)} instead."
    )
    assert y_pred.shape == (test_size, n_outputs), (
        f"{name}.decision_function is expected to provide a NumPy array "
        f"of shape (n_samples, n_outputs). Got {y_pred.shape} instead of "
        f"{(test_size, n_outputs)}."
    )
    assert y_pred.dtype.kind == "f", (
        f"{name}.decision_function is expected to output a floating dtype."
        f" Got {y_pred.dtype} instead."
    )


@ignore_warnings(category=FutureWarning)
def check_get_feature_names_out_error(name, estimator_orig):
    """Check the error raised by get_feature_names_out when called before fit.

    Unfitted estimators with get_feature_names_out should raise a NotFittedError.
    """

    estimator = clone(estimator_orig)
    err_msg = (
        f"Estimator {name} should have raised a NotFitted error when fit is called"
        " before get_feature_names_out"
    )
    with raises(NotFittedError, err_msg=err_msg):
        estimator.get_feature_names_out()


@ignore_warnings(category=FutureWarning)
def check_estimators_fit_returns_self(name, estimator_orig):
    """Check if self is returned when calling fit."""
    X, y = make_blobs(random_state=0, n_samples=21)
    X = _enforce_estimator_tags_X(estimator_orig, X)

    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    set_random_state(estimator)
    assert estimator.fit(X, y) is estimator


@ignore_warnings(category=FutureWarning)
def check_readonly_memmap_input(name, estimator_orig):
    """Check that the estimator can handle readonly memmap backed data.

    This is particularly needed to support joblib parallelisation.
    """
    X, y = make_blobs(random_state=0, n_samples=21)
    X = _enforce_estimator_tags_X(estimator_orig, X)

    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    X, y = create_memmap_backed_data([X, y])

    set_random_state(estimator)
    # This should not raise an error and should return self
    assert estimator.fit(X, y) is estimator


@ignore_warnings
def check_estimators_unfitted(name, estimator_orig):
    """Check that predict raises an exception in an unfitted estimator.

    Unfitted estimators should raise a NotFittedError.
    """
    err_msg = (
        "Estimator should raise a NotFittedError when calling `{method}` before fit. "
        "Either call `check_is_fitted(self)` at the beginning of `{method}` or "
        "set `tags.requires_fit=False` on estimator tags to disable this check.\n"
        "- `check_is_fitted`: https://scikit-learn.org/dev/modules/generated/sklearn."
        "utils.validation.check_is_fitted.html\n"
        "- Estimator Tags: https://scikit-learn.org/dev/developers/develop."
        "html#estimator-tags"
    )
    # Common test for Regressors, Classifiers and Outlier detection estimators
    X, y = _regression_dataset()

    estimator = clone(estimator_orig)
    for method in (
        "decision_function",
        "predict",
        "predict_proba",
        "predict_log_proba",
    ):
        if hasattr(estimator, method):
            with raises(NotFittedError, err_msg=err_msg.format(method=method)):
                getattr(estimator, method)(X)


@ignore_warnings(category=FutureWarning)
def check_supervised_y_2d(name, estimator_orig):
    tags = get_tags(estimator_orig)
    rnd = np.random.RandomState(0)
    n_samples = 30
    X = _enforce_estimator_tags_X(estimator_orig, rnd.uniform(size=(n_samples, 3)))
    y = np.arange(n_samples) % 3
    y = _enforce_estimator_tags_y(estimator_orig, y)
    estimator = clone(estimator_orig)
    set_random_state(estimator)
    # fit
    estimator.fit(X, y)
    y_pred = estimator.predict(X)

    set_random_state(estimator)
    # Check that when a 2D y is given, a DataConversionWarning is
    # raised
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", DataConversionWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        estimator.fit(X, y[:, np.newaxis])
    y_pred_2d = estimator.predict(X)
    msg = "expected 1 DataConversionWarning, got: %s" % ", ".join(
        [str(w_x) for w_x in w]
    )
    if not tags.target_tags.multi_output:
        # check that we warned if we don't support multi-output
        assert len(w) > 0, msg
        assert (
            "DataConversionWarning('A column-vector y"
            " was passed when a 1d array was expected" in msg
        )
    assert_allclose(y_pred.ravel(), y_pred_2d.ravel())


@ignore_warnings
def check_classifiers_predictions(X, y, name, classifier_orig):
    classes = np.unique(y)
    classifier = clone(classifier_orig)
    if name == "BernoulliNB":
        X = X > X.mean()
    set_random_state(classifier)

    classifier.fit(X, y)
    y_pred = classifier.predict(X)

    if hasattr(classifier, "decision_function"):
        decision = classifier.decision_function(X)
        assert isinstance(decision, np.ndarray)
        if len(classes) == 2:
            dec_pred = (decision.ravel() > 0).astype(int)
            dec_exp = classifier.classes_[dec_pred]
            assert_array_equal(
                dec_exp,
                y_pred,
                err_msg=(
                    "decision_function does not match "
                    "classifier for %r: expected '%s', got '%s'"
                )
                % (
                    classifier,
                    ", ".join(map(str, dec_exp)),
                    ", ".join(map(str, y_pred)),
                ),
            )
        elif getattr(classifier, "decision_function_shape", "ovr") == "ovr":
            decision_y = np.argmax(decision, axis=1).astype(int)
            y_exp = classifier.classes_[decision_y]
            assert_array_equal(
                y_exp,
                y_pred,
                err_msg=(
                    "decision_function does not match "
                    "classifier for %r: expected '%s', got '%s'"
                )
                % (
                    classifier,
                    ", ".join(map(str, y_exp)),
                    ", ".join(map(str, y_pred)),
                ),
            )

    assert_array_equal(
        classes,
        classifier.classes_,
        err_msg="Unexpected classes_ attribute for %r: expected '%s', got '%s'"
        % (
            classifier,
            ", ".join(map(str, classes)),
            ", ".join(map(str, classifier.classes_)),
        ),
    )


def _choose_check_classifiers_labels(name, y, y_names):
    # Semisupervised classifiers use -1 as the indicator for an unlabeled
    # sample.
    return (
        y
        if name in ["LabelPropagation", "LabelSpreading", "SelfTrainingClassifier"]
        else y_names
    )


def check_classifiers_classes(name, classifier_orig):
    X_multiclass, y_multiclass = make_blobs(
        n_samples=30, random_state=0, cluster_std=0.1
    )
    X_multiclass, y_multiclass = shuffle(X_multiclass, y_multiclass, random_state=7)
    X_multiclass = StandardScaler().fit_transform(X_multiclass)

    X_binary = X_multiclass[y_multiclass != 2]
    y_binary = y_multiclass[y_multiclass != 2]

    X_multiclass = _enforce_estimator_tags_X(classifier_orig, X_multiclass)
    X_binary = _enforce_estimator_tags_X(classifier_orig, X_binary)

    labels_multiclass = ["one", "two", "three"]
    labels_binary = ["one", "two"]

    y_names_multiclass = np.take(labels_multiclass, y_multiclass)
    y_names_binary = np.take(labels_binary, y_binary)

    problems = [(X_binary, y_binary, y_names_binary)]
    if get_tags(classifier_orig).classifier_tags.multi_class:
        problems.append((X_multiclass, y_multiclass, y_names_multiclass))

    for X, y, y_names in problems:
        for y_names_i in [y_names, y_names.astype("O")]:
            y_ = _choose_check_classifiers_labels(name, y, y_names_i)
            check_classifiers_predictions(X, y_, name, classifier_orig)

    labels_binary = [-1, 1]
    y_names_binary = np.take(labels_binary, y_binary)
    y_binary = _choose_check_classifiers_labels(name, y_binary, y_names_binary)
    check_classifiers_predictions(X_binary, y_binary, name, classifier_orig)


@ignore_warnings(category=FutureWarning)
def check_regressors_int(name, regressor_orig):
    X, _ = _regression_dataset()
    X = _enforce_estimator_tags_X(regressor_orig, X[:50])
    rnd = np.random.RandomState(0)
    y = rnd.randint(3, size=X.shape[0])
    y = _enforce_estimator_tags_y(regressor_orig, y)
    rnd = np.random.RandomState(0)
    # separate estimators to control random seeds
    regressor_1 = clone(regressor_orig)
    regressor_2 = clone(regressor_orig)
    set_random_state(regressor_1)
    set_random_state(regressor_2)

    if name in CROSS_DECOMPOSITION:
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y

    # fit
    regressor_1.fit(X, y_)
    pred1 = regressor_1.predict(X)
    regressor_2.fit(X, y_.astype(float))
    pred2 = regressor_2.predict(X)
    assert_allclose(pred1, pred2, atol=1e-2, err_msg=name)


@ignore_warnings(category=FutureWarning)
def check_regressors_train(
    name, regressor_orig, readonly_memmap=False, X_dtype=np.float64
):
    X, y = _regression_dataset()
    X = X.astype(X_dtype)
    y = scale(y)  # X is already scaled
    regressor = clone(regressor_orig)
    X = _enforce_estimator_tags_X(regressor, X)
    y = _enforce_estimator_tags_y(regressor, y)
    if name in CROSS_DECOMPOSITION:
        rnd = np.random.RandomState(0)
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y

    if readonly_memmap:
        X, y, y_ = create_memmap_backed_data([X, y, y_])

    if not hasattr(regressor, "alphas") and hasattr(regressor, "alpha"):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01
    if name == "PassiveAggressiveRegressor":
        regressor.C = 0.01

    # raises error on malformed input for fit
    with raises(
        ValueError,
        err_msg=(
            f"The classifier {name} does not raise an error when "
            "incorrect/malformed input data for fit is passed. The number of "
            "training examples is not the same as the number of labels. Perhaps "
            "use check_X_y in fit."
        ),
    ):
        regressor.fit(X, y[:-1])
    # fit
    set_random_state(regressor)
    regressor.fit(X, y_)
    regressor.fit(X.tolist(), y_.tolist())
    y_pred = regressor.predict(X)
    assert y_pred.shape == y_.shape

    # TODO: find out why PLS and CCA fail. RANSAC is random
    # and furthermore assumes the presence of outliers, hence
    # skipped
    if not get_tags(regressor).regressor_tags.poor_score:
        assert regressor.score(X, y_) > 0.5


@ignore_warnings
def check_regressors_no_decision_function(name, regressor_orig):
    # check that regressors don't have a decision_function, predict_proba, or
    # predict_log_proba method.
    rng = np.random.RandomState(0)
    regressor = clone(regressor_orig)

    X = rng.normal(size=(10, 4))
    X = _enforce_estimator_tags_X(regressor_orig, X)
    y = _enforce_estimator_tags_y(regressor, X[:, 0])

    regressor.fit(X, y)
    funcs = ["decision_function", "predict_proba", "predict_log_proba"]
    for func_name in funcs:
        assert not hasattr(regressor, func_name)


@ignore_warnings(category=FutureWarning)
def check_class_weight_classifiers(name, classifier_orig):
    if get_tags(classifier_orig).classifier_tags.multi_class:
        problems = [2, 3]
    else:  # binary classification only
        problems = [2]

    for n_centers in problems:
        # create a very noisy dataset
        X, y = make_blobs(centers=n_centers, random_state=0, cluster_std=20)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.5, random_state=0
        )

        # can't use gram_if_pairwise() here, setting up gram matrix manually
        if get_tags(classifier_orig).input_tags.pairwise:
            X_test = rbf_kernel(X_test, X_train)
            X_train = rbf_kernel(X_train, X_train)

        n_centers = len(np.unique(y_train))

        if n_centers == 2:
            class_weight = {0: 1000, 1: 0.0001}
        else:
            class_weight = {0: 1000, 1: 0.0001, 2: 0.0001}

        classifier = clone(classifier_orig).set_params(class_weight=class_weight)
        if hasattr(classifier, "n_iter"):
            classifier.set_params(n_iter=100)
        if hasattr(classifier, "max_iter"):
            classifier.set_params(max_iter=1000)
        if hasattr(classifier, "min_weight_fraction_leaf"):
            classifier.set_params(min_weight_fraction_leaf=0.01)
        if hasattr(classifier, "n_iter_no_change"):
            classifier.set_params(n_iter_no_change=20)

        set_random_state(classifier)
        classifier.fit(X_train, y_train)
        y_pred = classifier.predict(X_test)
        # XXX: Generally can use 0.89 here. On Windows, LinearSVC gets
        #      0.88 (Issue #9111)
        if not get_tags(classifier_orig).classifier_tags.poor_score:
            assert np.mean(y_pred == 0) > 0.87


@ignore_warnings(category=FutureWarning)
def check_class_weight_balanced_classifiers(
    name, classifier_orig, X_train, y_train, X_test, y_test, weights
):
    classifier = clone(classifier_orig)
    if hasattr(classifier, "n_iter"):
        classifier.set_params(n_iter=100)
    if hasattr(classifier, "max_iter"):
        classifier.set_params(max_iter=1000)

    set_random_state(classifier)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)

    classifier.set_params(class_weight="balanced")
    classifier.fit(X_train, y_train)
    y_pred_balanced = classifier.predict(X_test)
    assert f1_score(y_test, y_pred_balanced, average="weighted") > f1_score(
        y_test, y_pred, average="weighted"
    )


@ignore_warnings(category=FutureWarning)
def check_class_weight_balanced_linear_classifier(name, estimator_orig):
    """Test class weights with non-contiguous class labels."""
    # this is run on classes, not instances, though this should be changed
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-0.8, -1.0], [1.0, 1.0], [1.0, 0.0]])
    y = np.array([1, 1, 1, -1, -1])

    classifier = clone(estimator_orig)

    if hasattr(classifier, "n_iter"):
        # This is a very small dataset, default n_iter are likely to prevent
        # convergence
        classifier.set_params(n_iter=1000)
    if hasattr(classifier, "max_iter"):
        classifier.set_params(max_iter=1000)
    if hasattr(classifier, "cv"):
        classifier.set_params(cv=3)
    set_random_state(classifier)

    # Let the model compute the class frequencies
    classifier.set_params(class_weight="balanced")
    coef_balanced = classifier.fit(X, y).coef_.copy()

    # Count each label occurrence to reweight manually
    n_samples = len(y)
    n_classes = float(len(np.unique(y)))

    class_weight = {
        1: n_samples / (np.sum(y == 1) * n_classes),
        -1: n_samples / (np.sum(y == -1) * n_classes),
    }
    classifier.set_params(class_weight=class_weight)
    coef_manual = classifier.fit(X, y).coef_.copy()

    assert_allclose(
        coef_balanced,
        coef_manual,
        err_msg="Classifier %s is not computing class_weight=balanced properly." % name,
    )


@ignore_warnings(category=FutureWarning)
def check_estimators_overwrite_params(name, estimator_orig):
    X, y = make_blobs(random_state=0, n_samples=21)
    X = _enforce_estimator_tags_X(estimator_orig, X, kernel=rbf_kernel)
    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    set_random_state(estimator)

    # Make a physical copy of the original estimator parameters before fitting.
    params = estimator.get_params()
    original_params = deepcopy(params)

    # Fit the model
    estimator.fit(X, y)

    # Compare the state of the model parameters with the original parameters
    new_params = estimator.get_params()
    for param_name, original_value in original_params.items():
        new_value = new_params[param_name]

        # We should never change or mutate the internal state of input
        # parameters by default. To check this we use the joblib.hash function
        # that introspects recursively any subobjects to compute a checksum.
        # The only exception to this rule of immutable constructor parameters
        # is possible RandomState instance but in this check we explicitly
        # fixed the random_state params recursively to be integer seeds.
        assert joblib.hash(new_value) == joblib.hash(original_value), (
            "Estimator %s should not change or mutate "
            " the parameter %s from %s to %s during fit."
            % (name, param_name, original_value, new_value)
        )


@ignore_warnings(category=FutureWarning)
def check_no_attributes_set_in_init(name, estimator_orig):
    """Check setting during init."""
    try:
        # Clone fails if the estimator does not store
        # all parameters as an attribute during init
        estimator = clone(estimator_orig)
    except AttributeError:
        raise AttributeError(
            f"Estimator {name} should store all parameters as an attribute during init."
        )

    if hasattr(type(estimator).__init__, "deprecated_original"):
        return

    init_params = _get_args(type(estimator).__init__)
    parents_init_params = [
        param
        for params_parent in (_get_args(parent) for parent in type(estimator).__mro__)
        for param in params_parent
    ]

    # Test for no setting apart from parameters during init
    invalid_attr = set(vars(estimator)) - set(init_params) - set(parents_init_params)
    # Ignore private attributes
    invalid_attr = set([attr for attr in invalid_attr if not attr.startswith("_")])
    assert not invalid_attr, (
        "Estimator %s should not set any attribute apart"
        " from parameters during init. Found attributes %s."
        % (name, sorted(invalid_attr))
    )


@ignore_warnings(category=FutureWarning)
def check_sparsify_coefficients(name, estimator_orig):
    X = np.array(
        [
            [-2, -1],
            [-1, -1],
            [-1, -2],
            [1, 1],
            [1, 2],
            [2, 1],
            [-1, -2],
            [2, 2],
            [-2, -2],
        ]
    )
    y = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3])
    y = _enforce_estimator_tags_y(estimator_orig, y)
    est = clone(estimator_orig)

    est.fit(X, y)
    pred_orig = est.predict(X)

    # test sparsify with dense inputs
    est.sparsify()
    assert sparse.issparse(est.coef_)
    pred = est.predict(X)
    assert_array_equal(pred, pred_orig)

    # pickle and unpickle with sparse coef_
    est = pickle.loads(pickle.dumps(est))
    assert sparse.issparse(est.coef_)
    pred = est.predict(X)
    assert_array_equal(pred, pred_orig)


@ignore_warnings(category=FutureWarning)
def check_classifier_data_not_an_array(name, estimator_orig):
    X = np.array(
        [
            [3, 0],
            [0, 1],
            [0, 2],
            [1, 1],
            [1, 2],
            [2, 1],
            [0, 3],
            [1, 0],
            [2, 0],
            [4, 4],
            [2, 3],
            [3, 2],
        ]
    )
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = np.array([1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2])
    y = _enforce_estimator_tags_y(estimator_orig, y)
    for obj_type in ["NotAnArray", "PandasDataframe"]:
        check_estimators_data_not_an_array(name, estimator_orig, X, y, obj_type)


@ignore_warnings(category=FutureWarning)
def check_regressor_data_not_an_array(name, estimator_orig):
    X, y = _regression_dataset()
    X = _enforce_estimator_tags_X(estimator_orig, X)
    y = _enforce_estimator_tags_y(estimator_orig, y)
    for obj_type in ["NotAnArray", "PandasDataframe"]:
        check_estimators_data_not_an_array(name, estimator_orig, X, y, obj_type)


@ignore_warnings(category=FutureWarning)
def check_estimators_data_not_an_array(name, estimator_orig, X, y, obj_type):
    if name in CROSS_DECOMPOSITION:
        raise SkipTest(
            "Skipping check_estimators_data_not_an_array "
            "for cross decomposition module as estimators "
            "are not deterministic."
        )
    # separate estimators to control random seeds
    estimator_1 = clone(estimator_orig)
    estimator_2 = clone(estimator_orig)
    set_random_state(estimator_1)
    set_random_state(estimator_2)

    if obj_type not in ["NotAnArray", "PandasDataframe"]:
        raise ValueError("Data type {0} not supported".format(obj_type))

    if obj_type == "NotAnArray":
        y_ = _NotAnArray(np.asarray(y))
        X_ = _NotAnArray(np.asarray(X))
    else:
        # Here pandas objects (Series and DataFrame) are tested explicitly
        # because some estimators may handle them (especially their indexing)
        # specially.
        try:
            import pandas as pd

            y_ = np.asarray(y)
            if y_.ndim == 1:
                y_ = pd.Series(y_, copy=False)
            else:
                y_ = pd.DataFrame(y_, copy=False)
            X_ = pd.DataFrame(np.asarray(X), copy=False)

        except ImportError:
            raise SkipTest(
                "pandas is not installed: not checking estimators for pandas objects."
            )

    # fit
    estimator_1.fit(X_, y_)
    pred1 = estimator_1.predict(X_)
    estimator_2.fit(X, y)
    pred2 = estimator_2.predict(X)
    assert_allclose(pred1, pred2, atol=1e-2, err_msg=name)


def check_estimator_cloneable(name, estimator_orig):
    """Checks whether the estimator can be cloned."""
    try:
        clone(estimator_orig)
    except Exception as e:
        raise AssertionError(f"Cloning of {name} failed with error: {e}.") from e


def check_estimator_repr(name, estimator_orig):
    """Check that the estimator has a functioning repr."""
    estimator = clone(estimator_orig)
    try:
        repr(estimator)
    except Exception as e:
        raise AssertionError(f"Repr of {name} failed with error: {e}.") from e


def check_parameters_default_constructible(name, estimator_orig):
    # test default-constructibility
    # get rid of deprecation warnings

    Estimator = estimator_orig.__class__
    estimator = clone(estimator_orig)

    with ignore_warnings(category=FutureWarning):
        # test that set_params returns self
        # TODO(devtools): this should be a separate check.
        assert estimator.set_params() is estimator

        # test if init does nothing but set parameters
        # this is important for grid_search etc.
        # We get the default parameters from init and then
        # compare these against the actual values of the attributes.

        # this comes from getattr. Gets rid of deprecation decorator.
        init = getattr(estimator.__init__, "deprecated_original", estimator.__init__)

        try:

            def param_default_value(p):
                """Identify hyper parameters of an estimator."""
                return (
                    p.name != "self"
                    and p.kind != p.VAR_KEYWORD
                    and p.kind != p.VAR_POSITIONAL
                    # and it should have a default value for this test
                    and p.default != p.empty
                )

            def param_required(p):
                """Identify hyper parameters of an estimator."""
                return (
                    p.name != "self"
                    and p.kind != p.VAR_KEYWORD
                    # technically VAR_POSITIONAL is also required, but we don't have a
                    # nice way to check for it. We assume there's no VAR_POSITIONAL in
                    # the constructor parameters.
                    #
                    # TODO(devtools): separately check that the constructor doesn't
                    # have *args.
                    and p.kind != p.VAR_POSITIONAL
                    # these are parameters that don't have a default value and are
                    # required to construct the estimator.
                    and p.default == p.empty
                )

            required_params_names = [
                p.name for p in signature(init).parameters.values() if param_required(p)
            ]

            default_value_params = [
                p for p in signature(init).parameters.values() if param_default_value(p)
            ]

        except (TypeError, ValueError):
            # init is not a python function.
            # true for mixins
            return

        # here we construct an instance of the estimator using only the required
        # parameters.
        old_params = estimator.get_params()
        init_params = {
            param: old_params[param]
            for param in old_params
            if param in required_params_names
        }
        estimator = Estimator(**init_params)
        params = estimator.get_params()

        for init_param in default_value_params:
            allowed_types = {
                str,
                int,
                float,
                bool,
                tuple,
                type(None),
                type,
            }
            # Any numpy numeric such as np.int32.
            allowed_types.update(np.sctypeDict.values())

            allowed_value = (
                type(init_param.default) in allowed_types
                or
                # Although callables are mutable, we accept them as argument
                # default value and trust that neither the implementation of
                # the callable nor of the estimator changes the state of the
                # callable.
                callable(init_param.default)
            )

            assert allowed_value, (
                f"Parameter '{init_param.name}' of estimator "
                f"'{Estimator.__name__}' is of type "
                f"{type(init_param.default).__name__} which is not allowed. "
                f"'{init_param.name}' must be a callable or must be of type "
                f"{set(type.__name__ for type in allowed_types)}."
            )
            if init_param.name not in params.keys():
                # deprecated parameter, not in get_params
                assert init_param.default is None, (
                    f"Estimator parameter '{init_param.name}' of estimator "
                    f"'{Estimator.__name__}' is not returned by get_params. "
                    "If it is deprecated, set its default value to None."
                )
                continue

            param_value = params[init_param.name]
            if isinstance(param_value, np.ndarray):
                assert_array_equal(param_value, init_param.default)
            else:
                failure_text = (
                    f"Parameter {init_param.name} was mutated on init. All "
                    "parameters must be stored unchanged."
                )
                if is_scalar_nan(param_value):
                    # Allows to set default parameters to np.nan
                    assert param_value is init_param.default, failure_text
                else:
                    assert param_value == init_param.default, failure_text


def _enforce_estimator_tags_y(estimator, y):
    # Estimators with a `requires_positive_y` tag only accept strictly positive
    # data
    tags = get_tags(estimator)
    if tags.target_tags.positive_only:
        # Create strictly positive y. The minimal increment above 0 is 1, as
        # y could be of integer dtype.
        y += 1 + abs(y.min())
    if (
        tags.classifier_tags is not None
        and not tags.classifier_tags.multi_class
        and y.size > 0
    ):
        y = np.where(y == y.min(), y, y.min() + 1)
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if tags.target_tags.multi_output and not tags.target_tags.single_output:
        return np.reshape(y, (-1, 1))
    return y


def _enforce_estimator_tags_X(estimator, X, X_test=None, kernel=linear_kernel):
    # Estimators with `1darray` in `X_types` tag only accept
    # X of shape (`n_samples`,)
    if get_tags(estimator).input_tags.one_d_array:
        X = X[:, 0]
        if X_test is not None:
            X_test = X_test[:, 0]  # pragma: no cover
    # Estimators with a `requires_positive_X` tag only accept
    # strictly positive data
    if get_tags(estimator).input_tags.positive_only:
        X = X - X.min()
        if X_test is not None:
            X_test = X_test - X_test.min()  # pragma: no cover
    if get_tags(estimator).input_tags.categorical:
        dtype = np.float64 if get_tags(estimator).input_tags.allow_nan else np.int32
        X = np.round((X - X.min())).astype(dtype)
        if X_test is not None:
            X_test = np.round((X_test - X_test.min())).astype(dtype)  # pragma: no cover

    if estimator.__class__.__name__ == "SkewedChi2Sampler":
        # SkewedChi2Sampler requires X > -skewdness in transform
        X = X - X.min()
        if X_test is not None:
            X_test = X_test - X_test.min()  # pragma: no cover

    X_res = X

    # Pairwise estimators only accept
    # X of shape (`n_samples`, `n_samples`)
    if _is_pairwise_metric(estimator):
        X_res = pairwise_distances(X, metric="euclidean")
        if X_test is not None:
            X_test = pairwise_distances(
                X_test, X, metric="euclidean"
            )  # pragma: no cover
    elif get_tags(estimator).input_tags.pairwise:
        X_res = kernel(X, X)
        if X_test is not None:
            X_test = kernel(X_test, X)  # pragma: no cover
    if X_test is not None:
        return X_res, X_test
    return X_res


@ignore_warnings(category=FutureWarning)
def check_positive_only_tag_during_fit(name, estimator_orig):
    """Test that the estimator correctly sets the tags.input_tags.positive_only

    If the tag is False, the estimator should accept negative input regardless of the
    tags.input_tags.pairwise flag.
    """
    estimator = clone(estimator_orig)
    tags = get_tags(estimator)

    X, y = load_iris(return_X_y=True)
    y = _enforce_estimator_tags_y(estimator, y)
    set_random_state(estimator, 0)
    X = _enforce_estimator_tags_X(estimator, X)
    # Make sure that the dtype of X stays unchanged: for instance estimator
    # that expect categorical inputs typically expected integer-based encoded
    # categories.
    X -= X.mean().astype(X.dtype)

    if tags.input_tags.positive_only:
        with raises(ValueError, match="Negative values in data"):
            estimator.fit(X, y)
    else:
        # This should pass
        try:
            estimator.fit(X, y)
        except Exception as e:
            err_msg = (
                f"Estimator {name!r} raised {e.__class__.__name__} unexpectedly."
                " This happens when passing negative input values as X."
                " If negative values are not supported for this estimator instance,"
                " then the tags.input_tags.positive_only tag needs to be set to True."
            )
            raise AssertionError(err_msg) from e


@ignore_warnings(category=FutureWarning)
def check_non_transformer_estimators_n_iter(name, estimator_orig):
    # Test that estimators that are not transformers with a parameter
    # max_iter, return the attribute of n_iter_ at least 1.

    if not hasattr(estimator_orig, "max_iter"):
        return

    estimator = clone(estimator_orig)
    iris = load_iris()
    X, y_ = iris.data, iris.target
    y_ = _enforce_estimator_tags_y(estimator, y_)
    set_random_state(estimator, 0)
    X = _enforce_estimator_tags_X(estimator_orig, X)

    estimator.fit(X, y_)

    assert np.all(np.asarray(estimator.n_iter_) >= 1), (
        "Estimators with a `max_iter` parameter, should expose an `n_iter_` attribute,"
        " indicating the number of iterations that were executed. The values in the "
        "`n_iter_` attribute should be greater or equal to 1."
    )


@ignore_warnings(category=FutureWarning)
def check_transformer_n_iter(name, estimator_orig):
    # Test that transformers with a parameter max_iter, return the
    # attribute of n_iter_ at least 1.
    estimator = clone(estimator_orig)
    if hasattr(estimator, "max_iter"):
        if name in CROSS_DECOMPOSITION:
            # Check using default data
            X = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [2.0, 2.0, 2.0], [2.0, 5.0, 4.0]]
            y_ = [[0.1, -0.2], [0.9, 1.1], [0.1, -0.5], [0.3, -0.2]]

        else:
            X, y_ = make_blobs(
                n_samples=30,
                centers=[[0, 0, 0], [1, 1, 1]],
                random_state=0,
                n_features=2,
                cluster_std=0.1,
            )
            X = _enforce_estimator_tags_X(estimator_orig, X)
        set_random_state(estimator, 0)
        estimator.fit(X, y_)

        # These return a n_iter per component.
        if name in CROSS_DECOMPOSITION:
            for iter_ in estimator.n_iter_:
                assert iter_ >= 1
        else:
            assert estimator.n_iter_ >= 1


@ignore_warnings(category=FutureWarning)
def check_get_params_invariance(name, estimator_orig):
    # Checks if get_params(deep=False) is a subset of get_params(deep=True)
    e = clone(estimator_orig)

    shallow_params = e.get_params(deep=False)
    deep_params = e.get_params(deep=True)

    assert all(item in deep_params.items() for item in shallow_params.items())


@ignore_warnings(category=FutureWarning)
def check_set_params(name, estimator_orig):
    # Check that get_params() returns the same thing
    # before and after set_params() with some fuzz
    estimator = clone(estimator_orig)

    orig_params = estimator.get_params(deep=False)
    msg = "get_params result does not match what was passed to set_params"

    estimator.set_params(**orig_params)
    curr_params = estimator.get_params(deep=False)
    assert set(orig_params.keys()) == set(curr_params.keys()), msg
    for k, v in curr_params.items():
        assert orig_params[k] is v, msg

    # some fuzz values
    test_values = [-np.inf, np.inf, None]

    test_params = deepcopy(orig_params)
    for param_name in orig_params.keys():
        default_value = orig_params[param_name]
        for value in test_values:
            test_params[param_name] = value
            try:
                estimator.set_params(**test_params)
            except (TypeError, ValueError) as e:
                e_type = e.__class__.__name__
                # Exception occurred, possibly parameter validation
                warnings.warn(
                    "{0} occurred during set_params of param {1} on "
                    "{2}. It is recommended to delay parameter "
                    "validation until fit.".format(e_type, param_name, name)
                )

                change_warning_msg = (
                    "Estimator's parameters changed after set_params raised {}".format(
                        e_type
                    )
                )
                params_before_exception = curr_params
                curr_params = estimator.get_params(deep=False)
                try:
                    assert set(params_before_exception.keys()) == set(
                        curr_params.keys()
                    )
                    for k, v in curr_params.items():
                        assert params_before_exception[k] is v
                except AssertionError:
                    warnings.warn(change_warning_msg)
            else:
                curr_params = estimator.get_params(deep=False)
                assert set(test_params.keys()) == set(curr_params.keys()), msg
                for k, v in curr_params.items():
                    assert test_params[k] is v, msg
        test_params[param_name] = default_value


@ignore_warnings(category=FutureWarning)
def check_classifiers_regression_target(name, estimator_orig):
    # Check if classifier throws an exception when fed regression targets

    X, y = _regression_dataset()

    X = _enforce_estimator_tags_X(estimator_orig, X)
    e = clone(estimator_orig)
    err_msg = (
        "When a classifier is passed a continuous target, it should raise a ValueError"
        " with a message containing 'Unknown label type: ' or a message indicating that"
        " a continuous target is passed and the message should include the word"
        " 'continuous'"
    )
    msg = "Unknown label type: |continuous"
    if not get_tags(e).no_validation:
        with raises(ValueError, match=msg, err_msg=err_msg):
            e.fit(X, y)


@ignore_warnings(category=FutureWarning)
def check_decision_proba_consistency(name, estimator_orig):
    # Check whether an estimator having both decision_function and
    # predict_proba methods has outputs with perfect rank correlation.

    centers = [(2, 2), (4, 4)]
    X, y = make_blobs(
        n_samples=100,
        random_state=0,
        n_features=4,
        centers=centers,
        cluster_std=1.0,
        shuffle=True,
    )
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=0
    )
    estimator = clone(estimator_orig)

    if hasattr(estimator, "decision_function") and hasattr(estimator, "predict_proba"):
        estimator.fit(X_train, y_train)
        # Since the link function from decision_function() to predict_proba()
        # is sometimes not precise enough (typically expit), we round to the
        # 10th decimal to avoid numerical issues: we compare the rank
        # with deterministic ties rather than get platform specific rank
        # inversions in case of machine level differences.
        a = estimator.predict_proba(X_test)[:, 1].round(decimals=10)
        b = estimator.decision_function(X_test).round(decimals=10)

        rank_proba, rank_score = rankdata(a), rankdata(b)
        try:
            assert_array_almost_equal(rank_proba, rank_score)
        except AssertionError:
            # Sometimes, the rounding applied on the probabilities will have
            # ties that are not present in the scores because it is
            # numerically more precise. In this case, we relax the test by
            # grouping the decision function scores based on the probability
            # rank and check that the score is monotonically increasing.
            grouped_y_score = np.array(
                [b[rank_proba == group].mean() for group in np.unique(rank_proba)]
            )
            sorted_idx = np.argsort(grouped_y_score)
            assert_array_equal(sorted_idx, np.arange(len(sorted_idx)))


def check_outliers_fit_predict(name, estimator_orig):
    # Check fit_predict for outlier detectors.

    n_samples = 300
    X, _ = make_blobs(n_samples=n_samples, random_state=0)
    X = shuffle(X, random_state=7)
    n_samples, n_features = X.shape
    estimator = clone(estimator_orig)

    set_random_state(estimator)

    y_pred = estimator.fit_predict(X)
    assert y_pred.shape == (n_samples,)
    assert y_pred.dtype.kind == "i"
    assert_array_equal(np.unique(y_pred), np.array([-1, 1]))

    # check fit_predict = fit.predict when the estimator has both a predict and
    # a fit_predict method. recall that it is already assumed here that the
    # estimator has a fit_predict method
    if hasattr(estimator, "predict"):
        y_pred_2 = estimator.fit(X).predict(X)
        assert_array_equal(y_pred, y_pred_2)

    if hasattr(estimator, "contamination"):
        # proportion of outliers equal to contamination parameter when not
        # set to 'auto'
        expected_outliers = 30
        contamination = float(expected_outliers) / n_samples
        estimator.set_params(contamination=contamination)
        y_pred = estimator.fit_predict(X)

        num_outliers = np.sum(y_pred != 1)
        # num_outliers should be equal to expected_outliers unless
        # there are ties in the decision_function values. this can
        # only be tested for estimators with a decision_function
        # method
        if num_outliers != expected_outliers and hasattr(
            estimator, "decision_function"
        ):
            decision = estimator.decision_function(X)
            check_outlier_corruption(num_outliers, expected_outliers, decision)


def check_fit_non_negative(name, estimator_orig):
    # Check that proper warning is raised for non-negative X
    # when tag requires_positive_X is present
    X = np.array([[-1.0, 1], [-1.0, 1]])
    y = np.array([1, 2])
    estimator = clone(estimator_orig)
    with raises(ValueError):
        estimator.fit(X, y)


def check_fit_idempotent(name, estimator_orig):
    # Check that est.fit(X) is the same as est.fit(X).fit(X). Ideally we would
    # check that the estimated parameters during training (e.g. coefs_) are
    # the same, but having a universal comparison function for those
    # attributes is difficult and full of edge cases. So instead we check that
    # predict(), predict_proba(), decision_function() and transform() return
    # the same results.

    check_methods = ["predict", "transform", "decision_function", "predict_proba"]
    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)
    if "warm_start" in estimator.get_params().keys():
        estimator.set_params(warm_start=False)

    n_samples = 100
    X = rng.normal(loc=100, size=(n_samples, 2))
    X = _enforce_estimator_tags_X(estimator, X)
    if is_regressor(estimator_orig):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = _enforce_estimator_tags_y(estimator, y)

    train, test = next(ShuffleSplit(test_size=0.2, random_state=rng).split(X))
    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)

    # Fit for the first time
    estimator.fit(X_train, y_train)

    result = {
        method: getattr(estimator, method)(X_test)
        for method in check_methods
        if hasattr(estimator, method)
    }

    # Fit again
    set_random_state(estimator)
    estimator.fit(X_train, y_train)

    for method in check_methods:
        if hasattr(estimator, method):
            new_result = getattr(estimator, method)(X_test)
            if hasattr(new_result, "dtype") and np.issubdtype(
                new_result.dtype, np.floating
            ):
                tol = 2 * np.finfo(new_result.dtype).eps
            else:
                tol = 2 * np.finfo(np.float64).eps
            assert_allclose_dense_sparse(
                result[method],
                new_result,
                atol=max(tol, 1e-9),
                rtol=max(tol, 1e-7),
                err_msg="Idempotency check failed for method {}".format(method),
            )


def check_fit_check_is_fitted(name, estimator_orig):
    # Make sure that estimator doesn't pass check_is_fitted before calling fit
    # and that passes check_is_fitted once it's fit.

    rng = np.random.RandomState(42)

    estimator = clone(estimator_orig)
    set_random_state(estimator)
    if "warm_start" in estimator.get_params():
        estimator.set_params(warm_start=False)

    n_samples = 100
    X = rng.normal(loc=100, size=(n_samples, 2))
    X = _enforce_estimator_tags_X(estimator, X)
    if is_regressor(estimator_orig):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = _enforce_estimator_tags_y(estimator, y)

    if get_tags(estimator).requires_fit:
        # stateless estimators (such as FunctionTransformer) are always "fit"!
        try:
            check_is_fitted(estimator)
            raise AssertionError(
                f"{estimator.__class__.__name__} passes check_is_fitted before being"
                " fit!"
            )
        except NotFittedError:
            pass
    estimator.fit(X, y)
    try:
        check_is_fitted(estimator)
    except NotFittedError as e:
        raise NotFittedError(
            "Estimator fails to pass `check_is_fitted` even though it has been fit."
        ) from e


def check_n_features_in(name, estimator_orig):
    # Make sure that n_features_in_ attribute doesn't exist until fit is
    # called, and that its value is correct.

    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)
    if "warm_start" in estimator.get_params():
        estimator.set_params(warm_start=False)

    n_samples = 100
    X = rng.normal(loc=100, size=(n_samples, 2))
    X = _enforce_estimator_tags_X(estimator, X)
    if is_regressor(estimator_orig):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = _enforce_estimator_tags_y(estimator, y)

    assert not hasattr(estimator, "n_features_in_")
    estimator.fit(X, y)
    assert hasattr(estimator, "n_features_in_")
    assert estimator.n_features_in_ == X.shape[1]


def check_requires_y_none(name, estimator_orig):
    # Make sure that an estimator with requires_y=True fails gracefully when
    # given y=None

    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)

    n_samples = 100
    X = rng.normal(loc=100, size=(n_samples, 2))
    X = _enforce_estimator_tags_X(estimator, X)

    expected_err_msgs = (
        "requires y to be passed, but the target y is None",
        "Expected array-like (array or non-string sequence), got None",
        "y should be a 1d array",
    )

    try:
        estimator.fit(X, None)
    except ValueError as ve:
        if not any(msg in str(ve) for msg in expected_err_msgs):
            raise ve


@ignore_warnings(category=FutureWarning)
def check_n_features_in_after_fitting(name, estimator_orig):
    # Make sure that n_features_in are checked after fitting
    tags = get_tags(estimator_orig)

    is_supported_X_types = tags.input_tags.two_d_array or tags.input_tags.categorical

    if not is_supported_X_types or tags.no_validation:
        return

    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)
    if "warm_start" in estimator.get_params():
        estimator.set_params(warm_start=False)

    n_samples = 10
    X = rng.normal(size=(n_samples, 4))
    X = _enforce_estimator_tags_X(estimator, X)

    if is_regressor(estimator):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = _enforce_estimator_tags_y(estimator, y)

    err_msg = (
        "`{name}.fit()` does not set the `n_features_in_` attribute. "
        "You might want to use `sklearn.utils.validation.validate_data` instead "
        "of `check_array` in `{name}.fit()` which takes care of setting the "
        "attribute.".format(name=name)
    )

    estimator.fit(X, y)
    assert hasattr(estimator, "n_features_in_"), err_msg
    assert estimator.n_features_in_ == X.shape[1], err_msg

    # check methods will check n_features_in_
    check_methods = [
        "predict",
        "transform",
        "decision_function",
        "predict_proba",
        "score",
    ]
    X_bad = X[:, [1]]

    err_msg = """\
        `{name}.{method}()` does not check for consistency between input number
        of features with {name}.fit(), via the `n_features_in_` attribute.
        You might want to use `sklearn.utils.validation.validate_data` instead
        of `check_array` in `{name}.fit()` and {name}.{method}()`. This can be done
        like the following:
        from sklearn.utils.validation import validate_data
        ...
        class MyEstimator(BaseEstimator):
            ...
            def fit(self, X, y):
                X, y = validate_data(self, X, y, ...)
                ...
                return self
            ...
            def {method}(self, X):
                X = validate_data(self, X, ..., reset=False)
                ...
            return X
    """
    err_msg = textwrap.dedent(err_msg)

    msg = f"X has 1 features, but \\w+ is expecting {X.shape[1]} features as input"
    for method in check_methods:
        if not hasattr(estimator, method):
            continue

        callable_method = getattr(estimator, method)
        if method == "score":
            callable_method = partial(callable_method, y=y)

        with raises(
            ValueError, match=msg, err_msg=err_msg.format(name=name, method=method)
        ):
            callable_method(X_bad)

    # partial_fit will check in the second call
    if not hasattr(estimator, "partial_fit"):
        return

    estimator = clone(estimator_orig)
    if is_classifier(estimator):
        estimator.partial_fit(X, y, classes=np.unique(y))
    else:
        estimator.partial_fit(X, y)
    assert estimator.n_features_in_ == X.shape[1]

    with raises(ValueError, match=msg):
        estimator.partial_fit(X_bad, y)


def check_valid_tag_types(name, estimator):
    """Check that estimator tags are valid."""
    assert hasattr(estimator, "__sklearn_tags__"), (
        f"Estimator {name} does not have `__sklearn_tags__` method. This method is"
        " implemented in BaseEstimator and returns a sklearn.utils.Tags instance."
    )
    err_msg = (
        "Tag values need to be of a certain type. "
        "Please refer to the documentation of `sklearn.utils.Tags` for more details."
    )
    tags = get_tags(estimator)
    assert isinstance(tags.estimator_type, (str, type(None))), err_msg
    assert isinstance(tags.target_tags, TargetTags), err_msg
    assert isinstance(tags.classifier_tags, (ClassifierTags, type(None))), err_msg
    assert isinstance(tags.regressor_tags, (RegressorTags, type(None))), err_msg
    assert isinstance(tags.transformer_tags, (TransformerTags, type(None))), err_msg
    assert isinstance(tags.input_tags, InputTags), err_msg
    assert isinstance(tags.array_api_support, bool), err_msg
    assert isinstance(tags.no_validation, bool), err_msg
    assert isinstance(tags.non_deterministic, bool), err_msg
    assert isinstance(tags.requires_fit, bool), err_msg
    assert isinstance(tags._skip_test, bool), err_msg

    assert isinstance(tags.target_tags.required, bool), err_msg
    assert isinstance(tags.target_tags.one_d_labels, bool), err_msg
    assert isinstance(tags.target_tags.two_d_labels, bool), err_msg
    assert isinstance(tags.target_tags.positive_only, bool), err_msg
    assert isinstance(tags.target_tags.multi_output, bool), err_msg
    assert isinstance(tags.target_tags.single_output, bool), err_msg

    assert isinstance(tags.input_tags.pairwise, bool), err_msg
    assert isinstance(tags.input_tags.allow_nan, bool), err_msg
    assert isinstance(tags.input_tags.sparse, bool), err_msg
    assert isinstance(tags.input_tags.categorical, bool), err_msg
    assert isinstance(tags.input_tags.string, bool), err_msg
    assert isinstance(tags.input_tags.dict, bool), err_msg
    assert isinstance(tags.input_tags.one_d_array, bool), err_msg
    assert isinstance(tags.input_tags.two_d_array, bool), err_msg
    assert isinstance(tags.input_tags.three_d_array, bool), err_msg
    assert isinstance(tags.input_tags.positive_only, bool), err_msg

    if tags.classifier_tags is not None:
        assert isinstance(tags.classifier_tags.poor_score, bool), err_msg
        assert isinstance(tags.classifier_tags.multi_class, bool), err_msg
        assert isinstance(tags.classifier_tags.multi_label, bool), err_msg

    if tags.regressor_tags is not None:
        assert isinstance(tags.regressor_tags.poor_score, bool), err_msg

    if tags.transformer_tags is not None:
        assert isinstance(tags.transformer_tags.preserves_dtype, list), err_msg


def check_estimator_tags_renamed(name, estimator_orig):
    help = """{tags_func}() was removed in 1.6. Please use __sklearn_tags__ instead.
You can implement both __sklearn_tags__() and {tags_func}() to support multiple
scikit-learn versions.
"""

    for klass in type(estimator_orig).mro():
        if (
            # Here we check vars(...) because we want to check if the method is
            # explicitly defined in the class instead of inherited from a parent class.
            ("_more_tags" in vars(klass) or "_get_tags" in vars(klass))
            and "__sklearn_tags__" not in vars(klass)
        ):
            raise TypeError(
                f"Estimator {name} has defined either `_more_tags` or `_get_tags`,"
                " but not `__sklearn_tags__`. If you're customizing tags, and need to"
                " support multiple scikit-learn versions, you can implement both"
                " `__sklearn_tags__` and `_more_tags` or `_get_tags`. This change was"
                " introduced in scikit-learn=1.6"
            )


def check_dataframe_column_names_consistency(name, estimator_orig):
    try:
        import pandas as pd
    except ImportError:
        raise SkipTest(
            "pandas is not installed: not checking column name consistency for pandas"
        )

    tags = get_tags(estimator_orig)
    is_supported_X_types = tags.input_tags.two_d_array or tags.input_tags.categorical

    if not is_supported_X_types or tags.no_validation:
        return

    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)

    X_orig = rng.normal(size=(150, 8))

    X_orig = _enforce_estimator_tags_X(estimator, X_orig)
    n_samples, n_features = X_orig.shape

    names = np.array([f"col_{i}" for i in range(n_features)])
    X = pd.DataFrame(X_orig, columns=names, copy=False)

    if is_regressor(estimator):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = _enforce_estimator_tags_y(estimator, y)

    # Check that calling `fit` does not raise any warnings about feature names.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "error",
            message="X does not have valid feature names",
            category=UserWarning,
            module="sklearn",
        )
        estimator.fit(X, y)

    if not hasattr(estimator, "feature_names_in_"):
        raise ValueError(
            "Estimator does not have a feature_names_in_ "
            "attribute after fitting with a dataframe"
        )
    assert isinstance(estimator.feature_names_in_, np.ndarray)
    assert estimator.feature_names_in_.dtype == object
    assert_array_equal(estimator.feature_names_in_, names)

    # Only check sklearn estimators for feature_names_in_ in docstring
    module_name = estimator_orig.__module__
    if (
        module_name.startswith("sklearn.")
        and not ("test_" in module_name or module_name.endswith("_testing"))
        and ("feature_names_in_" not in (estimator_orig.__doc__))
    ):
        raise ValueError(
            f"Estimator {name} does not document its feature_names_in_ attribute"
        )

    check_methods = []
    for method in (
        "predict",
        "transform",
        "decision_function",
        "predict_proba",
        "score",
        "score_samples",
        "predict_log_proba",
    ):
        if not hasattr(estimator, method):
            continue

        callable_method = getattr(estimator, method)
        if method == "score":
            callable_method = partial(callable_method, y=y)
        check_methods.append((method, callable_method))

    for _, method in check_methods:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "error",
                message="X does not have valid feature names",
                category=UserWarning,
                module="sklearn",
            )
            method(X)  # works without UserWarning for valid features

    invalid_names = [
        (names[::-1], "Feature names must be in the same order as they were in fit."),
        (
            [f"another_prefix_{i}" for i in range(n_features)],
            (
                "Feature names unseen at fit time:\n- another_prefix_0\n-"
                " another_prefix_1\n"
            ),
        ),
        (
            names[:3],
            f"Feature names seen at fit time, yet now missing:\n- {min(names[3:])}\n",
        ),
    ]
    params = {
        key: value
        for key, value in estimator.get_params().items()
        if "early_stopping" in key
    }
    early_stopping_enabled = any(value is True for value in params.values())

    for invalid_name, additional_message in invalid_names:
        X_bad = pd.DataFrame(X, columns=invalid_name, copy=False)

        expected_msg = re.escape(
            "The feature names should match those that were passed during fit.\n"
            f"{additional_message}"
        )
        for name, method in check_methods:
            with raises(
                ValueError, match=expected_msg, err_msg=f"{name} did not raise"
            ):
                method(X_bad)

        # partial_fit checks on second call
        # Do not call partial fit if early_stopping is on
        if not hasattr(estimator, "partial_fit") or early_stopping_enabled:
            continue

        estimator = clone(estimator_orig)
        if is_classifier(estimator):
            classes = np.unique(y)
            estimator.partial_fit(X, y, classes=classes)
        else:
            estimator.partial_fit(X, y)

        with raises(ValueError, match=expected_msg):
            estimator.partial_fit(X_bad, y)


def check_transformer_get_feature_names_out(name, transformer_orig):
    tags = get_tags(transformer_orig)
    if not tags.input_tags.two_d_array or tags.no_validation:
        return

    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        n_features=2,
        cluster_std=0.1,
    )
    X = StandardScaler().fit_transform(X)

    transformer = clone(transformer_orig)
    X = _enforce_estimator_tags_X(transformer, X)

    n_features = X.shape[1]
    set_random_state(transformer)

    y_ = y
    if name in CROSS_DECOMPOSITION:
        y_ = np.c_[np.asarray(y), np.asarray(y)]
        y_[::2, 1] *= 2

    X_transform = transformer.fit_transform(X, y=y_)
    input_features = [f"feature{i}" for i in range(n_features)]

    # input_features names is not the same length as n_features_in_
    with raises(ValueError, match="input_features should have length equal"):
        transformer.get_feature_names_out(input_features[::2])

    feature_names_out = transformer.get_feature_names_out(input_features)
    assert feature_names_out is not None
    assert isinstance(feature_names_out, np.ndarray)
    assert feature_names_out.dtype == object
    assert all(isinstance(name, str) for name in feature_names_out)

    if isinstance(X_transform, tuple):
        n_features_out = X_transform[0].shape[1]
    else:
        n_features_out = X_transform.shape[1]

    assert len(feature_names_out) == n_features_out, (
        f"Expected {n_features_out} feature names, got {len(feature_names_out)}"
    )


def check_transformer_get_feature_names_out_pandas(name, transformer_orig):
    try:
        import pandas as pd
    except ImportError:
        raise SkipTest(
            "pandas is not installed: not checking column name consistency for pandas"
        )

    tags = get_tags(transformer_orig)
    if not tags.input_tags.two_d_array or tags.no_validation:
        return

    X, y = make_blobs(
        n_samples=30,
        centers=[[0, 0, 0], [1, 1, 1]],
        random_state=0,
        n_features=2,
        cluster_std=0.1,
    )
    X = StandardScaler().fit_transform(X)

    transformer = clone(transformer_orig)
    X = _enforce_estimator_tags_X(transformer, X)

    n_features = X.shape[1]
    set_random_state(transformer)

    y_ = y
    if name in CROSS_DECOMPOSITION:
        y_ = np.c_[np.asarray(y), np.asarray(y)]
        y_[::2, 1] *= 2

    feature_names_in = [f"col{i}" for i in range(n_features)]
    df = pd.DataFrame(X, columns=feature_names_in, copy=False)
    X_transform = transformer.fit_transform(df, y=y_)

    # error is raised when `input_features` do not match feature_names_in
    invalid_feature_names = [f"bad{i}" for i in range(n_features)]
    with raises(ValueError, match="input_features is not equal to feature_names_in_"):
        transformer.get_feature_names_out(invalid_feature_names)

    feature_names_out_default = transformer.get_feature_names_out()
    feature_names_in_explicit_names = transformer.get_feature_names_out(
        feature_names_in
    )
    assert_array_equal(feature_names_out_default, feature_names_in_explicit_names)

    if isinstance(X_transform, tuple):
        n_features_out = X_transform[0].shape[1]
    else:
        n_features_out = X_transform.shape[1]

    assert len(feature_names_out_default) == n_features_out, (
        f"Expected {n_features_out} feature names, got {len(feature_names_out_default)}"
    )


def check_param_validation(name, estimator_orig):
    # Check that an informative error is raised when the value of a constructor
    # parameter does not have an appropriate type or value.
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(20, 5))
    y = rng.randint(0, 2, size=20)
    y = _enforce_estimator_tags_y(estimator_orig, y)
    tags = get_tags(estimator_orig)

    estimator_params = estimator_orig.get_params(deep=False).keys()

    # check that there is a constraint for each parameter
    if estimator_params:
        validation_params = estimator_orig._parameter_constraints.keys()
        unexpected_params = set(validation_params) - set(estimator_params)
        missing_params = set(estimator_params) - set(validation_params)
        err_msg = (
            f"Mismatch between _parameter_constraints and the parameters of {name}."
            f"\nConsider the unexpected parameters {unexpected_params} and expected but"
            f" missing parameters {missing_params}"
        )
        assert validation_params == estimator_params, err_msg

    # this object does not have a valid type for sure for all params
    param_with_bad_type = type("BadType", (), {})()

    fit_methods = ["fit", "partial_fit", "fit_transform", "fit_predict"]

    for param_name in estimator_params:
        constraints = estimator_orig._parameter_constraints[param_name]

        if constraints == "no_validation":
            # This parameter is not validated
            continue

        # Mixing an interval of reals and an interval of integers must be avoided.
        if any(
            isinstance(constraint, Interval) and constraint.type == Integral
            for constraint in constraints
        ) and any(
            isinstance(constraint, Interval) and constraint.type == Real
            for constraint in constraints
        ):
            raise ValueError(
                f"The constraint for parameter {param_name} of {name} can't have a mix"
                " of intervals of Integral and Real types. Use the type RealNotInt"
                " instead of Real."
            )

        match = rf"The '{param_name}' parameter of {name} must be .* Got .* instead."
        err_msg = (
            f"{name} does not raise an informative error message when the "
            f"parameter {param_name} does not have a valid type or value."
        )

        estimator = clone(estimator_orig)

        # First, check that the error is raised if param doesn't match any valid type.
        estimator.set_params(**{param_name: param_with_bad_type})

        for method in fit_methods:
            if not hasattr(estimator, method):
                # the method is not accessible with the current set of parameters
                continue

            err_msg = (
                f"{name} does not raise an informative error message when the parameter"
                f" {param_name} does not have a valid type. If any Python type is"
                " valid, the constraint should be 'no_validation'."
            )

            with raises(InvalidParameterError, match=match, err_msg=err_msg):
                if tags.target_tags.one_d_labels or tags.target_tags.two_d_labels:
                    # The estimator is a label transformer and take only `y`
                    getattr(estimator, method)(y)
                else:
                    getattr(estimator, method)(X, y)

        # Then, for constraints that are more than a type constraint, check that the
        # error is raised if param does match a valid type but does not match any valid
        # value for this type.
        constraints = [make_constraint(constraint) for constraint in constraints]

        for constraint in constraints:
            try:
                bad_value = generate_invalid_param_val(constraint)
            except NotImplementedError:
                continue

            estimator.set_params(**{param_name: bad_value})

            for method in fit_methods:
                if not hasattr(estimator, method):
                    # the method is not accessible with the current set of parameters
                    continue

                err_msg = (
                    f"{name} does not raise an informative error message when the "
                    f"parameter {param_name} does not have a valid value.\n"
                    "Constraints should be disjoint. For instance "
                    "[StrOptions({'a_string'}), str] is not a acceptable set of "
                    "constraint because generating an invalid string for the first "
                    "constraint will always produce a valid string for the second "
                    "constraint."
                )

                with raises(InvalidParameterError, match=match, err_msg=err_msg):
                    if tags.target_tags.one_d_labels or tags.target_tags.two_d_labels:
                        # The estimator is a label transformer and take only `y`
                        getattr(estimator, method)(y)
                    else:
                        getattr(estimator, method)(X, y)


def check_set_output_transform(name, transformer_orig):
    # Check transformer.set_output with the default configuration does not
    # change the transform output.
    tags = get_tags(transformer_orig)
    if not tags.input_tags.two_d_array or tags.no_validation:
        return

    rng = np.random.RandomState(0)
    transformer = clone(transformer_orig)

    X = rng.uniform(size=(20, 5))
    X = _enforce_estimator_tags_X(transformer_orig, X)
    y = rng.randint(0, 2, size=20)
    y = _enforce_estimator_tags_y(transformer_orig, y)
    set_random_state(transformer)

    def fit_then_transform(est):
        if name in CROSS_DECOMPOSITION:
            return est.fit(X, y).transform(X, y)
        return est.fit(X, y).transform(X)

    def fit_transform(est):
        return est.fit_transform(X, y)

    transform_methods = {
        "transform": fit_then_transform,
        "fit_transform": fit_transform,
    }
    for name, transform_method in transform_methods.items():
        transformer = clone(transformer)
        if not hasattr(transformer, name):
            continue
        X_trans_no_setting = transform_method(transformer)

        # Auto wrapping only wraps the first array
        if name in CROSS_DECOMPOSITION:
            X_trans_no_setting = X_trans_no_setting[0]

        transformer.set_output(transform="default")
        X_trans_default = transform_method(transformer)

        if name in CROSS_DECOMPOSITION:
            X_trans_default = X_trans_default[0]

        # Default and no setting -> returns the same transformation
        assert_allclose_dense_sparse(X_trans_no_setting, X_trans_default)


def _output_from_fit_transform(transformer, name, X, df, y):
    """Generate output to test `set_output` for different configuration:

    - calling either `fit.transform` or `fit_transform`;
    - passing either a dataframe or a numpy array to fit;
    - passing either a dataframe or a numpy array to transform.
    """
    outputs = {}

    # fit then transform case:
    cases = [
        ("fit.transform/df/df", df, df),
        ("fit.transform/df/array", df, X),
        ("fit.transform/array/df", X, df),
        ("fit.transform/array/array", X, X),
    ]
    if all(hasattr(transformer, meth) for meth in ["fit", "transform"]):
        for (
            case,
            data_fit,
            data_transform,
        ) in cases:
            transformer.fit(data_fit, y)
            if name in CROSS_DECOMPOSITION:
                X_trans, _ = transformer.transform(data_transform, y)
            else:
                X_trans = transformer.transform(data_transform)
            outputs[case] = (X_trans, transformer.get_feature_names_out())

    # fit_transform case:
    cases = [
        ("fit_transform/df", df),
        ("fit_transform/array", X),
    ]
    if hasattr(transformer, "fit_transform"):
        for case, data in cases:
            if name in CROSS_DECOMPOSITION:
                X_trans, _ = transformer.fit_transform(data, y)
            else:
                X_trans = transformer.fit_transform(data, y)
            outputs[case] = (X_trans, transformer.get_feature_names_out())

    return outputs


def _check_generated_dataframe(
    name,
    case,
    index,
    outputs_default,
    outputs_dataframe_lib,
    is_supported_dataframe,
    create_dataframe,
    assert_frame_equal,
):
    """Check if the generated DataFrame by the transformer is valid.

    The DataFrame implementation is specified through the parameters of this function.

    Parameters
    ----------
    name : str
        The name of the transformer.
    case : str
        A single case from the cases generated by `_output_from_fit_transform`.
    index : index or None
        The index of the DataFrame. `None` if the library does not implement a DataFrame
        with an index.
    outputs_default : tuple
        A tuple containing the output data and feature names for the default output.
    outputs_dataframe_lib : tuple
        A tuple containing the output data and feature names for the pandas case.
    is_supported_dataframe : callable
        A callable that takes a DataFrame instance as input and return whether or
        E.g. `lambda X: isintance(X, pd.DataFrame)`.
    create_dataframe : callable
        A callable taking as parameters `data`, `columns`, and `index` and returns
        a callable. Be aware that `index` can be ignored. For example, polars dataframes
        would ignore the idnex.
    assert_frame_equal : callable
        A callable taking 2 dataframes to compare if they are equal.
    """
    X_trans, feature_names_default = outputs_default
    df_trans, feature_names_dataframe_lib = outputs_dataframe_lib

    assert is_supported_dataframe(df_trans)
    # We always rely on the output of `get_feature_names_out` of the
    # transformer used to generate the dataframe as a ground-truth of the
    # columns.
    # If a dataframe is passed into transform, then the output should have the same
    # index
    expected_index = index if case.endswith("df") else None
    expected_dataframe = create_dataframe(
        X_trans, columns=feature_names_dataframe_lib, index=expected_index
    )

    try:
        assert_frame_equal(df_trans, expected_dataframe)
    except AssertionError as e:
        raise AssertionError(
            f"{name} does not generate a valid dataframe in the {case} "
            "case. The generated dataframe is not equal to the expected "
            f"dataframe. The error message is: {e}"
        ) from e


def _check_set_output_transform_dataframe(
    name,
    transformer_orig,
    *,
    dataframe_lib,
    is_supported_dataframe,
    create_dataframe,
    assert_frame_equal,
    context,
):
    """Check that a transformer can output a DataFrame when requested.

    The DataFrame implementation is specified through the parameters of this function.

    Parameters
    ----------
    name : str
        The name of the transformer.
    transformer_orig : estimator
        The original transformer instance.
    dataframe_lib : str
        The name of the library implementing the DataFrame.
    is_supported_dataframe : callable
        A callable that takes a DataFrame instance as input and returns whether or
        not it is supported by the dataframe library.
        E.g. `lambda X: isintance(X, pd.DataFrame)`.
    create_dataframe : callable
        A callable taking as parameters `data`, `columns`, and `index` and returns
        a callable. Be aware that `index` can be ignored. For example, polars dataframes
        will ignore the index.
    assert_frame_equal : callable
        A callable taking 2 dataframes to compare if they are equal.
    context : {"local", "global"}
        Whether to use a local context by setting `set_output(...)` on the transformer
        or a global context by using the `with config_context(...)`
    """
    # Check transformer.set_output configures the output of transform="pandas".
    tags = get_tags(transformer_orig)
    if not tags.input_tags.two_d_array or tags.no_validation:
        return

    rng = np.random.RandomState(0)
    transformer = clone(transformer_orig)

    X = rng.uniform(size=(20, 5))
    X = _enforce_estimator_tags_X(transformer_orig, X)
    y = rng.randint(0, 2, size=20)
    y = _enforce_estimator_tags_y(transformer_orig, y)
    set_random_state(transformer)

    feature_names_in = [f"col{i}" for i in range(X.shape[1])]
    index = [f"index{i}" for i in range(X.shape[0])]
    df = create_dataframe(X, columns=feature_names_in, index=index)

    transformer_default = clone(transformer).set_output(transform="default")
    outputs_default = _output_from_fit_transform(transformer_default, name, X, df, y)

    if context == "local":
        transformer_df = clone(transformer).set_output(transform=dataframe_lib)
        context_to_use = nullcontext()
    else:  # global
        transformer_df = clone(transformer)
        context_to_use = config_context(transform_output=dataframe_lib)

    try:
        with context_to_use:
            outputs_df = _output_from_fit_transform(transformer_df, name, X, df, y)
    except ValueError as e:
        # transformer does not support sparse data
        capitalized_lib = dataframe_lib.capitalize()
        error_message = str(e)
        assert (
            f"{capitalized_lib} output does not support sparse data." in error_message
            or "The transformer outputs a scipy sparse matrix." in error_message
        ), e
        return

    for case in outputs_default:
        _check_generated_dataframe(
            name,
            case,
            index,
            outputs_default[case],
            outputs_df[case],
            is_supported_dataframe,
            create_dataframe,
            assert_frame_equal,
        )


def _check_set_output_transform_pandas_context(name, transformer_orig, context):
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover
        raise SkipTest("pandas is not installed: not checking set output")

    _check_set_output_transform_dataframe(
        name,
        transformer_orig,
        dataframe_lib="pandas",
        is_supported_dataframe=lambda X: isinstance(X, pd.DataFrame),
        create_dataframe=lambda X, columns, index: pd.DataFrame(
            X, columns=columns, copy=False, index=index
        ),
        assert_frame_equal=pd.testing.assert_frame_equal,
        context=context,
    )


def check_set_output_transform_pandas(name, transformer_orig):
    _check_set_output_transform_pandas_context(name, transformer_orig, "local")


def check_global_output_transform_pandas(name, transformer_orig):
    _check_set_output_transform_pandas_context(name, transformer_orig, "global")


def _check_set_output_transform_polars_context(name, transformer_orig, context):
    try:
        import polars as pl
        from polars.testing import assert_frame_equal
    except ImportError:  # pragma: no cover
        raise SkipTest("polars is not installed: not checking set output")

    def create_dataframe(X, columns, index):
        if isinstance(columns, np.ndarray):
            columns = columns.tolist()

        return pl.DataFrame(X, schema=columns, orient="row")

    _check_set_output_transform_dataframe(
        name,
        transformer_orig,
        dataframe_lib="polars",
        is_supported_dataframe=lambda X: isinstance(X, pl.DataFrame),
        create_dataframe=create_dataframe,
        assert_frame_equal=assert_frame_equal,
        context=context,
    )


def check_set_output_transform_polars(name, transformer_orig):
    _check_set_output_transform_polars_context(name, transformer_orig, "local")


def check_global_set_output_transform_polars(name, transformer_orig):
    _check_set_output_transform_polars_context(name, transformer_orig, "global")


@ignore_warnings(category=FutureWarning)
def check_inplace_ensure_writeable(name, estimator_orig):
    """Check that estimators able to do inplace operations can work on read-only
    input data even if a copy is not explicitly requested by the user.

    Make sure that a copy is made and consequently that the input array and its
    writeability are not modified by the estimator.
    """
    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)

    n_samples = 100

    X, _ = make_blobs(n_samples=n_samples, n_features=3, random_state=rng)
    X = _enforce_estimator_tags_X(estimator, X)

    # These estimators can only work inplace with fortran ordered input
    if name in ("Lasso", "ElasticNet", "MultiTaskElasticNet", "MultiTaskLasso"):
        X = np.asfortranarray(X)

    # Add a missing value for imputers so that transform has to do something
    if hasattr(estimator, "missing_values"):
        X[0, 0] = np.nan

    if is_regressor(estimator):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = _enforce_estimator_tags_y(estimator, y)

    X_copy = X.copy()

    # Make X read-only
    X.setflags(write=False)

    estimator.fit(X, y)

    if hasattr(estimator, "transform"):
        estimator.transform(X)

    assert not X.flags.writeable
    assert_allclose(X, X_copy)


def check_do_not_raise_errors_in_init_or_set_params(name, estimator_orig):
    """Check that init or set_param does not raise errors."""
    Estimator = type(estimator_orig)
    params = signature(Estimator).parameters

    smoke_test_values = [-1, 3.0, "helloworld", np.array([1.0, 4.0]), [1], {}, []]
    for value in smoke_test_values:
        new_params = {key: value for key in params}

        # Does not raise
        est = Estimator(**new_params)

        # Also do does not raise
        est.set_params(**new_params)


def check_classifier_not_supporting_multiclass(name, estimator_orig):
    """Check that if the classifier has tags.classifier_tags.multi_class=False,
    then it should raise a ValueError when calling fit with a multiclass dataset.

    This test is not yielded if the tag is not False.
    """
    estimator = clone(estimator_orig)
    set_random_state(estimator)

    X, y = make_classification(
        n_samples=100,
        n_classes=3,
        n_informative=3,
        n_clusters_per_class=1,
        random_state=0,
    )
    err_msg = """\
        The estimator tag `tags.classifier_tags.multi_class` is False for {name}
        which means it does not support multiclass classification. However, it does
        not raise the right `ValueError` when calling fit with a multiclass dataset,
        including the error message 'Only binary classification is supported.' This
        can be achieved by the following pattern:

        y_type = type_of_target(y, input_name='y', raise_unknown=True)
        if y_type != 'binary':
            raise ValueError(
                'Only binary classification is supported. The type of the target '
                f'is {{y_type}}.'
        )
    """.format(name=name)
    err_msg = textwrap.dedent(err_msg)

    with raises(
        ValueError, match="Only binary classification is supported.", err_msg=err_msg
    ):
        estimator.fit(X, y)
