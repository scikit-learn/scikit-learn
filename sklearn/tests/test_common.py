"""
General tests for all estimators in sklearn.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import os
import pkgutil
import re
import warnings
from functools import partial
from inspect import isgenerator
from itertools import chain

import pytest
from scipy.linalg import LinAlgWarning

import sklearn
from sklearn.base import BaseEstimator
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning

# make it possible to discover experimental estimators when calling `all_estimators`
from sklearn.experimental import (
    enable_halving_search_cv,  # noqa: F401
    enable_iterative_imputer,  # noqa: F401
)
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import FeatureUnion, make_pipeline
from sklearn.preprocessing import (
    FunctionTransformer,
    MinMaxScaler,
    OneHotEncoder,
    StandardScaler,
)
from sklearn.utils import all_estimators
from sklearn.utils._test_common.instance_generator import (
    _get_check_estimator_ids,
    _get_expected_failed_checks,
    _tested_estimators,
)
from sklearn.utils._testing import (
    SkipTest,
    ignore_warnings,
)
from sklearn.utils.estimator_checks import (
    check_dataframe_column_names_consistency,
    check_estimator,
    check_get_feature_names_out_error,
    check_global_output_transform_pandas,
    check_global_set_output_transform_polars,
    check_inplace_ensure_writeable,
    check_param_validation,
    check_set_output_transform,
    check_set_output_transform_pandas,
    check_set_output_transform_polars,
    check_transformer_get_feature_names_out,
    check_transformer_get_feature_names_out_pandas,
    parametrize_with_checks,
)


def test_all_estimator_no_base_class():
    # test that all_estimators doesn't find abstract classes.
    for name, Estimator in all_estimators():
        msg = (
            "Base estimators such as {0} should not be included in all_estimators"
        ).format(name)
        assert not name.lower().startswith("base"), msg


def _sample_func(x, y=1):
    pass


class CallableEstimator(BaseEstimator):
    """Dummy development stub for an estimator.

    This is to make sure a callable estimator passes common tests.
    """

    def __call__(self):
        pass  # pragma: nocover


@pytest.mark.parametrize(
    "val, expected",
    [
        (partial(_sample_func, y=1), "_sample_func(y=1)"),
        (_sample_func, "_sample_func"),
        (partial(_sample_func, "world"), "_sample_func"),
        (LogisticRegression(C=2.0), "LogisticRegression(C=2.0)"),
        (
            LogisticRegression(
                random_state=1,
                solver="newton-cg",
                class_weight="balanced",
                warm_start=True,
            ),
            (
                "LogisticRegression(class_weight='balanced',random_state=1,"
                "solver='newton-cg',warm_start=True)"
            ),
        ),
        (CallableEstimator(), "CallableEstimator()"),
    ],
)
def test_get_check_estimator_ids(val, expected):
    assert _get_check_estimator_ids(val) == expected


@parametrize_with_checks(
    list(_tested_estimators()), expected_failed_checks=_get_expected_failed_checks
)
def test_estimators(estimator, check, request):
    # Common tests for estimator instances
    with ignore_warnings(
        category=(FutureWarning, ConvergenceWarning, UserWarning, LinAlgWarning)
    ):
        check(estimator)


# TODO(1.8): remove test when generate_only is removed
def test_check_estimator_generate_only_deprecation():
    """Check that check_estimator with generate_only=True raises a deprecation
    warning."""
    with pytest.warns(FutureWarning, match="`generate_only` is deprecated in 1.6"):
        all_instance_gen_checks = check_estimator(
            LogisticRegression(), generate_only=True
        )
    assert isgenerator(all_instance_gen_checks)


@pytest.mark.filterwarnings(
    "ignore:Since version 1.0, it is not needed to import "
    "enable_hist_gradient_boosting anymore"
)
def test_import_all_consistency():
    sklearn_path = [os.path.dirname(sklearn.__file__)]
    # Smoke test to check that any name in a __all__ list is actually defined
    # in the namespace of the module or package.
    pkgs = pkgutil.walk_packages(
        path=sklearn_path, prefix="sklearn.", onerror=lambda _: None
    )
    submods = [modname for _, modname, _ in pkgs]
    for modname in submods + ["sklearn"]:
        if ".tests." in modname or "sklearn.externals" in modname:
            continue
        # Avoid test suite depending on build dependencies, for example Cython
        if "sklearn._build_utils" in modname:
            continue
        package = __import__(modname, fromlist="dummy")
        for name in getattr(package, "__all__", ()):
            assert hasattr(package, name), "Module '{0}' has no attribute '{1}'".format(
                modname, name
            )


def test_root_import_all_completeness():
    sklearn_path = [os.path.dirname(sklearn.__file__)]
    EXCEPTIONS = ("utils", "tests", "base", "conftest")
    for _, modname, _ in pkgutil.walk_packages(
        path=sklearn_path, onerror=lambda _: None
    ):
        if "." in modname or modname.startswith("_") or modname in EXCEPTIONS:
            continue
        assert modname in sklearn.__all__


def test_all_tests_are_importable():
    # Ensure that for each contentful subpackage, there is a test directory
    # within it that is also a subpackage (i.e. a directory with __init__.py)

    HAS_TESTS_EXCEPTIONS = re.compile(
        r"""(?x)
                                      \.externals(\.|$)|
                                      \.tests(\.|$)|
                                      \._
                                      """
    )
    resource_modules = {
        "sklearn.datasets.data",
        "sklearn.datasets.descr",
        "sklearn.datasets.images",
    }
    sklearn_path = [os.path.dirname(sklearn.__file__)]
    lookup = {
        name: ispkg
        for _, name, ispkg in pkgutil.walk_packages(sklearn_path, prefix="sklearn.")
    }
    missing_tests = [
        name
        for name, ispkg in lookup.items()
        if ispkg
        and name not in resource_modules
        and not HAS_TESTS_EXCEPTIONS.search(name)
        and name + ".tests" not in lookup
    ]
    assert missing_tests == [], (
        "{0} do not have `tests` subpackages. "
        "Perhaps they require "
        "__init__.py or a meson.build "
        "in the parent "
        "directory".format(missing_tests)
    )


def test_class_support_removed():
    # Make sure passing classes to check_estimator or parametrize_with_checks
    # raises an error

    msg = "Passing a class was deprecated.* isn't supported anymore"
    with pytest.raises(TypeError, match=msg):
        check_estimator(LogisticRegression)

    with pytest.raises(TypeError, match=msg):
        parametrize_with_checks([LogisticRegression])


def _estimators_that_predict_in_fit():
    for estimator in _tested_estimators():
        est_params = set(estimator.get_params())
        if "oob_score" in est_params:
            yield estimator.set_params(oob_score=True, bootstrap=True)
        elif "early_stopping" in est_params:
            est = estimator.set_params(early_stopping=True, n_iter_no_change=1)
            if est.__class__.__name__ in {"MLPClassifier", "MLPRegressor"}:
                # TODO: FIX MLP to not check validation set during MLP
                yield pytest.param(
                    est, marks=pytest.mark.xfail(msg="MLP still validates in fit")
                )
            else:
                yield est
        elif "n_iter_no_change" in est_params:
            yield estimator.set_params(n_iter_no_change=1)


# NOTE: When running `check_dataframe_column_names_consistency` on a meta-estimator that
# delegates validation to a base estimator, the check is testing that the base estimator
# is checking for column name consistency.
column_name_estimators = list(
    chain(
        _tested_estimators(),
        [make_pipeline(LogisticRegression(C=1))],
        _estimators_that_predict_in_fit(),
    )
)


@pytest.mark.parametrize(
    "estimator", column_name_estimators, ids=_get_check_estimator_ids
)
def test_pandas_column_name_consistency(estimator):
    if isinstance(estimator, ColumnTransformer):
        pytest.skip("ColumnTransformer is not tested here")
    if "check_dataframe_column_names_consistency" in _get_expected_failed_checks(
        estimator
    ):
        pytest.skip(
            "Estimator does not support check_dataframe_column_names_consistency"
        )
    with ignore_warnings(category=(FutureWarning)):
        with warnings.catch_warnings(record=True) as record:
            check_dataframe_column_names_consistency(
                estimator.__class__.__name__, estimator
            )
        for warning in record:
            assert "was fitted without feature names" not in str(warning.message)


# TODO: As more modules support get_feature_names_out they should be removed
# from this list to be tested
GET_FEATURES_OUT_MODULES_TO_IGNORE = [
    "ensemble",
    "kernel_approximation",
]


def _include_in_get_feature_names_out_check(transformer):
    if hasattr(transformer, "get_feature_names_out"):
        return True
    module = transformer.__module__.split(".")[1]
    return module not in GET_FEATURES_OUT_MODULES_TO_IGNORE


GET_FEATURES_OUT_ESTIMATORS = [
    est
    for est in _tested_estimators("transformer")
    if _include_in_get_feature_names_out_check(est)
]


@pytest.mark.parametrize(
    "transformer", GET_FEATURES_OUT_ESTIMATORS, ids=_get_check_estimator_ids
)
def test_transformers_get_feature_names_out(transformer):
    with ignore_warnings(category=(FutureWarning)):
        check_transformer_get_feature_names_out(
            transformer.__class__.__name__, transformer
        )
        check_transformer_get_feature_names_out_pandas(
            transformer.__class__.__name__, transformer
        )


ESTIMATORS_WITH_GET_FEATURE_NAMES_OUT = [
    est for est in _tested_estimators() if hasattr(est, "get_feature_names_out")
]


@pytest.mark.parametrize(
    "estimator", ESTIMATORS_WITH_GET_FEATURE_NAMES_OUT, ids=_get_check_estimator_ids
)
def test_estimators_get_feature_names_out_error(estimator):
    estimator_name = estimator.__class__.__name__
    check_get_feature_names_out_error(estimator_name, estimator)


@pytest.mark.parametrize(
    "estimator", list(_tested_estimators()), ids=_get_check_estimator_ids
)
def test_check_param_validation(estimator):
    if isinstance(estimator, FeatureUnion):
        pytest.skip("FeatureUnion is not tested here")
    name = estimator.__class__.__name__
    check_param_validation(name, estimator)


SET_OUTPUT_ESTIMATORS = list(
    chain(
        _tested_estimators("transformer"),
        [
            make_pipeline(StandardScaler(), MinMaxScaler()),
            OneHotEncoder(sparse_output=False),
            FunctionTransformer(feature_names_out="one-to-one"),
        ],
    )
)


@pytest.mark.parametrize(
    "estimator", SET_OUTPUT_ESTIMATORS, ids=_get_check_estimator_ids
)
def test_set_output_transform(estimator):
    name = estimator.__class__.__name__
    if not hasattr(estimator, "set_output"):
        pytest.skip(
            f"Skipping check_set_output_transform for {name}: Does not support"
            " set_output API"
        )
    with ignore_warnings(category=(FutureWarning)):
        check_set_output_transform(estimator.__class__.__name__, estimator)


@pytest.mark.parametrize(
    "estimator", SET_OUTPUT_ESTIMATORS, ids=_get_check_estimator_ids
)
@pytest.mark.parametrize(
    "check_func",
    [
        check_set_output_transform_pandas,
        check_global_output_transform_pandas,
        check_set_output_transform_polars,
        check_global_set_output_transform_polars,
    ],
)
def test_set_output_transform_configured(estimator, check_func):
    name = estimator.__class__.__name__
    if not hasattr(estimator, "set_output"):
        pytest.skip(
            f"Skipping {check_func.__name__} for {name}: Does not support"
            " set_output API yet"
        )
    with ignore_warnings(category=(FutureWarning)):
        check_func(estimator.__class__.__name__, estimator)


@pytest.mark.parametrize(
    "estimator", _tested_estimators(), ids=_get_check_estimator_ids
)
def test_check_inplace_ensure_writeable(estimator):
    name = estimator.__class__.__name__

    if hasattr(estimator, "copy"):
        estimator.set_params(copy=False)
    elif hasattr(estimator, "copy_X"):
        estimator.set_params(copy_X=False)
    else:
        raise SkipTest(f"{name} doesn't require writeable input.")

    # The following estimators can work inplace only with certain settings
    if name == "HDBSCAN":
        estimator.set_params(metric="precomputed", algorithm="brute")

    if name == "PCA":
        estimator.set_params(svd_solver="full")

    if name == "KernelPCA":
        estimator.set_params(kernel="precomputed")

    check_inplace_ensure_writeable(name, estimator)
