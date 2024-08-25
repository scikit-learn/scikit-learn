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
from inspect import isgenerator, signature
from itertools import chain

import numpy as np
import pytest
from scipy.linalg import LinAlgWarning

import sklearn
from sklearn.base import BaseEstimator
from sklearn.cluster import (
    OPTICS,
    AffinityPropagation,
    Birch,
    MeanShift,
    SpectralClustering,
)
from sklearn.datasets import make_blobs
from sklearn.exceptions import ConvergenceWarning, FitFailedWarning

# make it possible to discover experimental estimators when calling `all_estimators`
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model._base import LinearClassifierMixin
from sklearn.manifold import TSNE, Isomap, LocallyLinearEmbedding
from sklearn.neighbors import (
    KNeighborsClassifier,
    KNeighborsRegressor,
    LocalOutlierFactor,
    RadiusNeighborsClassifier,
    RadiusNeighborsRegressor,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import (
    FunctionTransformer,
    MinMaxScaler,
    OneHotEncoder,
    StandardScaler,
)
from sklearn.semi_supervised import LabelPropagation, LabelSpreading
from sklearn.utils import all_estimators
from sklearn.utils._tags import _DEFAULT_TAGS, _safe_tags
from sklearn.utils._test_common.instance_generator import (
    _generate_column_transformer_instances,
    _generate_pipeline,
    _generate_search_cv_instances,
    _get_check_estimator_ids,
    _set_checking_parameters,
    _tested_estimators,
)
from sklearn.utils._testing import (
    SkipTest,
    ignore_warnings,
)
from sklearn.utils.estimator_checks import (
    check_class_weight_balanced_linear_classifier,
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
from sklearn.utils.fixes import _IS_WASM


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


@parametrize_with_checks(list(chain(_tested_estimators(), _generate_pipeline())))
def test_estimators(estimator, check, request):
    # Common tests for estimator instances
    with ignore_warnings(
        category=(FutureWarning, ConvergenceWarning, UserWarning, LinAlgWarning)
    ):
        _set_checking_parameters(estimator)
        check(estimator)


def test_check_estimator_generate_only():
    all_instance_gen_checks = check_estimator(LogisticRegression(), generate_only=True)
    assert isgenerator(all_instance_gen_checks)


def _tested_linear_classifiers():
    classifiers = all_estimators(type_filter="classifier")

    with warnings.catch_warnings(record=True):
        for name, clazz in classifiers:
            required_parameters = getattr(clazz, "_required_parameters", [])
            if len(required_parameters):
                # FIXME
                continue

            if "class_weight" in clazz().get_params().keys() and issubclass(
                clazz, LinearClassifierMixin
            ):
                yield name, clazz


@pytest.mark.parametrize("name, Classifier", _tested_linear_classifiers())
def test_class_weight_balanced_linear_classifiers(name, Classifier):
    check_class_weight_balanced_linear_classifier(name, Classifier)


@pytest.mark.xfail(_IS_WASM, reason="importlib not supported for Pyodide packages")
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
        if ".tests." in modname:
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


@parametrize_with_checks(list(_generate_search_cv_instances()))
def test_search_cv(estimator, check, request):
    # Common tests for SearchCV instances
    # We have a separate test because those meta-estimators can accept a
    # wide range of base estimators (classifiers, regressors, pipelines)
    with ignore_warnings(
        category=(
            FutureWarning,
            ConvergenceWarning,
            UserWarning,
            FitFailedWarning,
        )
    ):
        check(estimator)


@pytest.mark.parametrize(
    "estimator", _tested_estimators(), ids=_get_check_estimator_ids
)
def test_valid_tag_types(estimator):
    """Check that estimator tags are valid."""
    tags = _safe_tags(estimator)

    for name, tag in tags.items():
        correct_tags = type(_DEFAULT_TAGS[name])
        if name == "_xfail_checks":
            # _xfail_checks can be a dictionary
            correct_tags = (correct_tags, dict)
        assert isinstance(tag, correct_tags)


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
    _set_checking_parameters(transformer)

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
    _set_checking_parameters(estimator)
    check_get_feature_names_out_error(estimator_name, estimator)


@pytest.mark.parametrize(
    "Estimator",
    [est for name, est in all_estimators()],
)
def test_estimators_do_not_raise_errors_in_init_or_set_params(Estimator):
    """Check that init or set_param does not raise errors."""
    params = signature(Estimator).parameters

    smoke_test_values = [-1, 3.0, "helloworld", np.array([1.0, 4.0]), [1], {}, []]
    for value in smoke_test_values:
        new_params = {key: value for key in params}

        # Does not raise
        est = Estimator(**new_params)

        # Also do does not raise
        est.set_params(**new_params)


@pytest.mark.parametrize(
    "estimator",
    chain(
        _tested_estimators(),
        _generate_pipeline(),
        _generate_column_transformer_instances(),
        _generate_search_cv_instances(),
    ),
    ids=_get_check_estimator_ids,
)
def test_check_param_validation(estimator):
    name = estimator.__class__.__name__
    _set_checking_parameters(estimator)
    check_param_validation(name, estimator)


@pytest.mark.parametrize(
    "Estimator",
    [
        AffinityPropagation,
        Birch,
        MeanShift,
        KNeighborsClassifier,
        KNeighborsRegressor,
        RadiusNeighborsClassifier,
        RadiusNeighborsRegressor,
        LabelPropagation,
        LabelSpreading,
        OPTICS,
        SpectralClustering,
        LocalOutlierFactor,
        LocallyLinearEmbedding,
        Isomap,
        TSNE,
    ],
)
def test_f_contiguous_array_estimator(Estimator):
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/23988
    # https://github.com/scikit-learn/scikit-learn/issues/24013

    X, _ = make_blobs(n_samples=80, n_features=4, random_state=0)
    X = np.asfortranarray(X)
    y = np.round(X[:, 0])

    est = Estimator()
    est.fit(X, y)

    if hasattr(est, "transform"):
        est.transform(X)

    if hasattr(est, "predict"):
        est.predict(X)


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
    _set_checking_parameters(estimator)
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
    _set_checking_parameters(estimator)
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

    _set_checking_parameters(estimator)

    # The following estimators can work inplace only with certain settings
    if name == "HDBSCAN":
        estimator.set_params(metric="precomputed", algorithm="brute")

    if name == "PCA":
        estimator.set_params(svd_solver="full")

    if name == "KernelPCA":
        estimator.set_params(kernel="precomputed")

    check_inplace_ensure_writeable(name, estimator)
