"""General API checks for estimators."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from copy import deepcopy
from inspect import signature

import joblib
import numpy as np

from sklearn import clone
from sklearn.datasets import make_blobs
from sklearn.exceptions import NotFittedError
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.utils._test_common._common import (
    _enforce_estimator_tags_X,
    _enforce_estimator_tags_y,
    _is_public_parameter,
    _regression_dataset,
)
from sklearn.utils._testing import (
    _get_args,
    create_memmap_backed_data,
    ignore_warnings,
    raises,
    set_random_state,
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


def check_estimators_do_not_raise_errors_in_init_or_set_params(name, estimator_orig):
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
def check_estimators_fit_returns_self(name, estimator_orig, readonly_memmap=False):
    """Check if self is returned when calling fit."""
    X, y = make_blobs(random_state=0, n_samples=21)
    X = _enforce_estimator_tags_X(estimator_orig, X)

    estimator = clone(estimator_orig)
    y = _enforce_estimator_tags_y(estimator, y)

    if readonly_memmap:
        X, y = create_memmap_backed_data([X, y])

    set_random_state(estimator)
    assert estimator.fit(X, y) is estimator


@ignore_warnings
def check_estimators_unfitted(name, estimator_orig):
    """Check that predict raises an exception in an unfitted estimator.

    Unfitted estimators should raise a NotFittedError.
    """
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
            with raises(NotFittedError):
                getattr(estimator, method)(X)
