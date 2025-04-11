"""
Unit tests for the SimpleMKL algorithm in the MKL module.
"""

from itertools import cycle, islice, permutations

import numpy as np
import pytest

from sklearn import mkl
from sklearn.datasets import make_classification, make_regression

# Sample data for classification and regression
X_class, y_class = make_classification(n_samples=50, n_features=10, random_state=42)
X_reg, y_reg = make_regression(n_samples=50, n_features=10, random_state=42)
X_oc, _ = make_classification(n_samples=50, n_features=10, n_classes=1, random_state=42)

# Parameters for the tests
six_kernels_params = {
    "kernels": ["rbf", "poly"],
    "kernels_scopes": ["all", "all"],
    "kernels_param_grids": [
        {"gamma": [0.1, 1, 10]},
        {"degree": [2, 3, 4]},
    ],
}


def test_mklc_kernel_order_independence():
    # Test that MKLC is independent of kernel order.
    # Train MKLC with different kernel orders
    obj_vals = []
    for perm in permutations(["rbf", "linear", "poly", "sigmoid"]):
        clf = mkl.MKLC(kernels=list(perm), algo="simple", C=1.0, random_state=42)
        clf.fit(X_class, y_class)
        obj_vals.append(clf._svm.objective_val_)

    # Check that all objective values are equal
    for i in range(len(obj_vals) - 1):
        assert obj_vals[i] == pytest.approx(obj_vals[-1], rel=0.01)


def test_mklr_kernel_order_independence():
    # Test that MKLR is independent of kernel order.
    # Train MKLR with different kernel orders
    obj_vals = []
    for perm in permutations(["rbf", "linear", "poly", "sigmoid"]):
        reg = mkl.MKLR(
            kernels=list(perm), algo="simple", C=1.0, epsilon=0.1, random_state=42
        )
        reg.fit(X_reg, y_reg)
        obj_vals.append(reg._svm.objective_val_)

    # Check that all objective values are equal
    for i in range(len(obj_vals) - 1):
        assert obj_vals[i] == pytest.approx(obj_vals[-1], rel=0.01)


def test_oneclassmkl_kernel_order_independance():
    # Test that OneClassMKL is independent of kernel order.
    # Train OneClassMKL with different kernel orders
    obj_vals = []
    for perm in permutations(["rbf", "linear", "poly", "sigmoid"]):
        oc = mkl.OneClassMKL(kernels=list(perm), algo="simple", nu=0.5, random_state=42)
        oc.fit(X_oc)
        obj_vals.append(oc._svm.objective_val_)

    # Check that all objective values are equal
    for i in range(len(obj_vals) - 1):
        assert obj_vals[i] == pytest.approx(obj_vals[-1], rel=0.01)


def test_mklc_binary_better_than_average():
    # Test that MKLC performs better than algo='average'.
    # Create an MKLC object with six kernels
    clf = mkl.MKLC(**six_kernels_params, C=1.0, random_state=42)

    # Get the objective values for both algorithms
    clf.set_params(algo="simple").fit(X_class, y_class)
    obj_simple = clf._svm.objective_val_

    clf.set_params(algo="average").fit(X_class, y_class)
    obj_average = clf._svm.objective_val_

    # Check that the simple algorithm performs better than the average algorithm
    assert obj_simple < obj_average


def test_mklc_multiclass_better_than_average():
    # Test that MKLC performs better than algo='average' for multiclass.
    # Create a multiclass dataset
    X_multclass, y_multclass = make_classification(
        n_samples=100, n_features=10, n_informative=3, n_classes=3, random_state=42
    )

    # Create an MKLC object with six kernels
    clf = mkl.MKLC(**six_kernels_params, C=1.0, random_state=42)

    # Get the objective values for both algorithms
    clf.set_params(algo="simple").fit(X_multclass, y_multclass)
    obj_simple = clf._svm.objective_val_

    clf.set_params(algo="average").fit(X_multclass, y_multclass)
    obj_average = clf._svm.objective_val_

    # Check that the simple algorithm performs better than the average algorithm
    assert obj_simple < obj_average


def test_mklr_better_than_average():
    # Test that MKLR performs better than algo='average'.
    # Create an MKLR object with six kernels
    reg = mkl.MKLR(**six_kernels_params, C=1.0, epsilon=0.1, random_state=42)

    # Get the objective values for both algorithms
    reg.set_params(algo="simple").fit(X_reg, y_reg)
    obj_simple = reg._svm.objective_val_

    reg.set_params(algo="average").fit(X_reg, y_reg)
    obj_average = reg._svm.objective_val_

    # Check that the simple algorithm performs better than the average algorithm
    assert obj_simple < obj_average


def test_oneclassmkl_better_than_average():
    # Test that OneClassMKL performs better than algo='average'.
    # Create an OneClassMKL object with six kernels
    oc = mkl.OneClassMKL(**six_kernels_params, nu=0.5, random_state=42)

    # Get the objective values for both algorithms
    oc.set_params(algo="simple").fit(X_oc)
    obj_simple = oc._svm.objective_val_

    oc.set_params(algo="average").fit(X_oc)
    obj_average = oc._svm.objective_val_

    # Check that the simple algorithm performs better than the average algorithm
    assert obj_simple < obj_average


def assert_weights_conditions(estimator):
    # Check that the weights are all positive and sum to 1
    assert hasattr(estimator, "weights_")
    assert np.all(estimator.weights_ >= 0)
    assert np.sum(estimator.weights_) == pytest.approx(1.0)


def test_mklc_weights_conditions():
    # Test MKLC weights conditions.
    # Create an MKLC object with six kernels
    clf = mkl.MKLC(
        **six_kernels_params,
        algo="simple",
        C=1.0,
        tol=2e-2,
        random_state=42,
    ).fit(X_class, y_class)

    assert_weights_conditions(clf)


def test_mklr_weights_conditions():
    # Test MKLR weights conditions.
    # Create an MKLR object with six kernels
    reg = mkl.MKLR(
        **six_kernels_params,
        algo="simple",
        C=1.0,
        epsilon=0.1,
        tol=2e-2,
        random_state=42,
    ).fit(X_reg, y_reg)

    assert_weights_conditions(reg)


def test_oneclassmkl_weights_conditions():
    # Test OneClassMKL weights conditions.
    # Create an OneClassMKL object with six kernels
    oc = mkl.OneClassMKL(
        **six_kernels_params,
        algo="simple",
        nu=0.5,
        tol=2e-2,
        random_state=42,
    ).fit(X_oc)

    assert_weights_conditions(oc)


def test_mklc_weights_permutations():
    # Test MKLC weights permutations.
    # Create an MKLC object with six kernels
    clf = mkl.MKLC(
        **six_kernels_params,
        algo="simple",
        C=1.0,
        tol=2e-2,
        random_state=42,
    ).fit(X_class, y_class)

    score_no_perm = clf.score(X_class, y_class)

    # Check that permuting the weights decreases the score
    weights = np.copy(clf.weights_)
    n = len(weights)
    for i in range(clf.n_kernels_):
        clf.weights_ = np.array(list(islice(cycle(weights), i, i + n)))
        score_perm = clf.score(X_class, y_class)
        assert score_no_perm >= score_perm


def test_mklr_weights_permutations():
    # Test MKLR weights permutations.
    # Create an MKLR object with six kernels
    reg = mkl.MKLR(
        **six_kernels_params,
        algo="simple",
        C=1.0,
        epsilon=0.1,
        tol=2e-2,
        random_state=42,
    ).fit(X_reg, y_reg)

    score_no_perm = reg.score(X_reg, y_reg)

    # Check that permuting the weights decreases the score
    weights = np.copy(reg.weights_)
    n = len(weights)
    for i in range(reg.n_kernels_):
        reg.weights_ = np.array(list(islice(cycle(weights), i, i + n)))
        score_perm = reg.score(X_reg, y_reg)
        assert score_no_perm >= score_perm


def test_oneclassmkl_weights_permutations():
    # Test OneClassMKL weights permutations.
    # Create a OneClassMKL object with six kernels
    oc = mkl.OneClassMKL(
        **six_kernels_params,
        algo="simple",
        nu=0.5,
        tol=2e-2,
        random_state=42,
    ).fit(X_oc)

    score_no_perm = oc.decision_function(X_oc).mean()

    # Check that permuting the weights decreases the score
    weights = np.copy(oc.weights_)
    n = len(weights)
    for i in range(oc.n_kernels_):
        oc.weights_ = np.array(list(islice(cycle(weights), i, i + n)))
        score_perm = oc.decision_function(X_oc).mean()
        assert score_no_perm >= score_perm
