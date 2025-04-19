"""
Unit tests for the AverageMKL algorithm in the MKL module.
"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from sklearn import mkl, svm
from sklearn.datasets import make_classification, make_regression
from sklearn.metrics.pairwise import linear_kernel, rbf_kernel

# Sample data for classification and regression
X_class, y_class = make_classification(n_samples=50, n_features=10, random_state=42)
K_linear_class = linear_kernel(X_class)
K_rbf_class = rbf_kernel(X_class)

X_reg, y_reg = make_regression(n_samples=50, n_features=10, random_state=42)
K_linear_reg = linear_kernel(X_reg)
K_rbf_reg = rbf_kernel(X_reg)

X_oc, _ = make_classification(n_samples=50, n_features=10, n_classes=1, random_state=42)
K_linear_oc = linear_kernel(X_oc)
K_rbf_oc = rbf_kernel(X_oc)

# Parameters for the tests
six_kernels_params = {
    "kernels": ["rbf", "poly"],
    "kernels_scopes": ["all", "all"],
    "kernels_param_grids": [
        {"gamma": [0.1, 1, 10]},
        {"degree": [2, 3, 4]},
    ],
}


def test_mklc_algo_average_equivalence_to_svc():
    # Test that MKLC is equivalent to SVC with manually averaged kernels.
    # Compute manually averaged kernel
    K_linear = K_linear_class
    K_rbf = K_rbf_class
    K_avg = (K_linear + K_rbf) / 2

    # Train SVC with manually averaged kernel
    svc = svm.SVC(kernel="precomputed", C=1.0, random_state=42)
    svc.fit(K_avg, y_class)

    # Train MKLC
    clf = mkl.MKLC(kernels=["linear", "rbf"], algo="average", C=1.0, random_state=42)
    clf.fit(X_class, y_class)

    # Same as above but with precomputed kernels
    clf_precomp = mkl.MKLC(algo="average", C=1.0)
    clf_precomp.fit([K_linear, K_rbf], y_class)

    assert_array_equal(clf.classes_, svc.classes_)
    assert_array_equal(clf_precomp.classes_, svc.classes_)

    assert_array_equal(clf.predict(X_class), svc.predict(K_avg))
    assert_array_equal(clf_precomp.predict([K_linear, K_rbf]), svc.predict(K_avg))

    assert_array_equal(clf.decision_function(X_class), svc.decision_function(K_avg))
    assert_array_equal(
        clf_precomp.decision_function([K_linear, K_rbf]), svc.decision_function(K_avg)
    )

    assert clf.score(X_class, y_class) == svc.score(K_avg, y_class)
    assert clf_precomp.score([K_linear, K_rbf], y_class) == svc.score(K_avg, y_class)


def test_mklr_algo_average_equivalence_to_svr():
    # Test that MKLR is equivalent to SVR with manually averaged kernels.
    # Compute manually averaged kernel
    K_linear = K_linear_reg
    K_rbf = K_rbf_reg
    K_avg = (K_linear + K_rbf) / 2

    # Train SVR with manually averaged kernel
    svr = svm.SVR(kernel="precomputed", C=1.0, epsilon=0.1)
    svr.fit(K_avg, y_reg)

    # Train MKLR
    reg = mkl.MKLR(
        kernels=["linear", "rbf"], algo="average", C=1.0, epsilon=0.1, random_state=42
    )
    reg.fit(X_reg, y_reg)

    # Same as above but with precomputed kernels
    reg_precomp = mkl.MKLR(algo="average", C=1.0, epsilon=0.1)
    reg_precomp.fit([K_linear, K_rbf], y_reg)

    assert_array_almost_equal(reg.predict(X_reg), svr.predict(K_avg))
    assert_array_almost_equal(
        reg_precomp.predict([K_linear, K_rbf]), svr.predict(K_avg)
    )

    assert reg.score(X_reg, y_reg) == pytest.approx(svr.score(K_avg, y_reg))
    assert reg_precomp.score([K_linear, K_rbf], y_reg) == pytest.approx(
        svr.score(K_avg, y_reg)
    )


def test_oneclassmkl_algo_average_equivalence_to_oneclasssvm():
    # Test that OneClassMKL is equivalent to OneClassSVM with manually averaged kernels.
    # Compute manually averaged kernel
    K_linear = K_linear_oc
    K_rbf = K_rbf_oc
    K_avg = (K_linear + K_rbf) / 2

    # Train OneClassSVM with manually averaged kernel
    ocsvm = svm.OneClassSVM(kernel="precomputed", nu=0.5)
    ocsvm.fit(K_avg)

    # Train OneClassMKL
    oc = mkl.OneClassMKL(
        kernels=["linear", "rbf"], algo="average", nu=0.5, random_state=42
    )
    oc.fit(X_oc)

    # Same as above but with precomputed kernels
    oc_precomp = mkl.OneClassMKL(algo="average", nu=0.5)
    oc_precomp.fit([K_linear, K_rbf])

    assert_array_equal(oc.predict(X_oc), ocsvm.predict(K_avg))
    assert_array_equal(oc_precomp.predict([K_linear, K_rbf]), ocsvm.predict(K_avg))

    assert_array_equal(oc.decision_function(X_oc), ocsvm.decision_function(K_avg))
    assert_array_equal(
        oc_precomp.decision_function([K_linear, K_rbf]), ocsvm.decision_function(K_avg)
    )

    assert_array_equal(oc.score_samples(X_oc), ocsvm.score_samples(K_avg))
    assert_array_equal(
        oc_precomp.score_samples([K_linear, K_rbf]), ocsvm.score_samples(K_avg)
    )


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
        algo="average",
        C=1.0,
        random_state=42,
    ).fit(X_class, y_class)

    assert_weights_conditions(clf)


def test_mklr_weights_conditions():
    # Test MKLR weights conditions.
    # Create an MKLR object with six kernels
    reg = mkl.MKLR(
        **six_kernels_params,
        algo="average",
        C=1.0,
        epsilon=0.1,
        random_state=42,
    ).fit(X_reg, y_reg)

    assert_weights_conditions(reg)


def test_oneclassmkl_weights_conditions():
    # Test OneClassMKL weights conditions.
    # Create an OneClassMKL object with six kernels
    oc = mkl.OneClassMKL(
        **six_kernels_params,
        algo="average",
        nu=0.5,
        random_state=42,
    ).fit(X_oc)

    assert_weights_conditions(oc)
