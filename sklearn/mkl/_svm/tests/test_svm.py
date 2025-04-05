"""
Unit tests for the internal SVM of the MKL module.
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn import svm
from sklearn.datasets import make_classification, make_regression
from sklearn.metrics.pairwise import linear_kernel
from sklearn.mkl import _svm as mkl_svm

# Sample data for classification and regression
X_class, y_class = make_classification(n_samples=50, n_features=10, random_state=42)
X_reg, y_reg = make_regression(n_samples=50, n_features=10, random_state=42)
X_oc, _ = make_classification(n_samples=50, n_features=10, n_classes=1, random_state=42)


def test_mkl_internal_svc_same_results_as_svc():
    # Test that MKL internal SVC gives the same results as sklearn SVC
    svc = svm.SVC(kernel="linear", C=1.0, random_state=42)
    svc.fit(X_class, y_class)

    # Train MKL internal SVC
    mkl_svc = mkl_svm.SVC(kernel="linear", C=1.0, random_state=42)
    mkl_svc.fit(X_class, y_class)

    X_class_pre = linear_kernel(X_class)
    mkl_svc_pre = mkl_svm.SVC(kernel="precomputed", C=1.0, random_state=42)
    mkl_svc_pre.fit(X_class_pre, y_class)

    mkl_svc_call = mkl_svm.SVC(kernel=linear_kernel, C=1.0, random_state=42)
    mkl_svc_call.fit(X_class, y_class)

    assert_array_equal(mkl_svc.classes_, svc.classes_)
    assert_array_equal(mkl_svc_pre.classes_, svc.classes_)
    assert_array_equal(mkl_svc_call.classes_, svc.classes_)

    assert_array_equal(mkl_svc.predict(X_class), svc.predict(X_class))
    assert_array_equal(mkl_svc_pre.predict(X_class_pre), svc.predict(X_class))
    assert_array_equal(mkl_svc_call.predict(X_class), svc.predict(X_class))


def test_mkl_internal_svr_same_results_as_svr():
    # Test that MKL internal SVR gives the same results as sklearn SVR
    svr = svm.SVR(kernel="linear", C=1.0, epsilon=0.1)
    svr.fit(X_reg, y_reg)

    # Train MKL internal SVR
    mkl_svr = mkl_svm.SVR(kernel="linear", C=1.0, epsilon=0.1)
    mkl_svr.fit(X_reg, y_reg)

    X_reg_pre = linear_kernel(X_reg)
    mkl_svr_pre = mkl_svm.SVR(kernel="precomputed", C=1.0, epsilon=0.1)
    mkl_svr_pre.fit(X_reg_pre, y_reg)

    mkl_svr_call = mkl_svm.SVR(kernel=linear_kernel, C=1.0, epsilon=0.1)
    mkl_svr_call.fit(X_reg, y_reg)

    assert_array_equal(mkl_svr.predict(X_reg), svr.predict(X_reg))
    assert_array_equal(mkl_svr_pre.predict(X_reg_pre), svr.predict(X_reg))
    assert_array_equal(mkl_svr_call.predict(X_reg), svr.predict(X_reg))


def test_mkl_internal_oneclasssvm_same_results_as_oneclasssvm():
    # Test that MKL internal OneClassSVM gives the same results as sklearn OneClassSVM
    oc_svm = svm.OneClassSVM(kernel="linear", nu=0.5)
    oc_svm.fit(X_oc)

    # Train MKL internal OneClassSVM
    mkl_oc_svm = mkl_svm.OneClassSVM(kernel="linear", nu=0.5)
    mkl_oc_svm.fit(X_oc)

    X_oc_pre = linear_kernel(X_oc)
    mkl_oc_svm_pre = mkl_svm.OneClassSVM(kernel="precomputed", nu=0.5)
    mkl_oc_svm_pre.fit(X_oc_pre)

    mkl_oc_svm_call = mkl_svm.OneClassSVM(kernel=linear_kernel, nu=0.5)
    mkl_oc_svm_call.fit(X_oc)

    assert_array_equal(mkl_oc_svm.predict(X_oc), oc_svm.predict(X_oc))
    assert_array_equal(mkl_oc_svm_pre.predict(X_oc_pre), oc_svm.predict(X_oc))
    assert_array_equal(mkl_oc_svm_call.predict(X_oc), oc_svm.predict(X_oc))


def test_mkl_internal_nusvc_same_results_as_nusvc():
    # Test that MKL internal NuSVC gives the same results as sklearn NuSVC
    nusvc = svm.NuSVC(kernel="linear", nu=0.5, random_state=42)
    nusvc.fit(X_class, y_class)

    # Train MKL internal NuSVC
    mkl_nusvc = mkl_svm.NuSVC(kernel="linear", nu=0.5, random_state=42)
    mkl_nusvc.fit(X_class, y_class)

    X_class_pre = linear_kernel(X_class)
    mkl_nusvc_pre = mkl_svm.NuSVC(kernel="precomputed", nu=0.5, random_state=42)
    mkl_nusvc_pre.fit(X_class_pre, y_class)

    mkl_nusvc_call = mkl_svm.NuSVC(kernel=linear_kernel, nu=0.5, random_state=42)
    mkl_nusvc_call.fit(X_class, y_class)

    assert_array_equal(mkl_nusvc.classes_, nusvc.classes_)
    assert_array_equal(mkl_nusvc_pre.classes_, nusvc.classes_)
    assert_array_equal(mkl_nusvc_call.classes_, nusvc.classes_)

    assert_array_equal(mkl_nusvc.predict(X_class), nusvc.predict(X_class))
    assert_array_equal(mkl_nusvc_pre.predict(X_class_pre), nusvc.predict(X_class))
    assert_array_equal(mkl_nusvc_call.predict(X_class), nusvc.predict(X_class))


def test_mkl_internal_nusvr_same_results_as_nusvr():
    # Test that MKL internal NuSVR gives the same results as sklearn NuSVR
    nusvr = svm.NuSVR(kernel="linear", nu=0.5, C=1.0)
    nusvr.fit(X_reg, y_reg)

    # Train MKL internal NuSVR
    mkl_nusvr = mkl_svm.NuSVR(kernel="linear", nu=0.5, C=1.0)
    mkl_nusvr.fit(X_reg, y_reg)

    X_reg_pre = linear_kernel(X_reg)
    mkl_nusvr_pre = mkl_svm.NuSVR(kernel="precomputed", nu=0.5, C=1.0)
    mkl_nusvr_pre.fit(X_reg_pre, y_reg)

    mkl_nusvr_call = mkl_svm.NuSVR(kernel=linear_kernel, nu=0.5, C=1.0)
    mkl_nusvr_call.fit(X_reg, y_reg)

    assert_array_equal(mkl_nusvr.predict(X_reg), nusvr.predict(X_reg))
    assert_array_equal(mkl_nusvr_pre.predict(X_reg_pre), nusvr.predict(X_reg))
    assert_array_equal(mkl_nusvr_call.predict(X_reg), nusvr.predict(X_reg))


def test_mkl_internal_svc_with_alpha_seeding_score():
    # Test that MKL internal SVC with alpha seeding gives at least same results as SVC
    svc = svm.SVC(kernel="linear", C=1.0, random_state=42, tol=1e-8)
    svc.fit(X_class, y_class)

    # Train MKL internal SVC with alpha seeding
    mkl_svc = mkl_svm.SVC(kernel="linear", C=1.0, random_state=42, tol=1e-8)
    mkl_svc.fit(X_class, y_class)

    # Alpha seeding using the optimal dual coefficients with noise
    rng = np.random.default_rng(seed=42)
    dual_coef_raw_noise = np.clip(
        mkl_svc.alpha_raw_
        + rng.uniform(-0.001, 0.001, mkl_svc.alpha_raw_.shape)
        * (mkl_svc.alpha_raw_ > 0),
        0,
        1,
    )
    mkl_svc.alpha_init_ = dual_coef_raw_noise  # Alpha seeding for the next fit
    mkl_svc.fit(X_class, y_class)

    assert mkl_svc.score(X_class, y_class) == pytest.approx(
        svc.score(X_class, y_class), rel=1e-3
    )


def test_mkl_internal_svr_with_alpha_seeding_score():
    # Test that MKL internal SVR with alpha seeding gives at least same results as SVR
    svr = svm.SVR(kernel="linear", C=1.0, epsilon=0.1, tol=1e-8)
    svr.fit(X_reg, y_reg)

    # Train MKL internal SVR with alpha seeding
    mkl_svr = mkl_svm.SVR(kernel="linear", C=1.0, epsilon=0.1, tol=1e-8)
    mkl_svr.fit(X_reg, y_reg)

    # Alpha seeding using the optimal dual coefficients with noise
    rng = np.random.default_rng(seed=42)
    dual_coef_raw_noise = np.clip(
        mkl_svr.alpha_raw_
        + rng.uniform(-0.001, 0.001, mkl_svr.alpha_raw_.shape)
        * (mkl_svr.alpha_raw_ > 0),
        0,
        1,
    )
    mkl_svr.alpha_init_ = dual_coef_raw_noise  # Alpha seeding for the next fit
    mkl_svr.fit(X_reg, y_reg)

    assert mkl_svr.score(X_reg, y_reg) == pytest.approx(
        svr.score(X_reg, y_reg), rel=1e-3
    )


def test_mkl_internal_nusvc_with_alpha_seeding_score():
    # Test that MKL internal NuSVC with alpha seeding gives
    # at least same results as sklearn NuSVC
    nusvc = svm.NuSVC(kernel="linear", nu=0.5, random_state=42, tol=1e-8)
    nusvc.fit(X_class, y_class)

    # Train MKL internal NuSVC with alpha seeding
    mkl_nusvc = mkl_svm.NuSVC(kernel="linear", nu=0.5, random_state=42, tol=1e-8)
    mkl_nusvc.fit(X_class, y_class)

    # Alpha seeding using the optimal dual coefficients with noise
    rng = np.random.default_rng(seed=42)
    dual_coef_raw_noise = np.clip(
        mkl_nusvc.alpha_raw_
        + rng.uniform(-0.001, 0.001, mkl_nusvc.alpha_raw_.shape)
        * (mkl_nusvc.alpha_raw_ > 0),
        0,
        1,
    )
    mkl_nusvc.alpha_init_ = dual_coef_raw_noise  # Alpha seeding for the next fit
    mkl_nusvc.fit(X_class, y_class)

    assert mkl_nusvc.score(X_class, y_class) == pytest.approx(
        nusvc.score(X_class, y_class), rel=1e-3
    )


def test_mkl_internal_nusvr_with_alpha_seeding_score():
    # Test that MKL internal NuSVR with alpha seeding gives
    # at least same results as sklearn NuSVR
    nusvr = svm.NuSVR(kernel="linear", nu=0.5, C=1.0, tol=1e-8)
    nusvr.fit(X_reg, y_reg)

    # Train MKL internal NuSVR with alpha seeding
    mkl_nusvr = mkl_svm.NuSVR(kernel="linear", nu=0.5, C=1.0, tol=1e-8)
    mkl_nusvr.fit(X_reg, y_reg)

    # Alpha seeding using the optimal dual coefficients with noise
    rng = np.random.default_rng(seed=42)
    dual_coef_raw_noise = np.clip(
        mkl_nusvr.alpha_raw_
        + rng.uniform(-0.001, 0.001, mkl_nusvr.alpha_raw_.shape)
        * (mkl_nusvr.alpha_raw_ > 0),
        0,
        1,
    )
    mkl_nusvr.alpha_init_ = dual_coef_raw_noise  # Alpha seeding for the next fit
    mkl_nusvr.fit(X_reg, y_reg)

    assert mkl_nusvr.score(X_reg, y_reg) == pytest.approx(
        nusvr.score(X_reg, y_reg), rel=1e-3
    )
