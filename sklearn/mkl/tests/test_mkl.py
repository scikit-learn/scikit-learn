"""
Unit tests for the Multiple Kernel Learning (MKL) module.
These tests cover the fundamentals of MKL.
For algorithm-specific tests, see the corresponding test files in "_aglo" directory.
"""

import pytest

from sklearn import mkl, svm
from sklearn.datasets import make_classification, make_regression
from sklearn.metrics.pairwise import linear_kernel, rbf_kernel
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

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

ten_kernels_single_params = {
    "kernels": ["rbf"],
    "kernels_scopes": ["single"],
}


def assert_fit_predict_works_for_all_algorithms(clf, X_ok, X_not_ok, y, n_kernels=2):
    for algo in mkl.MKL_ALGORITHMS:
        clf.set_params(algo=algo)

        with pytest.raises(ValueError):
            clf.fit(X_not_ok, y)
        clf.fit(X_ok, y)

        assert hasattr(clf, "weights_")
        assert hasattr(clf, "n_kernels_")
        assert clf.n_kernels_ == n_kernels

        with pytest.raises(ValueError):
            clf.predict(X_not_ok)
        clf.predict(X_ok)


def test_mklc_precomputed_kernels():
    # Test MKLC with precomputed kernels.
    clf = mkl.MKLC(C=1.0, random_state=42)

    assert_fit_predict_works_for_all_algorithms(
        clf, [K_linear_class, K_rbf_class], X_class, y_class
    )


def test_mklr_precomputed_kernels():
    # Test MKLR with precomputed kernels.
    reg = mkl.MKLR(C=1.0, epsilon=0.1, random_state=42)

    assert_fit_predict_works_for_all_algorithms(
        reg, [K_linear_reg, K_rbf_reg], X_reg, y_reg
    )


def test_oneclassmkl_precomputed_kernels():
    # Test OneClassMKL with precomputed kernels.
    oc = mkl.OneClassMKL(nu=0.5, random_state=42)

    assert_fit_predict_works_for_all_algorithms(oc, [K_linear_oc, K_rbf_oc], X_oc, None)


def test_mklc_string_list_kernels():
    # Test MKLC with string list kernels.
    clf = mkl.MKLC(kernels=["linear", "rbf"], C=1.0, random_state=42)

    assert_fit_predict_works_for_all_algorithms(
        clf, X_class, [K_linear_class, K_rbf_class], y_class
    )


def test_mklr_string_list_kernels():
    # Test MKLR with string list kernels.
    reg = mkl.MKLR(kernels=["linear", "rbf"], C=1.0, epsilon=0.1, random_state=42)

    assert_fit_predict_works_for_all_algorithms(
        reg, X_reg, [K_linear_reg, K_rbf_reg], y_reg
    )


def test_oneclassmkl_string_list_kernels():
    # Test OneClassMKL with string list kernels.
    oc = mkl.OneClassMKL(kernels=["linear", "rbf"], nu=0.5, random_state=42)

    assert_fit_predict_works_for_all_algorithms(oc, X_oc, [K_linear_oc, K_rbf_oc], None)


def test_mklc_callable_list_kernels():
    # Test MKLC with callable list kernels.
    clf = mkl.MKLC(kernels=[linear_kernel, rbf_kernel], C=1.0, random_state=42)

    assert_fit_predict_works_for_all_algorithms(
        clf, X_class, [K_linear_class, K_rbf_class], y_class
    )


def test_mklr_callable_list_kernels():
    # Test MKLR with callable list kernels.
    reg = mkl.MKLR(
        kernels=[linear_kernel, rbf_kernel], C=1.0, epsilon=0.1, random_state=42
    )

    assert_fit_predict_works_for_all_algorithms(
        reg, X_reg, [K_linear_reg, K_rbf_reg], y_reg
    )


def test_oneclassmkl_callable_list_kernels():
    # Test OneClassMKL with callable list kernels.
    oc = mkl.OneClassMKL(kernels=[linear_kernel, rbf_kernel], nu=0.5, random_state=42)

    assert_fit_predict_works_for_all_algorithms(oc, X_oc, [K_linear_oc, K_rbf_oc], None)


def test_mkl_invalid_string_in_kernel_list():
    # Test MKLR with invalid string kernel.
    clf = mkl.MKLC(kernels=["invalid_kernel"], C=1.0, random_state=42)
    with pytest.raises(ValueError, match="Invalid kernel"):
        clf.fit(X_class, y_class)

    reg = mkl.MKLR(kernels=["invalid_kernel"], C=1.0, epsilon=0.1, random_state=42)
    with pytest.raises(ValueError, match="Invalid kernel"):
        reg.fit(X_reg, y_reg)

    oc = mkl.OneClassMKL(kernels=["invalid_kernel"], nu=0.5, random_state=42)
    with pytest.raises(ValueError, match="Invalid kernel"):
        oc.fit(X_oc)


def test_mkl_invalid_type_in_kernel_list():
    # Test MKLC with invalid type in kernel list.
    clf = mkl.MKLC(kernels=[112_145_163_165_163_040_151_163_040_114_157_166_145])
    with pytest.raises(ValueError, match="Invalid kernel"):
        clf.fit(X_class, y_class)

    reg = mkl.MKLR(kernels=[112_145_163_165_163_040_151_163_40_114_151_146_145])
    with pytest.raises(ValueError, match="Invalid kernel"):
        reg.fit(X_reg, y_reg)

    oc = mkl.OneClassMKL(kernels=[112_145_163_165_163_040_151_163_040_107_157_144])
    with pytest.raises(ValueError, match="Invalid kernel"):
        oc.fit(X_oc)


def test_mkl_kernel_scope_warning_when_precomputed():
    # Test MKLC with kernel scope when precomputed.
    with pytest.warns(Warning):
        mkl.MKLC(kernels_scopes=["all"], C=1.0, random_state=42)

    with pytest.warns(Warning):
        mkl.MKLR(kernels_scopes=["all"], C=1.0, epsilon=0.1, random_state=42)

    with pytest.warns(Warning):
        mkl.OneClassMKL(kernels_scopes=["all"], nu=0.5, random_state=42)


def test_mkl_kernel_param_grids_warning_when_precomputed():
    # Test MKLC with kernel param grids when precomputed.
    with pytest.warns(Warning):
        mkl.MKLC(kernels_param_grids=[{"gamma": [0.1, 1, 10]}], C=1.0)

    with pytest.warns(Warning):
        mkl.MKLR(kernels_param_grids=[{"gamma": [0.1, 1, 10]}], C=1.0, epsilon=0.1)

    with pytest.warns(Warning):
        mkl.OneClassMKL(kernels_param_grids=[{"gamma": [0.1, 1, 10]}], nu=0.5)


def test_mkl_kernel_scope_wrong_length():
    # Test MKLC with wrong length of kernel scope.
    clf = mkl.MKLC(kernels=["linear", "rbf"], kernels_scopes=["all"], C=1.0)
    with pytest.raises(ValueError, match="Invalid 'kernels_scopes'"):
        clf.fit(X_class, y_class)

    reg = mkl.MKLR(
        kernels=["linear", "rbf"], kernels_scopes=["all"], C=1.0, epsilon=0.1
    )
    with pytest.raises(ValueError, match="Invalid 'kernels_scopes'"):
        reg.fit(X_reg, y_reg)

    oc = mkl.OneClassMKL(kernels=["linear", "rbf"], kernels_scopes=["all"], nu=0.5)
    with pytest.raises(ValueError, match="Invalid 'kernels_scopes'"):
        oc.fit(X_oc)


def test_mkl_kernels_param_grids_wrong_length():
    # Test MKLC with wrong length of kernel param grids.
    clf = mkl.MKLC(
        kernels=["linear", "rbf"], kernels_param_grids=[{"gamma": [0.1, 1, 10]}]
    )
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        clf.fit(X_class, y_class)

    reg = mkl.MKLR(
        kernels=["linear", "rbf"], kernels_param_grids=[{"gamma": [0.1, 1, 10]}]
    )
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        reg.fit(X_reg, y_reg)

    oc = mkl.OneClassMKL(
        kernels=["linear", "rbf"], kernels_param_grids=[{"gamma": [0.1, 1, 10]}]
    )
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        oc.fit(X_oc)


def test_mkl_kernel_param_grids_wrong_key_type():
    # Test MKLC with wrong keys in kernel param grids.
    clf = mkl.MKLC(kernels=["linear"], kernels_param_grids=[{1: [0.1, 1, 10]}])
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        clf.fit(X_class, y_class)

    reg = mkl.MKLR(kernels=["linear"], kernels_param_grids=[{2: [0.1, 1, 10]}])
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        reg.fit(X_reg, y_reg)

    oc = mkl.OneClassMKL(kernels=["linear"], kernels_param_grids=[{3: [0.1, 1, 10]}])
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        oc.fit(X_oc)


def test_mkl_kernel_param_grids_wrong_value_type():
    # Test MKLC with wrong values in kernel param grids.
    clf = mkl.MKLC(kernels=["rbf"], kernels_param_grids=[{"gamma": "invalid"}])
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        clf.fit(X_class, y_class)

    reg = mkl.MKLR(kernels=["rbf"], kernels_param_grids=[{"gamma": "invalid"}])
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        reg.fit(X_reg, y_reg)

    oc = mkl.OneClassMKL(kernels=["rbf"], kernels_param_grids=[{"gamma": "invalid"}])
    with pytest.raises(ValueError, match="Invalid 'kernels_param_grids'"):
        oc.fit(X_oc)


def test_mkl_svm_parameters_through_svm_params():
    # Test MKL with SVM parameters through svm_params.
    with pytest.warns(Warning):
        mkl.MKLC(svm_params={"kernel": "linear", "C": 1.0})

    with pytest.warns(Warning):
        mkl.MKLR(svm_params={"kernel": "linear", "C": 1.0, "epsilon": 0.1})

    with pytest.warns(Warning):
        mkl.OneClassMKL(svm_params={"kernel": "linear", "nu": 0.5})


def assert_valid_number_of_kernels_for_all_algorithms(clf, n_kernels, X, y):
    for algo in mkl.MKL_ALGORITHMS:
        clf.set_params(algo=algo)
        clf.fit(X, y)

        assert hasattr(clf, "n_kernels_")
        assert clf.n_kernels_ == n_kernels


def test_mklc_valid_number_of_kernels():
    # Test MKLC with valid number of kernels.
    # Six kernels ("all" scopes)
    clf = mkl.MKLC(**six_kernels_params, C=1.0, random_state=42)
    assert_valid_number_of_kernels_for_all_algorithms(clf, 6, X_class, y_class)

    # Ten kernels ("single" scope)
    clf = mkl.MKLC(**ten_kernels_single_params, C=1.0, random_state=42)
    assert_valid_number_of_kernels_for_all_algorithms(clf, 10, X_class, y_class)


def test_mklr_valid_number_of_kernels():
    # Test MKLR with valid number of kernels.
    # Six kernels ("all" scopes)
    reg = mkl.MKLR(**six_kernels_params, C=1.0, epsilon=0.1, random_state=42)
    assert_valid_number_of_kernels_for_all_algorithms(reg, 6, X_reg, y_reg)

    # Ten kernels ("single" scope)
    reg = mkl.MKLR(**ten_kernels_single_params, C=1.0, epsilon=0.1, random_state=42)
    assert_valid_number_of_kernels_for_all_algorithms(reg, 10, X_reg, y_reg)


def test_oneclassmkl_valid_number_of_kernels():
    # Test OneClassMKL with valid number of kernels.
    # Six kernels ("all" scopes)
    oc = mkl.OneClassMKL(**six_kernels_params, nu=0.5, random_state=42)
    assert_valid_number_of_kernels_for_all_algorithms(oc, 6, X_oc, None)

    # Ten kernels ("single" scope)
    oc = mkl.OneClassMKL(**ten_kernels_single_params, nu=0.5, random_state=42)
    assert_valid_number_of_kernels_for_all_algorithms(oc, 10, X_oc, None)


def test_mklc_wrong_number_of_kernels_when_predicting():
    # Test MKLC with wrong number of kernels when predicting.
    clf = mkl.MKLC(C=1.0, random_state=42)

    for algo in mkl.MKL_ALGORITHMS:
        clf.set_params(algo=algo)
        clf.fit([K_linear_class, K_rbf_class], y_class)

        with pytest.raises(ValueError):
            clf.predict([K_linear_class, K_rbf_class, K_rbf_class])


def test_mklr_wrong_number_of_kernels_when_predicting():
    # Test MKLR with wrong number of kernels when predicting.
    reg = mkl.MKLR(C=1.0, epsilon=0.1, random_state=42)

    for algo in mkl.MKL_ALGORITHMS:
        reg.set_params(algo=algo)
        reg.fit([K_linear_reg, K_rbf_reg], y_reg)

        with pytest.raises(ValueError):
            reg.predict([K_linear_reg, K_rbf_reg, K_rbf_reg])


def test_oneclassmkl_wrong_number_of_kernels_when_predicting():
    # Test OneClassMKL with wrong number of kernels when predicting.
    oc = mkl.OneClassMKL(nu=0.5, random_state=42)

    for algo in mkl.MKL_ALGORITHMS:
        oc.set_params(algo=algo)
        oc.fit([K_linear_oc, K_rbf_oc])

        with pytest.raises(ValueError):
            oc.predict([K_linear_oc, K_rbf_oc, K_rbf_oc])


def test_mklc_pipeline_as_a_transformer():
    # Test MKLC as a transformer in a pipeline.
    clf_mkl = mkl.MKLC(C=1.0, random_state=42)
    clf_svm = svm.SVC(C=1.0, kernel="precomputed")
    pipe = Pipeline([("mkl", clf_mkl), ("svm", clf_svm)])

    for algo in mkl.MKL_ALGORITHMS:
        pipe.set_params(mkl__algo=algo)
        pipe.fit([K_linear_class, K_rbf_class], y_class)
        pipe.predict([K_linear_class, K_rbf_class])


def test_mklr_pipeline_as_a_transformer():
    # Test MKLR as a transformer in a pipeline.
    reg_mkl = mkl.MKLR(C=1.0, epsilon=0.1, random_state=42)
    reg_svm = svm.SVR(C=1.0, kernel="precomputed")
    pipe = Pipeline([("mkl", reg_mkl), ("svm", reg_svm)])

    for algo in mkl.MKL_ALGORITHMS:
        pipe.set_params(mkl__algo=algo)
        pipe.fit([K_linear_reg, K_rbf_reg], y_reg)
        pipe.predict([K_linear_reg, K_rbf_reg])


def test_oneclassmkl_pipeline_as_a_transformer():
    # Test OneClassMKL as a transformer in a pipeline.
    oc_mkl = mkl.OneClassMKL(nu=0.5, random_state=42)
    oc_svm = svm.OneClassSVM(nu=0.5, kernel="precomputed")
    pipe = Pipeline([("mkl", oc_mkl), ("svm", oc_svm)])

    for algo in mkl.MKL_ALGORITHMS:
        pipe.set_params(mkl__algo=algo)
        pipe.fit([K_linear_oc, K_rbf_oc])
        pipe.predict([K_linear_oc, K_rbf_oc])


def test_mklc_pipeline_as_an_estimator():
    # Test MKLC as an estimator in a pipeline.
    sc = StandardScaler()
    clf = mkl.MKLC(kernels=["linear", "rbf"], C=1.0, random_state=42)
    pipe = Pipeline([("sc", sc), ("mkl", clf)])

    for algo in mkl.MKL_ALGORITHMS:
        pipe.set_params(mkl__algo=algo)
        pipe.fit(X_class, y_class)
        pipe.predict(X_class)


def test_mklr_pipeline_as_an_estimator():
    # Test MKLR as an estimator in a pipeline.
    sc = StandardScaler()
    reg = mkl.MKLR(kernels=["linear", "rbf"], C=1.0, epsilon=0.1, random_state=42)
    pipe = Pipeline([("sc", sc), ("mkl", reg)])

    for algo in mkl.MKL_ALGORITHMS:
        pipe.set_params(mkl__algo=algo)
        pipe.fit(X_reg, y_reg)
        pipe.predict(X_reg)


def test_oneclassmkl_pipeline_as_an_estimator():
    # Test OneClassMKL as an estimator in a pipeline.
    sc = StandardScaler()
    oc = mkl.OneClassMKL(kernels=["linear", "rbf"], nu=0.5, random_state=42)
    pipe = Pipeline([("sc", sc), ("mkl", oc)])

    for algo in mkl.MKL_ALGORITHMS:
        pipe.set_params(mkl__algo=algo)
        pipe.fit(X_oc)
        pipe.predict(X_oc)


def test_mklc_as_a_callable_kernel():
    # Test MKLC as a callable kernel.
    clf_mkl = mkl.MKLC(kernels=["linear", "rbf"], C=1.0, random_state=42)

    for algo in mkl.MKL_ALGORITHMS:
        clf_mkl.set_params(algo=algo)
        clf_mkl.fit(X_class, y_class)

        # Check if the MKL instance can be used as a kernel
        clf_svm = svm.SVC(kernel=clf_mkl, C=1.0)
        clf_svm.fit(X_class, y_class)
        clf_svm.predict(X_class)


def test_mklr_as_a_callable_kernel():
    # Test MKLR as a callable kernel.
    reg_mkl = mkl.MKLR(kernels=["linear", "rbf"], C=1.0, epsilon=0.1, random_state=42)

    for algo in mkl.MKL_ALGORITHMS:
        reg_mkl.set_params(algo=algo)
        reg_mkl.fit(X_reg, y_reg)

        # Check if the MKL instance can be used as a kernel
        reg_svm = svm.SVR(kernel=reg_mkl, C=1.0)
        reg_svm.fit(X_reg, y_reg)
        reg_svm.predict(X_reg)


def test_oneclassmkl_as_a_callable_kernel():
    # Test OneClassMKL as a callable kernel.
    oc_mkl = mkl.OneClassMKL(kernels=["linear", "rbf"], nu=0.5, random_state=42)

    for algo in mkl.MKL_ALGORITHMS:
        oc_mkl.set_params(algo=algo)
        oc_mkl.fit(X_oc)

        # Check if the MKL instance can be used as a kernel
        oc_svm = svm.OneClassSVM(kernel=oc_mkl, nu=0.5)
        oc_svm.fit(X_oc)
        oc_svm.predict(X_oc)
