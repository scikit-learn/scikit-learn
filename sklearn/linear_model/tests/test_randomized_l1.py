# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause
from tempfile import mkdtemp
import shutil

import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex

from sklearn.linear_model.randomized_l1 import (lasso_stability_path,
                                                RandomizedLasso,
                                                RandomizedLogisticRegression)
from sklearn.datasets import load_diabetes, load_iris
from sklearn.feature_selection import f_regression, f_classif
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model.base import _preprocess_data

diabetes = load_diabetes()
X = diabetes.data
y = diabetes.target
X = StandardScaler().fit_transform(X)
X = X[:, [2, 3, 6, 7, 8]]

# test that the feature score of the best features
F, _ = f_regression(X, y)


def test_lasso_stability_path():
    # Check lasso stability path
    # Load diabetes data and add noisy features
    scaling = 0.3
    coef_grid, scores_path = lasso_stability_path(X, y, scaling=scaling,
                                                  random_state=42,
                                                  n_resampling=30)

    assert_array_equal(np.argsort(F)[-3:],
                       np.argsort(np.sum(scores_path, axis=1))[-3:])


def test_randomized_lasso_error_memory():
    scaling = 0.3
    selection_threshold = 0.5
    tempdir = 5
    clf = RandomizedLasso(verbose=False, alpha=[1, 0.8], random_state=42,
                          scaling=scaling,
                          selection_threshold=selection_threshold,
                          memory=tempdir)
    assert_raises_regex(ValueError, "'memory' should either be a string or"
                        " a sklearn.externals.joblib.Memory instance",
                        clf.fit, X, y)


def test_randomized_lasso():
    # Check randomized lasso
    scaling = 0.3
    selection_threshold = 0.5

    # or with 1 alpha
    clf = RandomizedLasso(verbose=False, alpha=1, random_state=42,
                          scaling=scaling,
                          selection_threshold=selection_threshold)
    feature_scores = clf.fit(X, y).scores_
    assert_array_equal(np.argsort(F)[-3:], np.argsort(feature_scores)[-3:])

    # or with many alphas
    clf = RandomizedLasso(verbose=False, alpha=[1, 0.8], random_state=42,
                          scaling=scaling,
                          selection_threshold=selection_threshold)
    feature_scores = clf.fit(X, y).scores_
    assert_equal(clf.all_scores_.shape, (X.shape[1], 2))
    assert_array_equal(np.argsort(F)[-3:], np.argsort(feature_scores)[-3:])
    # test caching
    try:
        tempdir = mkdtemp()
        clf = RandomizedLasso(verbose=False, alpha=[1, 0.8], random_state=42,
                              scaling=scaling,
                              selection_threshold=selection_threshold,
                              memory=tempdir)
        feature_scores = clf.fit(X, y).scores_
        assert_equal(clf.all_scores_.shape, (X.shape[1], 2))
        assert_array_equal(np.argsort(F)[-3:], np.argsort(feature_scores)[-3:])
    finally:
        shutil.rmtree(tempdir)

    X_r = clf.transform(X)
    X_full = clf.inverse_transform(X_r)
    assert_equal(X_r.shape[1], np.sum(feature_scores > selection_threshold))
    assert_equal(X_full.shape, X.shape)

    clf = RandomizedLasso(verbose=False, alpha='aic', random_state=42,
                          scaling=scaling)
    feature_scores = clf.fit(X, y).scores_
    assert_array_equal(feature_scores, X.shape[1] * [1.])

    clf = RandomizedLasso(verbose=False, scaling=-0.1)
    assert_raises(ValueError, clf.fit, X, y)

    clf = RandomizedLasso(verbose=False, scaling=1.1)
    assert_raises(ValueError, clf.fit, X, y)


def test_dtype_match():
    # Check randomized lasso
    scaling = 0.3
    selection_threshold = 0.5

    X_64 = X.astype(np.float64)
    y_64 = y.astype(np.float64)
    X_32 = X.astype(np.float32)
    y_32 = y.astype(np.float32)

    # check type consistency 32 bits
    lasso_32 = RandomizedLasso(verbose=False, alpha=1, random_state=42,
                          scaling=scaling,
                          selection_threshold=selection_threshold)
    lasso_32.fit(X_32, y_32)
    assert_equal(lasso_32.coef_.dtype, X_32.dtype)

    # check type consistency 64 bits
    lasso_64 = RandomizedLasso(verbose=False, alpha=1, random_state=42,
                          scaling=scaling,
                          selection_threshold=selection_threshold)
    lasso_64.fit(X_64, y_64)
    assert_equal(lasso_64.coef_.dtype, X_64.dtype)

    assert_almost_equal(lasso_32.coef_, lasso_64.coef_.astype(np.float32))


def test_randomized_logistic():
    # Check randomized sparse logistic regression
    iris = load_iris()
    X = iris.data[:, [0, 2]]
    y = iris.target
    X = X[y != 2]
    y = y[y != 2]

    F, _ = f_classif(X, y)

    scaling = 0.3
    clf = RandomizedLogisticRegression(verbose=False, C=1., random_state=42,
                                       scaling=scaling, n_resampling=50,
                                       tol=1e-3)
    X_orig = X.copy()
    feature_scores = clf.fit(X, y).scores_
    assert_array_equal(X, X_orig)   # fit does not modify X
    assert_array_equal(np.argsort(F), np.argsort(feature_scores))

    clf = RandomizedLogisticRegression(verbose=False, C=[1., 0.5],
                                       random_state=42, scaling=scaling,
                                       n_resampling=50, tol=1e-3)
    feature_scores = clf.fit(X, y).scores_
    assert_array_equal(np.argsort(F), np.argsort(feature_scores))

    clf = RandomizedLogisticRegression(verbose=False, C=[[1., 0.5]])
    assert_raises(ValueError, clf.fit, X, y)


def test_randomized_logistic_sparse():
    # Check randomized sparse logistic regression on sparse data
    iris = load_iris()
    X = iris.data[:, [0, 2]]
    y = iris.target
    X = X[y != 2]
    y = y[y != 2]

    # center here because sparse matrices are usually not centered
    # labels should not be centered
    X, _, _, _, _ = _preprocess_data(X, y, True, True)

    X_sp = sparse.csr_matrix(X)

    F, _ = f_classif(X, y)

    scaling = 0.3
    clf = RandomizedLogisticRegression(verbose=False, C=1., random_state=42,
                                       scaling=scaling, n_resampling=50,
                                       tol=1e-3)
    feature_scores = clf.fit(X, y).scores_
    clf = RandomizedLogisticRegression(verbose=False, C=1., random_state=42,
                                       scaling=scaling, n_resampling=50,
                                       tol=1e-3)
    feature_scores_sp = clf.fit(X_sp, y).scores_
    assert_array_equal(feature_scores, feature_scores_sp)
