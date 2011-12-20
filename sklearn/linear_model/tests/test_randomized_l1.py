# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np
from numpy.testing import assert_equal, assert_raises

from sklearn.linear_model.randomized_l1 import lasso_stability_path, \
                                        RandomizedLasso, RandomizedLogistic
from sklearn.datasets import load_diabetes, load_iris
from sklearn.feature_selection import f_regression, f_classif
from sklearn.preprocessing import Scaler

diabetes = load_diabetes()
X = diabetes.data
y = diabetes.target
X = Scaler().fit_transform(X)
X = X[:, [2, 3, 6, 7, 8]]

# test that the feature score of the best features
F, _ = f_regression(X, y)


def test_lasso_stability_path():
    """Check lasso stability path"""
    # Load diabetes data and add noisy features
    a = 0.3
    coef_grid, scores_path = lasso_stability_path(X, y, a=a, random_state=42,
                                                  n_resampling=30)

    assert_equal(np.argsort(F)[-3:],
                 np.argsort(np.sum(scores_path, axis=1))[-3:])


def test_randomized_lasso():
    """Check randomized lasso"""
    a = 0.3
    selection_threshold = 0.5
    clf = RandomizedLasso(verbose=False, alpha=1, random_state=42, a=a,
                          selection_threshold=selection_threshold)
    feature_scores = clf.fit(X, y).scores_
    assert_equal(np.argsort(F)[-3:], np.argsort(feature_scores)[-3:])

    X_r = clf.transform(X)
    X_full = clf.inverse_transform(X_r)
    assert_equal(X_r.shape[1], np.sum(feature_scores > selection_threshold))
    assert_equal(X_full.shape, X.shape)

    clf = RandomizedLasso(verbose=False, alpha='aic', random_state=42, a=a)
    feature_scores = clf.fit(X, y).scores_
    assert_equal(feature_scores, X.shape[1] * [1.])

    clf = RandomizedLasso(verbose=False, a=-0.1)
    assert_raises(ValueError, clf.fit, X, y)

    clf = RandomizedLasso(verbose=False, a=1.1)
    assert_raises(ValueError, clf.fit, X, y)


def test_randomized_logistic():
    """Check randomized sparse logistic regression"""
    iris = load_iris()
    X = iris.data[:, :3]
    y = iris.target
    X = X[y != 2]
    y = y[y != 2]

    F, _ = f_classif(X, y)

    a = 0.3
    clf = RandomizedLogistic(verbose=False, C=1, random_state=42, a=a,
                             n_resampling=50)
    feature_scores = clf.fit(X, y).scores_
    assert_equal(np.argsort(F), np.argsort(feature_scores))
