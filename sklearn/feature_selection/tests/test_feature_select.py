"""
Todo: cross-check the F-value with stats model
"""

from sklearn.feature_selection import (chi2, f_classif, f_oneway, f_regression,
                                       SelectPercentile, SelectKBest,
                                       SelectFpr, SelectFdr, SelectFwe,
                                       GenericUnivariateSelect)
from nose.tools import assert_equal, assert_true
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from scipy import stats, sparse
from sklearn.datasets.samples_generator import make_classification, \
                                                     make_regression


##############################################################################
# Test the score functions

def test_f_oneway_vs_scipy_stats():
    """Test that our f_oneway gives the same result as scipy.stats"""
    rng = np.random.RandomState(0)
    X1 = rng.randn(10, 3)
    X2 = 1 + rng.randn(10, 3)
    f, pv = stats.f_oneway(X1, X2)
    f2, pv2 = f_oneway(X1, X2)
    assert_true(np.allclose(f, f2))
    assert_true(np.allclose(pv, pv2))


def test_f_oneway_ints():
    # Smoke test f_oneway on integers: that it does raise casting errors
    # with recent numpys
    f_oneway(np.random.randint(10, size=(10, 10)), np.arange(10))


def test_f_classif():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    F, pv = f_classif(X, y)
    F_sparse,  pv_sparse = f_classif(sparse.csr_matrix(X), y)
    assert(F > 0).all()
    assert(pv > 0).all()
    assert(pv < 1).all()
    assert(pv[:5] < 0.05).all()
    assert(pv[5:] > 1.e-4).all()
    assert_array_almost_equal(F_sparse, F)
    assert_array_almost_equal(pv_sparse, pv)


def test_f_regression():
    """
    Test whether the F test yields meaningful results
    on a simple simulated regression problem
    """
    X, y = make_regression(n_samples=200, n_features=20,
        n_informative=5, shuffle=False, random_state=0)

    F, pv = f_regression(X, y)
    assert(F > 0).all()
    assert(pv > 0).all()
    assert(pv < 1).all()
    assert(pv[:5] < 0.05).all()
    assert(pv[5:] > 1.e-4).all()

    # again without centering, compare with sparse
    F, pv = f_regression(X, y, center=False)
    F_sparse, pv_sparse = f_regression(sparse.csr_matrix(X), y, center=False)
    assert_array_almost_equal(F_sparse, F)
    assert_array_almost_equal(pv_sparse, pv)


def test_f_regression_input_dtype():
    """
    Test whether f_regression returns the same value
    for any numeric data_type
    """
    rng = np.random.RandomState(0)
    X = rng.rand(10, 20)
    y = np.arange(10).astype(np.int)

    F1, pv1 = f_regression(X, y)
    F2, pv2 = f_regression(X, y.astype(np.float))
    assert_array_almost_equal(F1, F2, 5)
    assert_array_almost_equal(pv1, pv2, 5)


def test_f_classif_multi_class():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    F, pv = f_classif(X, y)
    assert(F > 0).all()
    assert(pv > 0).all()
    assert(pv < 1).all()
    assert(pv[:5] < 0.05).all()
    assert(pv[5:] > 1.e-5).all()


def test_select_percentile_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the percentile heuristic
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectPercentile(f_classif, percentile=25)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='percentile',
                    param=25).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


def test_select_percentile_classif_sparse():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the percentile heuristic
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)
    X = sparse.csr_matrix(X)
    univariate_filter = SelectPercentile(f_classif, percentile=25)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='percentile',
                    param=25).fit(X, y).transform(X)
    assert_array_equal(X_r.toarray(), X_r2.toarray())
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


##############################################################################
# Test univariate selection in classification settings

def test_select_kbest_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the k best heuristic
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectKBest(f_classif, k=5)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='k_best',
                    param=5).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


def test_select_fpr_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectFpr(f_classif, alpha=0.0001)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='fpr',
                    param=0.0001).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


def test_select_fdr_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectFdr(f_classif, alpha=0.0001)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='fdr',
                    param=0.0001).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


def test_select_fwe_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    X, y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectFwe(f_classif, alpha=0.01)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='fwe',
                    param=0.01).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert(np.sum(np.abs(support - gtruth)) < 2)


##############################################################################
# Test univariate selection in regression settings

def test_select_percentile_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the percentile heuristic
    """
    X, y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectPercentile(f_regression, percentile=25)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='percentile',
                    param=25).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)
    X_2 = X.copy()
    X_2[:, np.logical_not(support)] = 0
    assert_array_equal(X_2, univariate_filter.inverse_transform(X_r))


def test_select_percentile_regression_full():
    """
    Test whether the relative univariate feature selection
    selects all features when '100%' is asked.
    """
    X, y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectPercentile(f_regression, percentile=100)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='percentile',
                    param=100).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.ones(20)
    assert_array_equal(support, gtruth)


def test_select_kbest_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the k best heuristic
    """
    X, y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectKBest(f_regression, k=5)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='k_best',
                    param=5).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


def test_select_fpr_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    X, y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectFpr(f_regression, alpha=0.01)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='fpr',
                    param=0.01).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert(support[:5] == 1).all()
    assert(np.sum(support[5:] == 1) < 3)


def test_select_fdr_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fdr heuristic
    """
    X, y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectFdr(f_regression, alpha=0.01)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='fdr',
                    param=0.01).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert_array_equal(support, gtruth)


def test_select_fwe_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fwe heuristic
    """
    X, y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectFwe(f_regression, alpha=0.01)
    X_r = univariate_filter.fit(X, y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='fwe',
                    param=0.01).fit(X, y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert(support[:5] == 1).all()
    assert(np.sum(support[5:] == 1) < 2)


def test_selectkbest_tiebreaking():
    """Test whether SelectKBest actually selects k features in case of ties.

    Prior to 0.11, SelectKBest would return more features than requested.
    """
    X = [[1, 0, 0], [0, 1, 1]]
    y = [0, 1]

    X1 = SelectKBest(chi2, k=1).fit_transform(X, y)
    assert_equal(X1.shape[1], 1)

    X2 = SelectKBest(chi2, k=2).fit_transform(X, y)
    assert_equal(X2.shape[1], 2)
