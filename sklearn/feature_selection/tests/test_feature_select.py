"""
Todo: cross-check the F-value with stats model
"""

from ..univariate_selection import (f_classif, f_regression, f_oneway,
                                    SelectPercentile, SelectKBest,
                                    SelectFpr, SelectFdr, SelectFwe,
                                    GenericUnivariateSelect)
import numpy as np
from numpy.testing import assert_array_equal
from scipy import stats
from sklearn.datasets.samples_generator import make_classification, \
                                                     make_regression


##############################################################################
# Test the score functions

def test_f_oneway_vs_scipy_stats():
    """Test that our f_oneway gives the same result as scipy.stats"""
    X1 = np.random.randn(10, 3)
    X2 = 1 + np.random.randn(10, 3)
    f, pv = stats.f_oneway(X1, X2)
    f2, pv2 = f_oneway(X1, X2)
    assert np.allclose(f, f2)
    assert np.allclose(pv, pv2)


def test_f_classif():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    F, pv = f_classif(X, Y)
    assert(F > 0).all()
    assert(pv > 0).all()
    assert(pv < 1).all()
    assert(pv[:5] < 0.05).all()
    assert(pv[5:] > 1.e-4).all()


def test_f_regression():
    """
    Test whether the F test yields meaningful results
    on a simple simulated regression problem
    """
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    F, pv = f_regression(X, Y)
    assert(F > 0).all()
    assert(pv > 0).all()
    assert(pv < 1).all()
    assert(pv[:5] < 0.05).all()
    assert(pv[5:] > 1.e-4).all()


def test_f_classif_multi_class():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    F, pv = f_classif(X, Y)
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
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectPercentile(f_classif, percentile=25)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='percentile',
                    param=25).fit(X, Y).transform(X)
    assert_array_equal(X_r, X_r2)
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
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectKBest(f_classif, k=5)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='k_best',
                    param=5).fit(X, Y).transform(X)
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
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectFpr(f_classif, alpha=0.0001)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='fpr',
                    param=0.0001).fit(X, Y).transform(X)
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
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectFdr(f_classif, alpha=0.0001)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='fdr',
                    param=0.0001).fit(X, Y).transform(X)
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
    X, Y = make_classification(n_samples=200, n_features=20,
                               n_informative=3, n_redundant=2,
                               n_repeated=0, n_classes=8,
                               n_clusters_per_class=1, flip_y=0.0,
                               class_sep=10, shuffle=False, random_state=0)

    univariate_filter = SelectFwe(f_classif, alpha=0.01)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_classif, mode='fwe',
                    param=0.01).fit(X, Y).transform(X)
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
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectPercentile(f_regression, percentile=25)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='percentile',
                    param=25).fit(X, Y).transform(X)
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
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectPercentile(f_regression, percentile=100)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='percentile',
                    param=100).fit(X, Y).transform(X)
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
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectKBest(f_regression, k=5)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='k_best',
                    param=5).fit(X, Y).transform(X)
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
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectFpr(f_regression, alpha=0.01)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='fpr',
                    param=0.01).fit(X, Y).transform(X)
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
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectFdr(f_regression, alpha=0.01)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='fdr',
                    param=0.01).fit(X, Y).transform(X)
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
    X, Y = make_regression(n_samples=200, n_features=20,
                           n_informative=5, shuffle=False, random_state=0)

    univariate_filter = SelectFwe(f_regression, alpha=0.01)
    X_r = univariate_filter.fit(X, Y).transform(X)
    X_r2 = GenericUnivariateSelect(f_regression, mode='fwe',
                    param=0.01).fit(X, Y).transform(X)
    assert_array_equal(X_r, X_r2)
    support = univariate_filter.get_support()
    gtruth = np.zeros(20)
    gtruth[:5] = 1
    assert(support[:5] == 1).all()
    assert(np.sum(support[5:] == 1) < 2)
