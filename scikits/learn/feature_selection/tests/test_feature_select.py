"""
Todo: cross-check the F-value with stats model
"""

from scikits.learn.feature_selection import univariate_selection  as us
import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises
import scikits.learn.datasets.samples_generator as sg


def test_F_test_classif():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None)
    F, pv = us.f_classif(X, Y)
    assert(F>0).all()
    assert(pv>0).all()
    assert(pv<1).all()
    assert(pv[:5]<0.05).all()
    assert(pv[5:]>1.e-4).all()

def test_F_test_reg():
    """
    Test whether the F test yields meaningful results
    on a simple simulated regression problem
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None)
    F, pv = us.f_regression(X, Y)
    assert(F>0).all()
    assert(pv>0).all()
    assert(pv<1).all()
    assert(pv[:5]<0.05).all()
    assert(pv[5:]>1.e-4).all()

def test_F_test_multi_class():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None,param=[1,1,1])
    F, pv = us.f_classif(X, Y)
    assert(F>0).all()
    assert(pv>0).all()
    assert(pv<1).all()
    assert(pv[:5]<0.05).all()
    assert(pv[5:]>1.e-4).all()

def test_univ_fs_percentile_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the percentile heuristic
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectPercentile(us.f_classif)
    X_r = univariate_filter.fit(X, Y).transform(X, percentile=25)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_kbest_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the k best heuristic
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectKBest(us.f_classif)
    X_r = univariate_filter.fit(X, Y).transform(X, k=5)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fpr_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectFpr(us.f_classif)
    X_r = univariate_filter.fit(X, Y).transform(X, alpha=0.0001)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fdr_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=3)
    univariate_filter = us.SelectFdr(us.f_classif)
    X_r = univariate_filter.fit(X, Y).transform(X, alpha=0.01)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fwe_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    X, Y = sg.test_dataset_classif(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectFwe(us.f_classif)
    X_r = univariate_filter.fit(X, Y).transform(X, alpha=0.01)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert(np.sum(np.abs(result-gtruth))<2)







def test_univ_fs_percentile_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the percentile heuristic
    """
    X, Y = sg.test_dataset_reg(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectPercentile(us.f_regression)
    X_r = univariate_filter.fit(X, Y).transform(X, percentile=25)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_full_percentile_regression():
    """
    Test whether the relative univariate feature selection
    selects all features when '100%' is asked.
    """
    X, Y = sg.test_dataset_reg(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectPercentile(us.f_regression)
    X_r = univariate_filter.fit(X, Y).transform(X, percentile=100)
    result = univariate_filter.support_.astype(int)
    gtruth = np.ones(20)
    assert_array_equal(result, gtruth)

def test_univ_fs_kbest_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the k best heuristic
    """
    X, Y = sg.test_dataset_reg(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectKBest(us.f_regression)
    X_r = univariate_filter.fit(X, Y).transform(X, k=5)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fpr_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    X, Y = sg.test_dataset_reg(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectFpr(us.f_regression)
    X_r = univariate_filter.fit(X, Y).transform(X, alpha=0.01)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert(result[:5]==1).all()
    assert(np.sum(result[5:]==1)<3)

def test_univ_fs_fdr_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    X, Y = sg.test_dataset_reg(n_samples=50, n_features=20, k=5,
                                           seed=2)
    univariate_filter = us.SelectFdr(us.f_regression)
    X_r = univariate_filter.fit(X, Y).transform(X, alpha=0.01)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fwe_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    X, Y = sg.test_dataset_reg(n_samples=50, n_features=20, k=5,
                                           seed=None)
    univariate_filter = us.SelectFwe(us.f_regression)
    X_r = univariate_filter.fit(X, Y).transform(X, alpha=0.01)
    result = univariate_filter.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert(result[:5]==1).all()
    assert(np.sum(result[5:]==1)<2)
