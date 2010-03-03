"""
Todo: cross-check the F-value with stats model
"""

from scikits.learn.feature_select import univ_selection  as fs
import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

def make_dataset(n_samples=50, n_features=20, k=5, seed=None, classif=True,
                 param=[1,1]):
    """
    Create a generic dataset for various tests
    """
    if classif:
        # classification
        x, y = fs.generate_dataset_classif(n_samples, n_features, k=k,
                                           seed=seed)
    else:
        # regression
        x, y = fs.generate_dataset_reg(n_samples, n_features, k=k, seed=seed)
        
    return x, y

def test_F_test_classif():
    """
    Test whether the F test yields meaningful results
    on a simple simulated classification problem
    """
    x, y = make_dataset()
    F, pv = fs.f_classif(x, y)
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
    x, y = make_dataset(classif=False)
    F, pv = fs.f_regression(x, y)
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
    x, y = make_dataset(param=[1,1,1])
    F, pv = fs.f_classif(x, y)
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
    x, y = make_dataset()
    univ_selection = fs.UnivSelection(score_func=fs.f_classif,
                                      select_args=(25,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_percentile_classif2():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the percentile heuristic
    """
    x, y = make_dataset()
    univ_selection = fs.UnivSelection(score_func=fs.f_classif,
                                      select_func=fs.select_percentile,
                                      select_args=(25,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_kbest_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the k best heuristic
    """
    x, y = make_dataset()
    univ_selection = fs.UnivSelection(score_func=fs.f_classif,
                                      select_func=fs.select_k_best,
                                      select_args=(5,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fpr_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    x, y = make_dataset()
    univ_selection = fs.UnivSelection(score_func=fs.f_classif,
                                      select_func=fs.select_fpr,
                                      select_args=(0.0001,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fdr_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    x, y = make_dataset(seed=3)
    univ_selection = fs.UnivSelection(score_func=fs.f_classif,
                                      select_func=fs.select_fdr,
                                      select_args=(0.01,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fwe_classif():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple classification problem
    with the fpr heuristic
    """
    x, y = make_dataset()
    univ_selection = fs.UnivSelection(score_func=fs.f_classif,
                                      select_func=fs.select_fwe,
                                      select_args=(0.01,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert(np.sum(np.abs(result-gtruth))<2)


def test_univ_fs_percentile_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the percentile heuristic
    """
    x, y = make_dataset(classif=False)
    univ_selection = fs.UnivSelection(score_func=fs.f_regression,
                                      select_args=(25,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_percentile_regression2():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the percentile heuristic
    """
    x, y = make_dataset(classif=False)
    univ_selection = fs.UnivSelection(score_func=fs.f_regression,
                                      select_func=fs.select_percentile,
                                      select_args=(25,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_kbest_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the k best heuristic
    """
    x, y = make_dataset(classif=False)
    univ_selection = fs.UnivSelection(score_func=fs.f_regression,
                                      select_func=fs.select_k_best,
                                      select_args=(5,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fpr_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    x, y = make_dataset(classif=False)
    univ_selection = fs.UnivSelection(score_func=fs.f_regression,
                                      select_func=fs.select_fpr,
                                      select_args=(0.001,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert(result[:5]==1).all()
    assert(np.sum(result[5:]==1)<2)

def test_univ_fs_fdr_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    x, y = make_dataset(seed=2, classif=False)
    univ_selection = fs.UnivSelection(score_func=fs.f_regression,
                                      select_func=fs.select_fdr,
                                      select_args=(0.01,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert_array_equal(result, gtruth)

def test_univ_fs_fwe_regression():
    """
    Test whether the relative univariate feature selection
    gets the correct items in a simple regression problem
    with the fpr heuristic
    """
    x, y = make_dataset(classif=False)
    univ_selection = fs.UnivSelection(score_func=fs.f_regression,
                                      select_func=fs.select_fwe,
                                      select_args=(0.01,))
    univ_selection.fit(x, y)
    result = univ_selection.support_.astype(int)
    gtruth = np.zeros(20)
    gtruth[:5]=1
    assert(result[:5]==1).all()
    assert(np.sum(result[5:]==1)<2)
