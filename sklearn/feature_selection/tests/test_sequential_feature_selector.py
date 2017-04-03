"""
Testing Sequential feature selection
"""
import numpy as np
from numpy.testing import assert_almost_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LinearRegression
from sklearn.datasets import load_iris
from sklearn.datasets import load_boston
from sklearn.feature_selection import SequentialFeatureSelector as SFS


def dict_compare_utility(d1, d2):
    assert d1.keys() == d2.keys(), "%s != %s" % (d1, d2)
    for i in d1:
        err_msg = ("d1[%s]['feature_idx']"
                   " != d2[%s]['feature_idx']" % (i, i))
        assert d1[i]['feature_idx'] == d1[i]["feature_idx"], err_msg
        assert_almost_equal(d1[i]['avg_score'],
                            d2[i]['avg_score'],
                            decimal=3,
                            err_msg=("d1[%s]['avg_score']"
                                     " != d2[%s]['avg_score']" % (i, i)))
        assert_almost_equal(d1[i]['cv_scores'],
                            d2[i]['cv_scores'],
                            decimal=3,
                            err_msg=("d1[%s]['cv_scores']"
                                     " != d2[%s]['cv_scores']" % (i, i)))


def test_run_default():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier()
    sfs = SFS(estimator=knn)
    sfs.fit(X, y)
    assert sfs.k_feature_idx_ == (3, )


def test_kfeatures_type_1():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier()
    expect = ('n_features_to_select must be a positive integer'
              ' between 1 and X.shape[1], got 0')
    sfs = SFS(estimator=knn,
              n_features_to_select=0)
    assert_raise_message(ValueError, expect, sfs.fit, X, y)


def test_kfeatures_type_2():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier()
    expect = 'n_features_to_select must be a positive integer or tuple'
    sfs = SFS(estimator=knn,
              n_features_to_select='abc')
    assert_raise_message(ValueError, expect, sfs.fit, X, y)


def test_kfeatures_type_3():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier()
    expect = ('n_features_to_select tuple min value must be in'
              ' range(1, X.shape[1]+1).')
    sfs = SFS(estimator=knn,
              n_features_to_select=(0, 5))
    assert_raise_message(ValueError, expect, sfs.fit, X, y)


def test_kfeatures_type_4():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier()
    expect = ('n_features_to_select tuple max value must be in'
              ' range(1, X.shape[1]+1).')
    sfs = SFS(estimator=knn,
              n_features_to_select=(1, 5))
    assert_raise_message(ValueError, expect, sfs.fit, X, y)


def test_kfeatures_type_5():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier()
    expect = ('he min n_features_to_select value must be'
              ' larger than the max n_features_to_select value.')
    sfs = SFS(estimator=knn,
              n_features_to_select=(3, 1))
    assert_raise_message(ValueError, expect, sfs.fit, X, y)


def test_knn_wo_cv():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier(n_neighbors=4)
    sfs1 = SFS(knn,
               n_features_to_select=3,
               forward=True,
               cv=0)
    sfs1 = sfs1.fit(X, y)
    expect = {1: {'avg_score': 0.95999999999999996,
                  'cv_scores': np.array([0.96]),
                  'feature_idx': (3,)},
              2: {'avg_score': 0.97333333333333338,
                  'cv_scores': np.array([0.97333333]),
                  'feature_idx': (2, 3)},
              3: {'avg_score': 0.97333333333333338,
                  'cv_scores': np.array([0.97333333]),
                  'feature_idx': (1, 2, 3)}}
    dict_compare_utility(d1=expect, d2=sfs1.subsets_)


def test_knn_cv3():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier(n_neighbors=4)
    sfs1 = SFS(knn,
               n_features_to_select=3,
               forward=True,
               cv=4)
    sfs1 = sfs1.fit(X, y)
    sfs1.subsets_
    expect = {1: {'avg_score': 0.95299145299145294,
                  'cv_scores': np.array([0.97435897,
                                         0.94871795,
                                         0.88888889,
                                         1.0]),
                  'feature_idx': (3,)},
              2: {'avg_score': 0.95993589743589736,
                  'cv_scores': np.array([0.97435897,
                                         0.94871795,
                                         0.91666667,
                                         1.0]),
                  'feature_idx': (2, 3)},
              3: {'avg_score': 0.97275641025641035,
                  'cv_scores': np.array([0.97435897,
                                         1.0,
                                         0.94444444,
                                         0.97222222]),
                  'feature_idx': (1, 2, 3)}}
    dict_compare_utility(d1=expect, d2=sfs1.subsets_)


def test_knn_option_sbs():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier(n_neighbors=4)
    sfs3 = SFS(knn,
               n_features_to_select=3,
               forward=False,
               cv=4)
    sfs3 = sfs3.fit(X, y)
    assert sfs3.k_feature_idx_ == (1, 2, 3)


def test_knn_option_sfs_tuplerange():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier(n_neighbors=3)
    sfs4 = SFS(knn,
               n_features_to_select=(1, 3),
               forward=True,
               cv=4)
    sfs4 = sfs4.fit(X, y)
    assert round(sfs4.k_score_, 3) == 0.967
    assert sfs4.k_feature_idx_ == (0, 2, 3)


def test_knn_scoring_metric():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier(n_neighbors=4)
    sfs5 = SFS(knn,
               n_features_to_select=3,
               forward=False,
               cv=4)
    sfs5 = sfs5.fit(X, y)
    assert round(sfs5.k_score_, 4) == 0.9728

    sfs6 = SFS(knn,
               n_features_to_select=3,
               forward=False,
               cv=4)
    sfs6 = sfs6.fit(X, y)
    assert round(sfs6.k_score_, 4) == 0.9728

    sfs7 = SFS(knn,
               n_features_to_select=3,
               forward=False,
               scoring='f1_macro',
               cv=4)
    sfs7 = sfs7.fit(X, y)
    assert round(sfs7.k_score_, 4) == 0.9727


def test_regression():
    boston = load_boston()
    X, y = boston.data, boston.target
    lr = LinearRegression()
    sfs_r = SFS(lr,
                n_features_to_select=13,
                forward=True,
                cv=10)
    sfs_r = sfs_r.fit(X, y)
    assert len(sfs_r.k_feature_idx_) == 13
    assert round(sfs_r.k_score_, 4) == 0.2001


def test_regression_in_tuplerange():
    boston = load_boston()
    X, y = boston.data, boston.target
    lr = LinearRegression()
    sfs_r = SFS(lr,
                n_features_to_select=(1, 13),
                forward=True,
                cv=10)
    sfs_r = sfs_r.fit(X, y)
    assert len(sfs_r.k_feature_idx_) == 9
    assert round(sfs_r.k_score_, 4) == 0.2991, sfs_r.k_score_


def test_transform_not_fitted():
    iris = load_iris()
    X = iris.data
    y = iris.target
    knn = KNeighborsClassifier(n_neighbors=4)

    sfs1 = SFS(knn,
               n_features_to_select=2,
               forward=True,
               cv=0)

    expect = ("This SequentialFeatureSelector instance is not fitted yet.")

    assert_raise_message(ValueError, expect, sfs1.transform, X)
