__author__ = "Henry Lin <hlin117@gmail.com>"

import scipy.sparse as sp
import numpy as np
from sklearn.preprocessing.discretization import Discretizer
from sklearn.preprocessing._discretization import binsearch
from sklearn.utils.testing import assert_array_almost_equal, assert_array_equal
from sklearn.utils.testing import assert_equal
from nose.tools import assert_raises

X = [[0, 0, 1],
     [1, 1, 2],
     [2, 0, 3]]

def test_binary_search():
    continuous = np.array([-10, -5])
    zero_interval = (-6, np.inf)
    search_points = np.array([-9, -8, -7])
    output = binsearch(continuous, zero_interval, search_points)
    # TODO: More test cases

def test_discretizer_mixed_signs():
    X_neg = [[-10, -2, 4],
             [-5,  8,  14]]

    dis = Discretizer(n_bins=5)
    discretized = dis.fit_transform(X_neg)

    expected_cut_points = [[-9, 0, 6],
                           [-8, 2, 8],
                           [-7, 4, 10],
                           [-6, 6, 12]]

    expected_zero_intervals = [(-6, np.inf), (0, 2), (-np.inf, 6)]
    assert_equal(expected_zero_intervals, dis.zero_intervals_)

    expected_searched_points = [[-9, 2, 8 ],
                                [-8, 4, 10],
                                [-7, 6, 12]]

    assert_array_equal(expected_searched_points, dis.searched_points_)
    assert_array_equal(expected_cut_points, dis.cut_points_)

def test_discretizer_sparse_with_categorical():
    # Also tests whether it can discretize with mixed signs.

    dis = Discretizer(n_bins=3, categorical_features=[2, 4])

    sparse_formats = [sp.bsr_matrix, sp.coo_matrix, sp.csc_matrix,
                      sp.csr_matrix, sp.dia_matrix, sp.dok_matrix, sp.lil_matrix]

    X = [[0, 2, 10, -20, 10, 0, 0],
         [0, 3, 11, -18, 11, 0, -8],
         [1, 5, 12, -16, 0, 3, 0],
         [9, 7, 13, -14, 0, 0, 0],
         [5, 9, 14, -12, 0, 0, 0],
         [7, 11, 15, -11, 6, -3, -4]]

    for format in sparse_formats:
        sparse_X = format(X)

        discretized = dis.fit_transform(sparse_X)

        expected_cut_points = [[3, 5, -17, -1, -5.333333],
                               [6, 8, -14, 1, -2.666666]]

        assert_array_almost_equal(expected_cut_points, dis.cut_points_)

        expected_output = [[0, 0, 1, 0, 0, 10, 10],
                           [0, 0, 1, 0, 1, 11, 11],
                           [0, 1, 2, 2, 0, 12, 0],
                           [2, 1, 0, 0, 0, 13, 0],
                           [1, 2, 0, 0, 0, 14, 0],
                           [2, 2, 0, 1, 2, 15, 6]]

        # The toarray() ensures that the output is sparse
        assert_array_equal(expected_output, discretized.toarray())

        # Checking to preserve sparsity
        if format is not sp.dia_matrix:
            assert_equal(sparse_X.size, discretized.size, "Failing with format {2},"
                                                          "sizes {0} and {1} differ" \
                         .format(sparse_X.size, discretized.size, format.__name__))


def test_discretizer_1d():
    pass # TODO

def test_discretizer_bad_n_bins():
    assert_raises(ValueError, Discretizer(n_bins=1).fit, X)

def test_discretizer_fit_transform():
    dis = Discretizer(n_bins=2)
    discretized = dis.fit_transform(X)

    expected_cut_points = [[1.0, 0.5, 2.0]]
    assert_array_equal(expected_cut_points, dis.cut_points_)
    expected = [[0, 0, 0],
                [1, 1, 1],
                [1, 0, 1]]
    assert_array_equal(expected, discretized)
    assert_equal(0, dis.searched_points_.size)  # Only when 2 bins

def test_discretizer_fit_transform_sparse():
    sparse_formats = [sp.bsr_matrix, sp.coo_matrix, sp.csc_matrix,
                      sp.csr_matrix, sp.dia_matrix, sp.dok_matrix, sp.lil_matrix]
    for format in sparse_formats:
        sparse_X = format(X)
        dis = Discretizer(n_bins=2)
        discretized = dis.fit_transform(sparse_X)

        expected_cut_points = [[1.0, 0.5, 2.0]]
        assert_array_equal(expected_cut_points, dis.cut_points_, "Failing with "
                           "format {0}".format(format.__name__))

        expected = [[0, 0, 0],
                    [1, 1, 1],
                    [1, 0, 1]]
        assert_array_equal(expected, discretized.toarray(), "Failing with format {0}"\
                           .format(format.__name__))

        # Checking to preserve sparsity
        if format is not sp.dia_matrix:
            assert_equal(sparse_X.size, discretized.size, "Failing with format {0}"\
                         .format(format.__name__))

def test_discretizer_fit_transform_single_cat_sparse():
    dis = Discretizer(n_bins=3, categorical_features=0)
    _categorical_helper(dis, sparse=True)

def test_discretizer_fit_transform_single_cat():
    dis = Discretizer(n_bins=3, categorical_features=0)
    _categorical_helper(dis)

def test_discretizer_fit_transform_cat():
    dis = Discretizer(n_bins=3, categorical_features=[0])
    _categorical_helper(dis)

def _categorical_helper(dis, sparse=False):
    if sparse:
        X_sparse = sp.csc_matrix(X)
        discretized = dis.fit_transform(X_sparse)
    else:
        discretized = dis.fit_transform(X)
    expected_cut_points = [[0.333333, 1.666666],
                           [0.666666, 2.333333]]

    assert_array_almost_equal(expected_cut_points, dis.cut_points_)

    expected = [[0, 0, 0],
                [2, 1, 1],
                [0, 2, 2]]

    if sparse:
        assert_array_equal(expected, discretized.toarray())
        assert_equal(X_sparse.size, discretized.size)
        return
    assert_array_equal(expected, discretized)

def test_discretizer_fit():
    dis = Discretizer(n_bins=3, categorical_features=[1]).fit(X)

    expected = [0, 2]
    assert_equal(expected, dis.continuous_features_)

    expected_cut_points = [[0.666666, 1.666666],
                           [1.333333, 2.333333]]

    assert_array_almost_equal(expected_cut_points, dis.cut_points_)

    expected = [0.0, 1.0]
    assert_array_almost_equal(expected, dis.min_)

    expected = [2.0, 3.0]
    assert_array_almost_equal(expected, dis.max_)

def test_discretizer_bad_cat_features():
    bad_cats = [[-1, 2], [0, 0], [0, 1, 2, 3]]
    for cats in bad_cats:
        try:
            dis = Discretizer(categorical_features=cats).fit(X)
        except:
            continue
        else:
            raise ValueError("Discretizer should have failed to fit with bad "
                             "input categories. Failed to fail on {0}" \
                             .format(cats))
