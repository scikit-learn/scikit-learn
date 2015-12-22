__author__ = "Henry Lin <hlin117@gmail.com>"

import scipy.sparse as sp
import numpy as np
from sklearn.preprocessing.discretization import Discretizer
from sklearn.preprocessing._discretization import binsearch
from sklearn.utils.testing import assert_array_almost_equal, assert_array_equal
from sklearn.utils.testing import assert_equal

X = [[0, 0, 1],
     [1, 1, 2],
     [2, 0, 3]]

def test_binary_search():
    continuous = np.array([-10, -5])
    zero_interval = (-6, np.inf)
    search_points = np.array([-9, -8, -7])
    output = binsearch(continuous, zero_interval, search_points)

def test_discretizer_negative():
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
    pass # TODO

def test_discretizer_1d():
    pass # TODO

def test_discretizer_bad_n_bins():
    try:
        dis = Discretizer(n_bins=1).fit(X)
    except:
        return
    else:
        raise ValueError("Discretizer needs to raise error when number of "
                         "bins is too low.")

def test_discretizer_fit_transform():
    dis = Discretizer(n_bins=2)
    discretized = dis.fit_transform(X)

    expected_cut_points = [[1.0, 0.5, 2.0]]
    assert_array_equal(expected_cut_points, dis.cut_points_)
    expected = [[0, 0, 0],
                [0, 1, 0],
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
                    [0, 1, 0],
                    [1, 0, 1]]
        assert_array_equal(expected, discretized, "Failing with format {0}"\
                           .format(format.__name__))

def test_discretizer_fit_transform_cat():
    dis = Discretizer(n_bins=3, categorical_features=[0])
    discretized = dis.fit_transform(X)

    expected = [[0, 0, 0],
                [2, 1, 1],
                [0, 2, 2]]

    assert_array_equal(expected, discretized)

def test_discretizer_fit():
    dis = Discretizer(n_bins=3, categorical_features=[1]).fit(X)

    expected = [0, 2]
    assert_equal(expected, dis.continuous_features_)

    expected = [[0.666666, 1.666666],
                [1.333333, 2.333333]]

    assert_array_almost_equal(expected, dis.cut_points_)

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

if __name__ == "__main__":
    test_binary_search()
