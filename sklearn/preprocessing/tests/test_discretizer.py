__author__ = "Henry Lin <hlin117@gmail.com>"

import scipy.sparse as sp
from sklearn.preprocessing.discretization import Discretizer
from sklearn.utils.testing import assert_array_almost_equal, assert_array_equal

X = [[0, 0, 1],
     [1, 1, 2],
     [2, 0, 3]]

def test_discretizer_bad_n_bins():
    try:
        dis = Discretizer(n_bins=1)
    except:
        return
    else:
        raise ValueError("Discretizer needs to raise error when number of "
                         "bins is too low.")

def test_discretizer_fit_transform():
    discretized = Discretizer(n_bins=2).fit_transform(X)
    expected = [[0, 0, 0],
                [0, 1, 0],
                [1, 0, 1]]
    assert_array_equal(expected, discretized)

def test_discretizer_fit_transform_sparse():
    sparse_formats = [sp.bsr_matrix, sp.coo_matrix, sp.csc_matrix,
                      sp.csr_matrix, sp.dia_matrix, sp.dok_matrix, sp.lil_matrix]
    for format in sparse_formats:
        sparse_X = format(X)
        discretized = Discretizer(n_bins=2).fit_transform(sparse_X)
        expected = [[0, 0, 0],
                    [0, 1, 0],
                    [1, 0, 1]]
        assert_array_equal(expected, discretized)

def test_discretizer_fit_transform_cat():
    dis = Discretizer(n_bins=3, categorical_features=[0])
    discretized = dis.fit_transform(X)

    expected = [[0, 0, 0],
                [2, 1, 1],
                [0, 2, 2]]

    assert_array_equal(expected, discretized)

def test_discretizer_fit():
    dis = Discretizer(n_bins=3, categorical_features=[1]).fit(X)

    expected = [True, False, True]
    assert_array_almost_equal(expected, dis.continuous_mask_)

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

