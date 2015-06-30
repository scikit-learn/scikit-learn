"""Test the rcv1 loader.

Skipped if rcv1 is not already downloaded to data_home.
"""

import errno
import scipy.sparse as sp
import numpy as np
from sklearn.datasets import fetch_rcv1
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import SkipTest


def test_fetch_rcv1():
    try:
        data1 = fetch_rcv1(shuffle=False, download_if_missing=False)
    except IOError as e:
        if e.errno == errno.ENOENT:
            raise SkipTest("Download RCV1 dataset to run this test.")

    X1, Y1 = data1.data, data1.target
    cat_list, s1 = data1.target_names.tolist(), data1.sample_id

    # test sparsity
    assert_true(sp.issparse(X1))
    assert_true(sp.issparse(Y1))
    assert_equal(60915113, X1.data.size)
    assert_equal(2606875, Y1.data.size)

    # test shapes
    assert_equal((804414, 47236), X1.shape)
    assert_equal((804414, 103), Y1.shape)
    assert_equal((804414,), s1.shape)
    assert_equal(103, len(cat_list))

    # test number of sample for some categories
    some_categories = ('GMIL', 'E143', 'CCAT')
    number_non_zero_in_cat = (5, 1206, 381327)
    for num, cat in zip(number_non_zero_in_cat, some_categories):
        j = cat_list.index(cat)
        assert_equal(num, Y1[:, j].data.size)

    # test shuffling
    data2 = fetch_rcv1(shuffle=True, random_state=77,
                       download_if_missing=False)
    X2, Y2 = data2.data, data2.target
    s2 = data2.sample_id

    assert_true((s1 != s2).any())
    assert_array_equal(np.sort(s1), np.sort(s2))

    # test some precise values
    some_sample_id = (2286, 333274, 810593)
    # indice of first nonzero feature
    indices = (863, 863, 814)
    # value of first nonzero feature
    feature = (0.04973993, 0.12272136, 0.14245221)
    for i, j, v in zip(some_sample_id, indices, feature):
        i1 = np.nonzero(s1 == i)[0][0]
        i2 = np.nonzero(s2 == i)[0][0]

        sp_1 = X1[i1].sorted_indices()
        sp_2 = X2[i2].sorted_indices()

        assert_almost_equal(sp_1[0, j], v)
        assert_almost_equal(sp_2[0, j], v)

        assert_array_equal(np.sort(Y1[i1].indices), np.sort(Y2[i2].indices))
