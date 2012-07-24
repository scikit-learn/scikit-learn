from sklearn.naive_bayes import MultinomialNB
from sklearn.semi_supervised import OneDNFTransformer

import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_equal


X = [[4, 5, 1, 0, 0],
     [6, 7, 1, 0, 0],
     [0, 0, 9, 0, 0],
     [0, 0, 0, 1, 6],
     [1, 0, 0, 4, 5]]

y = [1, 1, -1, -1, -1]


def test_1dnf():
    onednf = OneDNFTransformer()
    y_l = onednf.fit_transform(X, y)
    assert_equal(onednf.pos_support_, [0, 1])
    assert_equal(y_l, [1, 1, 0, 0, -1])


def test_sparse_1dnf():
    onednf = OneDNFTransformer()
    y_l = onednf.fit_transform(sp.lil_matrix(X), y)
    assert_equal(onednf.pos_support_, [0, 1])
    assert_equal(y_l, [1, 1, 0, 0, -1])
