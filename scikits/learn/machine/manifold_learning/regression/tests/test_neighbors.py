import numpy as np
from ..neighbors import Neighbors

from numpy.testing import assert_array_equal

def test_neighbors_1D():
    """
    Nearest Neighbors in a line.
    """
    # some constants
    n = 10
    n_2 = n/2
    samples = [[x] for x in range(0, n)]
    labels  = [0]*n_2 + [1]*n_2
    zeros   = np.zeros(n_2)
    ones    = np.ones(n_2)

    # k = 1
    nn = Neighbors(samples, labels=labels, k=1)
    assert_array_equal( nn.predict([ [i +0.01] for i in range(0,n_2)]),
                        zeros)
    assert_array_equal( nn.predict([ [i-0.01] for i in range(n_2, 10)]),
                        ones)

    # k = 3
    nn = Neighbors(samples, labels=labels, k=3)
    assert_array_equal( nn.predict([ [i +0.01] for i in range(0,5)]),
                        zeros)
    assert_array_equal( nn.predict([ [i-0.01] for i in range(5, 10)]),
                        ones)
