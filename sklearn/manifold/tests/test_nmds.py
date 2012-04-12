import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn.manifold import nmds


def test_pva():
    distances = np.array([10., 8, 11, 5, 13, 11, 9, 14, 6, 16])
    similarities = np.arange(10)
    distances_fit = nmds.PVA(distances, similarities)

    assert_array_almost_equal(distances_fit,
                       np.array([8.5, 8.5, 8.5, 8.5, 10.6, 10.6, 10.6, 10.6,
                                 10.6, 16]))
