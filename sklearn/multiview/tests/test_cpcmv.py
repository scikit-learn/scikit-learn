import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
import sklearn.multiview.cpcmv as cpcmv


def test_cpc():
    data = np.arange(50, dtype=float).reshape((2, 5, 5))
    eigenvalues, eigenvectors = cpcmv.cpc(data, k=2)

    real_eigenvalues = np.array([[66.60954262, 176.89090037],
                                 [-6.60954262, 8.10909963]])
    real_eigenvectors = np.array([[0.20303518, 0.74751369],
                                  [0.31154745, 0.45048661],
                                  [0.42005972, 0.15345953],
                                  [0.528572, -0.14356755],
                                  [0.63708427, -0.44059463]])

    assert_array_almost_equal(eigenvalues, real_eigenvalues, decimal=4)
    assert_array_almost_equal(eigenvectors, real_eigenvectors, decimal=4)


def test_cpc_error():
    # k value cannot be negative
    data = np.arange(40, dtype=float).reshape((2, 5, 4))
    cpc_est = cpcmv.MVCPC(k=-2)
    assert_raises(ValueError, cpc_est.fit, data)

    # Second and third dimensions must be equal
    cpc_est = cpcmv.MVCPC(k=2)
    assert_raises(ValueError, cpc_est.fit, data)


def test_cpcmv():
    data = np.arange(50, dtype=float).reshape((2, 5, 5))

    cpc_est = cpcmv.MVCPC(k=2)
    cpc_est.fit(data)
