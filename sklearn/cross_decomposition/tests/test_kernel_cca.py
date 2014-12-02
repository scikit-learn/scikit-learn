import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.cross_decomposition import KernelCCA


def test_kernel_cca():
    # test against matlab implementation by David Hardoon
    X = np.array([[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [3.,5.,4.]])
    Y = np.array([[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]])
    kcca = KernelCCA(kernel="linear",n_components=2,
                     kapa=0.1,eta=0.1, pgso=True)
    kcca.fit(X, Y)
    matlab_lambdas = np.array([0.9998, 0.7698])
    matlab_alphas = np.array(
        [[-0.0868,  1.3628],
         [ 0.0181,  0.1309],
         [-0.0090,  0.5048],
         [ 0.0281, -0.3625]])
    matlab_betas = np.array(
        [[-0.0207, -3.2597],
         [ 0.0128,  1.9876],
         [-0.0269, -4.4307],
         [ 0.0154,  2.0363]])
    assert_array_equal(matlab_lambdas, np.round(kcca.lambdas_, decimals=4))
    assert_array_equal(matlab_alphas, np.round(kcca.alphas_, decimals=4))
