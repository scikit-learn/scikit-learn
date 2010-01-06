import os

import numpy as np
from numpy.testing import *
from scipy.io import loadmat

from scikits.learn.machine.em2.gmm import Parameters, EM

def load_dataset(filename):
    from scipy.io import loadmat
    dic = loadmat(os.path.join(os.path.dirname(__file__), filename), 
                  squeeze_me=False)
    dic['w0'] = dic['w0'].squeeze()
    dic['w'] = dic['w'].squeeze()
    dic['tw'] = dic['tw'].squeeze()
    return dic

class BasicEmTest(TestCase):
    """This class tests whether the EM algorithms works using pre-computed
    datasets."""
    def _test(self, dataset):
        dic = load_dataset(dataset)

        params = Parameters.fromvalues(dic['w0'], dic['mu0'], dic['va0'])
        em = EM()
        em.train(dic['data'], params)

        assert_array_almost_equal(params.w, dic['w'])
        assert_array_almost_equal(params.mu, dic['mu'])
        assert_array_almost_equal(params.va, dic['va'])

    def test_2d_diag(self, level = 1):
        filename = 'diag_2d_3k.mat'
        self._test(filename)

if __name__ == '__main__':
    run_module_suite()
