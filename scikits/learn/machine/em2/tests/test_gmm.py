import os

import numpy as np
from numpy.testing import *
from scipy.io import loadmat

from scikits.learn.machine.em2.gmm import Parameters, EM
from scikits.learn.machine.em2.gm import GM
from scikits.learn.machine.em2.gmm import logresp as plogresp

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

class TestLogresp(TestCase):
    def ref(self, x, w, mu, va):
        k, d = mu.shape
        tmp = np.ones((x.shape[0], k))
        for c in range(k):
            inva = 1/va[c]
            fac = (2 * np.pi) ** (-d/2.) * np.prod(np.sqrt(inva))
            tmp[:, c] = np.log(fac) + np.dot((x-mu[c]) ** 2, -0.5 * inva.T)

        tmp += np.log(w)
        return tmp - np.log(np.sum(np.exp(tmp), axis=-1))[:, np.newaxis]

    def _test(self, d, k, mode):
        from scikits.learn.machine.em2.gmm import logresp
        x = np.random.randn(100, d)
        w, mu, va = GM.genparams(d, k, mode)
        yr = self.ref(x, w, mu, va)

        y = plogresp(x, w, mu, va)
        assert_array_almost_equal(y, yr)

    def test_1d_1k(self):
        return self._test(1, 1, 'diag')

    def test_1d_2k(self):
        return self._test(1, 2, 'diag')

    def test_1d_10k(self):
        return self._test(1, 10, 'diag')

if __name__ == '__main__':
    run_module_suite()
