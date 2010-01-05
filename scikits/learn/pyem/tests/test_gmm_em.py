#! /usr/bin/env python
# Last Change: Wed Jun 13 07:00 PM 2007 J

# For now, just test that all mode/dim execute correctly

import sys
import os 
from numpy.testing import *

import numpy as N

set_package_path()
from pyem import GMM, GM, EM
restore_path()

set_local_path()
# import modules that are located in the same directory as this file.
from testcommon import DEF_DEC
curpath = sys.path[0]
restore_path()

def load_dataset(filename):
    from scipy.io import loadmat
    dic = loadmat(os.path.join(curpath, filename), squeeze_me = False)
    dic['w0'] = dic['w0'].squeeze()
    dic['w'] = dic['w'].squeeze()
    dic['tw'] = dic['tw'].squeeze()
    return dic

class EmTest(NumpyTestCase):
    def _create_model_and_run_em(self, d, k, mode, nframes):
        #+++++++++++++++++++++++++++++++++++++++++++++++++
        # Generate a model with k components, d dimensions
        #+++++++++++++++++++++++++++++++++++++++++++++++++
        w, mu, va   = GM.gen_param(d, k, mode, spread = 1.5)
        gm          = GM.fromvalues(w, mu, va)
        # Sample nframes frames  from the model
        data        = gm.sample(nframes)

        #++++++++++++++++++++++++++++++++++++++++++
        # Approximate the models with classical EM
        #++++++++++++++++++++++++++++++++++++++++++
        # Init the model
        lgm = GM(d, k, mode)
        gmm = GMM(lgm, 'kmean')

        em  = EM()
        lk  = em.train(data, gmm)

class test_full_run(EmTest):
    """This class only tests whether the algorithms runs. Do not check the
    results."""
    def test_1d(self, level = 1):
        d       = 1
        k       = 2
        mode    = 'full'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def test_2d(self, level = 1):
        d       = 2
        k       = 2
        mode    = 'full'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def test_5d(self, level = 1):
        d       = 5
        k       = 3
        mode    = 'full'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

class test_diag_run(EmTest):
    """This class only tests whether the algorithms runs. Do not test the
    results."""
    def test_1d(self, level = 1):
        d       = 1
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def test_2d(self, level = 1):
        d       = 2
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def test_5d(self, level = 1):
        d       = 5
        k       = 3
        mode    = 'diag'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

class test_datasets(EmTest):
    """This class tests whether the EM algorithms works using pre-computed
    datasets."""
    def test_1d_full(self, level = 1):
        d = 1
        k = 4
        mode = 'full'
        # Data are exactly the same than in diagonal mode, just test that
        # calling full mode works even in 1d, even if it is kind of stupid to
        # do so
        dic = load_dataset('diag_1d_4k.mat')

        gm = GM.fromvalues(dic['w0'], dic['mu0'], dic['va0'])
        gmm = GMM(gm, 'test')
        EM().train(dic['data'], gmm)

        assert_array_almost_equal(gmm.gm.w, dic['w'], DEF_DEC)
        assert_array_almost_equal(gmm.gm.mu, dic['mu'], DEF_DEC)
        assert_array_almost_equal(gmm.gm.va, dic['va'], DEF_DEC)

    def test_1d_diag(self, level = 1):
        d = 1
        k = 4
        mode = 'diag'
        dic = load_dataset('diag_1d_4k.mat')

        gm = GM.fromvalues(dic['w0'], dic['mu0'], dic['va0'])
        gmm = GMM(gm, 'test')
        EM().train(dic['data'], gmm)

        assert_array_equal(gmm.gm.w, dic['w'], DEF_DEC)
        assert_array_equal(gmm.gm.mu, dic['mu'], DEF_DEC)
        assert_array_equal(gmm.gm.va, dic['va'], DEF_DEC)

    def test_2d_full(self, level = 1):
        d = 2
        k = 3
        mode = 'full'
        dic = load_dataset('full_2d_3k.mat')

        gm = GM.fromvalues(dic['w0'], dic['mu0'], dic['va0'])
        gmm = GMM(gm, 'test')
        EM().train(dic['data'], gmm)

        assert_array_equal(gmm.gm.w, dic['w'], DEF_DEC)
        assert_array_equal(gmm.gm.mu, dic['mu'], DEF_DEC)
        assert_array_equal(gmm.gm.va, dic['va'], DEF_DEC)

    def test_2d_diag(self, level = 1):
        d = 2
        k = 3
        mode = 'diag'
        dic = load_dataset('diag_2d_3k.mat')

        gm = GM.fromvalues(dic['w0'], dic['mu0'], dic['va0'])
        gmm = GMM(gm, 'test')
        EM().train(dic['data'], gmm)

        assert_array__almost_equal(gmm.gm.w, dic['w'], DEF_DEC)
        assert_array__almost_equal(gmm.gm.mu, dic['mu'], DEF_DEC)
        assert_array__almost_equal(gmm.gm.va, dic['va'], DEF_DEC)

class test_log_domain(EmTest):
    """This class tests whether the GMM works in log domain."""
    def _test_common(self, d, k, mode):
        dic = load_dataset('%s_%dd_%dk.mat' % (mode, d, k))

        gm = GM.fromvalues(dic['w0'], dic['mu0'], dic['va0'])
        gmm = GMM(gm, 'test')

        a, na = gmm.compute_responsabilities(dic['data'])
        la, nla = gmm.compute_log_responsabilities(dic['data'])

        ta = N.log(a)
        tna = N.log(na)
        if not N.all(N.isfinite(ta)):
            print "precision problem for %s, %dd, %dk, test need fixing" % (mode, d, k)
        else:
            assert_array_almost_equal(ta, la, DEF_DEC)

        if not N.all(N.isfinite(tna)):
            print "precision problem for %s, %dd, %dk, test need fixing" % (mode, d, k)
        else:
            assert_array_almost_equal(tna, nla, DEF_DEC)

    def test_2d_diag(self, level = 1):
        d = 2
        k = 3
        mode = 'diag'
        self._test_common(d, k, mode)

    def test_1d_full(self, level = 1):
        d = 1
        k = 4
        mode = 'diag'
        self._test_common(d, k, mode)

    def test_2d_full(self, level = 1):
        d = 2
        k = 3
        mode = 'full'
        self._test_common(d, k, mode)

if __name__ == "__main__":
    NumpyTest().run()
