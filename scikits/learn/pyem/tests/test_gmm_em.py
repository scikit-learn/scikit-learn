#! /usr/bin/env python
# Last Change: Tue Oct 24 06:00 PM 2006 J

# For now, just test that all mode/dim execute correctly

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from pyem import GMM, GM, EM
restore_path()

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

class test_full(EmTest):
    def check_1d(self, level = 1):
        d       = 1
        k       = 2
        mode    = 'full'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def check_2d(self, level = 1):
        d       = 2
        k       = 2
        mode    = 'full'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def check_5d(self, level = 1):
        d       = 5
        k       = 3
        mode    = 'full'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

class test_diag(EmTest):
    def check_1d(self, level = 1):
        d       = 1
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def check_2d(self, level = 1):
        d       = 2
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

    def check_5d(self, level = 1):
        d       = 5
        k       = 3
        mode    = 'diag'
        nframes = int(1e2)

        #seed(1)
        self._create_model_and_run_em(d, k, mode, nframes)

if __name__ == "__main__":
    NumpyTest().run()
