#! /usr/bin/env python
# Last Change: Fri Oct 20 12:00 PM 2006 J

import copy

import sys
from numpy.testing import *

import numpy as N
from numpy.random import seed

set_package_path()
from pyem import GM, GMM
from pyem.online_em import OnGMM
restore_path()

# #Optional:
# set_local_path()
# # import modules that are located in the same directory as this file.
# restore_path()

# Error precision allowed (nb of decimals)
AR_AS_PREC  = 12

class OnlineEmTest(NumpyTestCase):
    def _create_model(self, d, k, mode, nframes, emiter):
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
        gmm.init(data)

        self.gm0    = copy.copy(gmm.gm)
        # The actual EM, with likelihood computation
        for i in range(emiter):
            g, tgd  = gmm.sufficient_statistics(data)
            gmm.update_em(data, g)

        self.data   = data
        self.gm     = lgm
    
class test_on_off_eq(OnlineEmTest):
    def check_1d(self, level = 1):
        d       = 1
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)
        emiter  = 3

        seed(1)
        self._create_model(d, k, mode, nframes, emiter)
        self._check(d, k, mode, nframes, emiter)

    def check_2d(self, level = 2):
        d       = 2
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)
        emiter  = 3

        seed(1)
        self._create_model(d, k, mode, nframes, emiter)
        self._check(d, k, mode, nframes, emiter)

    def check_5d(self, level = 2):
        d       = 5
        k       = 2
        mode    = 'diag'
        nframes = int(1e2)
        emiter  = 3

        seed(1)
        self._create_model(d, k, mode, nframes, emiter)
        self._check(d, k, mode, nframes, emiter)

    def _check(self, d, k, mode, nframes, emiter):
        #++++++++++++++++++++++++++++++++++++++++
        # Approximate the models with online EM
        #++++++++++++++++++++++++++++++++++++++++
        # Learn the model with Online EM
        ogm         = GM(d, k, mode)
        ogmm        = OnGMM(ogm, 'kmean')
        init_data   = self.data
        ogmm.init(init_data)

        # Check that online kmean init is the same than kmean offline init
        ogm0    = copy.copy(ogm)
        assert_array_equal(ogm0.w, self.gm0.w)
        assert_array_equal(ogm0.mu, self.gm0.mu)
        assert_array_equal(ogm0.va, self.gm0.va)

        # Forgetting param
        lamb	= N.ones((nframes, 1))
        lamb[0] = 0
        nu0		= 1.0
        nu		= N.zeros((len(lamb), 1))
        nu[0]	= nu0
        for i in range(1, len(lamb)):
            nu[i]	= 1./(1 + lamb[i] / nu[i-1])

        # object version of online EM: the p* arguments are updated only at each 
        # epoch, which is equivalent to on full EM iteration on the 
        # classic EM algorithm
        ogmm.pw    = ogmm.cw.copy()
        ogmm.pmu   = ogmm.cmu.copy()
        ogmm.pva   = ogmm.cva.copy()
        for e in range(emiter):
            for t in range(nframes):
                gamma   = ogmm.sufficient_statistics(self.data[t:t+1, :], nu[t])
                ogmm.update_em(self.data[t, :], gamma, nu[t])

            # Change pw args only a each epoch 
            ogmm.pw  = ogmm.cw.copy()
            ogmm.pmu = ogmm.cmu.copy()
            ogmm.pva = ogmm.cva.copy()

        # For equivalence between off and on, we allow a margin of error,
        # because of round-off errors.
        print " Checking precision of equivalence with offline EM trainer "
        maxtestprec = 18
        try :
            for i in range(maxtestprec):
                    assert_array_almost_equal(self.gm.w, ogmm.pw, decimal = i)
                    assert_array_almost_equal(self.gm.mu, ogmm.pmu, decimal = i)
                    assert_array_almost_equal(self.gm.va, ogmm.pva, decimal = i)
            print "\t !! Precision up to %d decimals !! " % i
        except AssertionError:
            if i < AR_AS_PREC:
                print """\t !!NOT OK: Precision up to %d decimals only, 
                    outside the allowed range (%d) !! """ % (i, AR_AS_PREC)
                raise AssertionError
            else:
                print "\t !!OK: Precision up to %d decimals !! " % i

class test_on(OnlineEmTest):
    def check_consistency(self):
        d       = 1
        k       = 2
        mode    = 'diag'
        nframes = int(1e3)
        emiter  = 4

        self._create_model(d, k, mode, nframes, emiter)
        self._run_pure_online(d, k, mode, nframes)
    
    def _run_pure_online(self, d, k, mode, nframes):
        #++++++++++++++++++++++++++++++++++++++++
        # Approximate the models with online EM
        #++++++++++++++++++++++++++++++++++++++++
        ogm     = GM(d, k, mode)
        ogmm    = OnGMM(ogm, 'kmean')
        init_data   = self.data[0:nframes / 20, :]
        ogmm.init(init_data)

        # Forgetting param
        ku		= 0.005
        t0		= 200
        lamb	= 1 - 1/(N.arange(-1, nframes-1) * ku + t0)
        nu0		= 0.2
        nu		= N.zeros((len(lamb), 1))
        nu[0]	= nu0
        for i in range(1, len(lamb)):
            nu[i]	= 1./(1 + lamb[i] / nu[i-1])

        # object version of online EM
        for t in range(nframes):
            # the assert are here to check we do not create copies
            # unvoluntary for parameters
            assert ogmm.pw is ogmm.cw
            assert ogmm.pmu is ogmm.cmu
            assert ogmm.pva is ogmm.cva
            gamma   = ogmm.sufficient_statistics(self.data[t:t+1, :], nu[t])
            ogmm.update_em(self.data[t, :], gamma, nu[t])

        ogmm.gm.set_param(ogmm.cw, ogmm.cmu, ogmm.cva)

if __name__ == "__main__":
    NumpyTest().run()
