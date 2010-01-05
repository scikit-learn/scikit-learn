# /usr/bin/python
# Last Change: Mon Oct 16 03:00 PM 2006 J

#---------------------------------------------
# This is not meant to be used yet !!!! I am 
# not sure how to integrate this stuff inside
# the package yet. The cases are:
#   - we have a set of data, and we want to test online EM 
#   compared to normal EM 
#   - we do not have all the data before putting them in online EM:
#   eg current frame depends on previous frame in some way.

import numpy as N
from numpy import mean
from numpy.testing import assert_array_almost_equal, assert_array_equal

from gmm_em import ExpMixtureModel, GMM, EM, multiple_gauss_den
from gauss_mix import GM
from kmean import kmean

print "======================================================"
import traceback
f = traceback.extract_stack()
print f[0][0] + " This is beta code, don't use it !        "
print "======================================================"

# Error classes
class OnGmmError(Exception):
    """Base class for exceptions in this module."""
    pass

class OnGmmParamError:
    """Exception raised for errors in gmm params

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message    = message
    
    def __str__(self):
        return self.message

class OnGMM(ExpMixtureModel):
    def init_kmean(self, init_data, niter = 5):
        """ Init the model at random."""
        k   = self.gm.k
        d   = self.gm.d
        if mode == 'diag':
            w           = N.ones(k) / k

            # Init the internal state of EM
            self.cx     = N.outer(w, mean(init_data, 0))
            self.cxx    = N.outer(w, mean(init_data ** 2, 0))

            # w, mu and va init is the same that in the standard case
            (code, label)   = kmean(init_data, init_data[0:k, :], niter)
            mu          = code.copy()
            va          = N.zeros((k, d))
            for i in range(k):
                for j in range(d):
                    va [i,j] = N.cov(init_data[N.where(label==i), j], rowvar = 0)
        else:
            raise OnGmmParamError("""init_online not implemented for
                    mode %s yet""", mode)

        self.gm.set_param(w, mu, va)
        self.cw     = w
        self.cmu    = mu
        self.cva    = va

        self.pw     = self.cw.copy()
        self.pmu    = self.cmu.copy()
        self.pva    = self.cva.copy()

    def __init__(self, gm, init_data, init = 'kmean'):
        self.gm = gm
        
        # Possible init methods
        init_methods = {'kmean' : self.init_kmean}

        self.init   = init_methods[init]

    def sufficient_statistics(self, frame, nu):
        """ sufficient_statistics(frame, nu)
        
        frame has to be rank 2 !"""
        gamma   = multiple_gauss_den(frame, self.pmu, self.pva)[0]
        gamma   *= self.pw
        gamma   /= N.sum(gamma)
        # <1>(t) = cw(t), each column is one component cw = (cw1, ..., cwK);
        self.cw	= (1 - nu) * self.cw + nu * gamma

        return gamma

    def update_em(self, frame, gamma, nu):
        for k in range(self.gm.k):
            self.cx[k]   = (1 - nu) * self.cx[k] + nu * frame * gamma[k]
            self.cxx[k]  = (1 - nu) * self.cxx[k] + nu * frame ** 2 * gamma[k]

            self.cmu[k]  = self.cx[k] / self.cw[k]
            self.cva[k]  = self.cxx[k] / self.cw[k] - self.cmu[k] ** 2
    
if __name__ == "__main__":
    import copy
    #=============================
    # Simple GMM with 2 components
    #=============================

    #+++++++++++++++++++++++++++++
    # Meta parameters of the model
    #   - k: Number of components
    #   - d: dimension of each Gaussian
    #   - mode: Mode of covariance matrix: full or diag
    #   - nframes: number of frames (frame = one data point = one
    #   row of d elements
    k       = 2 
    d       = 1
    mode    = 'diag'
    nframes = int(1e3)
    emiter  = 10

    #+++++++++++++++++++++++++++++++++++++++++++
    # Create an artificial GMM model, samples it
    #+++++++++++++++++++++++++++++++++++++++++++
    print "Generating the mixture"
    # Generate a model with k components, d dimensions
    w, mu, va   = GM.gen_param(d, k, mode, spread = 1.5)
    gm          = GM.fromvalues(w, mu, va)

    # Sample nframes frames  from the model
    data    = gm.sample(nframes, )

    #++++++++++++++++++++++++
    # Learn the model with EM
    #++++++++++++++++++++++++

    # Init the model
    print "Init a model for learning, with kmean for initialization"
    lgm = GM(d, k, mode)
    gmm = GMM(lgm, 'kmean')
    gmm.init(data)

    # Keep the initialized model for drawing
    gm0 = copy.copy(lgm)

    # The actual EM, with likelihood computation
    like    = N.zeros(emiter)
    print "computing..."
    for i in range(emiter):
        g, tgd  = gmm.sufficient_statistics(data)
        like[i] = N.sum(N.log(N.sum(tgd, 1)), axis = 0)
        gmm.update_em(data, g)

    #+++++++++++++++++++++++++++++++
    # Learn the model with Online EM
    #+++++++++++++++++++++++++++++++
    ogm     = GM(d, k, mode)
    ogmm    = OnGMM(ogm, 'kmean')
    #init_data   = data[0:nframes / 20, :]
    init_data   = data
    ogmm.init(init_data)

    ogm2    = GM(d, k, mode)
    ogmm2   = OnGMM(ogm2, 'kmean')
    #init_data   = data[0:nframes / 20, :]
    init_data   = data
    ogmm2.init(init_data)

    ogm0    = copy.copy(ogm)
    assert_array_equal(ogm0.w, gm0.w)
    assert_array_equal(ogm0.mu, gm0.mu)
    assert_array_equal(ogm0.va, gm0.va)

    ogm02   = copy.copy(ogm2)
    assert_array_equal(ogm02.w, gm0.w)
    assert_array_equal(ogm02.mu, gm0.mu)
    assert_array_equal(ogm02.va, gm0.va)

    w0  = gm0.w.copy()
    mu0 = gm0.mu.copy()
    va0 = gm0.va.copy()

    cx  = ogmm.cx
    cxx = ogmm.cxx
    
    cw  = w0.copy()
    cmu = mu0.copy()
    cva = va0.copy()
    
    pw  = cw.copy()
    pmu = cmu.copy()
    pva = cva.copy()

    # Forgetting param
    ku		= 0.005
    t0		= 200
    lamb	= N.ones((nframes, 1))
    lamb[0] = 0
    nu0		= 1.0
    #lamb	= 1 - 1/(N.arange(-1, nframes-1) * ku + t0)
    #nu0		= 0.2
    nu		= N.zeros((len(lamb), 1))
    nu[0]	= nu0
    for i in range(1, len(lamb)):
        nu[i]	= 1./(1 + lamb[i] / nu[i-1])

    gamma   = N.zeros((nframes, k))
    agamma  = []
    apw     = []
    apmu    = []
    apva    = []
    print "========== Online Manual ==========="
    for e in range(emiter):
        print "online manual: estep %d, printing p* state " % e
        apw.append(pw.copy())
        apmu.append(pmu.copy())
        apva.append(pva.copy())
        for t in range(nframes):
            gamma[t]    = multiple_gauss_den(data[t:t+1, :], pmu, pva)[0]
            gamma[t]    *= pw
            gamma[t]    /= N.sum(gamma[t])
            # <1>(t) = cw(t), each column is one component cw = (cw1, ..., cwK);
            # <x>(t) = cx(t)
            cw	= (1 - nu[t]) * cw + nu[t] * gamma[t]
            # loop through each component
            for i in range(k):
                cx[i]   = (1 - nu[t]) * cx[i] + nu[t] * data[t, :] * gamma[t, i]
                cxx[i]  = (1 - nu[t]) * cxx[i] + nu[t] * data[t, :] ** 2 * gamma[t, i]

                cmu[i]  = cx[i] / cw[i]
                cva[i]  = cxx[i] / cw[i] - cmu[i] ** 2

            #pw  = cw.copy()
            #pmu = cmu.copy()
            #pva = cva.copy()
        print "gamma[end]: " + str(gamma[-1])
        pw  = cw.copy()
        pmu = cmu.copy()
        pva = cva.copy()
        agamma.append(gamma.copy())

    gamma2  = N.zeros((nframes, k))
    agamma2 = []
    apw2    = []
    apmu2   = []
    apva2   = []
    print "========== Online Automatic ==========="
    for e in range(emiter):
        print "online automatic: estep %d, printing p* state " % e
        apw2.append(ogmm2.pw.copy())
        apmu2.append(ogmm2.pmu.copy())
        apva2.append(ogmm2.pva.copy())
        for t in range(nframes):
            gamma2[t]   = ogmm2.sufficient_statistics(data[t:t+1, :], nu[t])
            #gamma2[t]   = multiple_gauss_den(data[t:t+1, :], ogmm2.pmu, ogmm2.pva)[0]
            #gamma2[t]   *= ogmm2.pw
            #gamma2[t]   /= N.sum(gamma2[t])
            #try:
            #    assert_array_equal(agamma, gamma2[t])
            #except AssertionError:
            #    print "error for estep %d, step %d" % (e, t)
            #    print ogmm2.pw
            #    print ogmm2.pmu
            #    print ogmm2.pva
            #    raise 
            ogmm2.update_em(data[t, :], gamma2[t], nu[t])
            #ogmm2.cw	= (1 - nu[t]) * ogmm2.cw + nu[t] * agamma
            ## loop through each component
            #for i in range(k):
            #    ogmm2.cx[i]   = (1 - nu[t]) * ogmm2.cx[i] + nu[t] * data[t, :] * agamma[i]
            #    ogmm2.cxx[i]  = (1 - nu[t]) * ogmm2.cxx[i] + nu[t] * data[t, :] ** 2 * agamma[i]

            #    ogmm2.cmu[i]  = ogmm2.cx[i] / ogmm2.cw[i]
            #    ogmm2.cva[i]  = ogmm2.cxx[i] / ogmm2.cw[i] - ogmm2.cmu[i] ** 2

        print "gamma[end]: " + str(gamma2[-1])
        agamma2.append(gamma2.copy())
        ogmm2.pw  = ogmm2.cw.copy()
        ogmm2.pmu = ogmm2.cmu.copy()
        ogmm2.pva = ogmm2.cva.copy()

    # #ogm.set_param(pw, pmu, pva)
    # print "weights off vs on: \n" + str(lgm.w) + "\n " + str(cw)
    # print "mean off vs on: \n" + str(lgm.mu) + "\n " + str(cmu)
    # print "variances off vs on: \n" + str(lgm.va) + "\n " + str(cva)
    # print "weights off vs on2: \n" + str(lgm.w) + "\n " + str(ogmm2.cw)
    # print "mean off vs on2: \n" + str(lgm.mu) + "\n " + str(ogmm2.cmu)
    # print "variances off vs on2: \n" + str(lgm.va) + "\n " + str(ogmm2.cva)
    # assert_array_almost_equal(cw, lgm.w)
    # assert_array_almost_equal(cmu, lgm.mu)
    # assert_array_almost_equal(cva, lgm.va)
    assert_array_equal(ogmm2.pw, pw)
    assert_array_equal(ogmm2.pmu, pmu)
    assert_array_equal(ogmm2.pva, pva)
    assert_array_equal(agamma[0], agamma2[0])
    #assert_array_almost_equal(ogmm2.cw, lgm.w)
    #assert_array_almost_equal(ogmm2.cmu, lgm.mu)
    #assert_array_almost_equal(ogmm2.cva, lgm.va)
    # #+++++++++++++++
    # # Draw the model
    # #+++++++++++++++
    # print "drawing..."
    # import pylab as P
    # P.subplot(2, 1, 1)

    # if d == 1:
    #     # Real confidence ellipses
    #     h   = gm.plot1d()
    #     [i.set_color('g') for i in h['pdf']]
    #     h['pdf'][0].set_label('true pdf')

    #     # Initial confidence ellipses as found by kmean
    #     h0  = gm0.plot1d()
    #     [i.set_color('k') for i in h0['pdf']]
    #     h0['pdf'][0].set_label('initial pdf')

    #     # Values found by EM
    #     hl  = lgm.plot1d(fill = 1, level = 0.66)
    #     [i.set_color('r') for i in hl['pdf']]
    #     hl['pdf'][0].set_label('pdf found by EM')

    #     # Values found by Online EM
    #     ho  = ogm.plot1d(fill = 1, level = 0.66)
    #     [i.set_color('b') for i in ho['pdf']]
    #     ho['pdf'][0].set_label('pdf found by online EM')

    #     P.legend(loc = 0)
    # else:
    #     # Draw what is happening
    #     P.plot(data[:, 0], data[:, 1], '.', label = '_nolegend_')

    #     # Real confidence ellipses
    #     Xre, Yre  = gm.conf_ellipses()
    #     P.plot(Xre[0], Yre[0], 'g', label = 'true confidence ellipsoides')
    #     for i in range(1,k):
    #         P.plot(Xre[i], Yre[i], 'g', label = '_nolegend_')

    #     # Initial confidence ellipses as found by kmean
    #     X0e, Y0e  = gm0.conf_ellipses()
    #     P.plot(X0e[0], Y0e[0], 'k', label = 'initial confidence ellipsoides')
    #     for i in range(1,k):
    #         P.plot(X0e[i], Y0e[i], 'k', label = '_nolegend_')

    #     # Values found by EM
    #     Xe, Ye  = lgm.conf_ellipses()
    #     P.plot(Xe[0], Ye[0], 'r', label = 'confidence ellipsoides found by EM')
    #     for i in range(1,k):
    #         P.plot(Xe[i], Ye[i], 'r', label = '_nolegend_')
    #     P.legend(loc = 0)

    #     # Values found by Online EM
    #     Xe, Ye  = ogm.conf_ellipses()
    #     P.plot(Xe[0], Ye[0], 'm', label = 'confidence ellipsoides found by Online EM')
    #     for i in range(1,k):
    #         P.plot(Xe[i], Ye[i], 'm', label = '_nolegend_')
    #     P.legend(loc = 0)


    # P.subplot(2, 1, 2)
    # P.plot(like)
    # P.title('log likelihood')

    # P.show()
