# /usr/bin/python
# Last Change: Mon Aug 28 05:00 PM 2006 J

# TODO:
#   - which methods to avoid va shrinking to 0 ?
#   - online EM

import numpy as N
import numpy.linalg as lin
from numpy.random import randn
#import _c_densities as densities
import densities
from kmean import kmean
from gauss_mix import GM

# Error classes
class GmmError(Exception):
    """Base class for exceptions in this module."""
    pass

class GmmParamError:
    """Exception raised for errors in gmm params

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message    = message
    
    def __str__(self):
        return self.message

# Not sure yet about how to design different mixture models. Most of the code 
# is different # (pdf, update part of EM, etc...) and I am not sure it makes 
# sense to use inheritance for # interface specification in python, since its 
# dynamic type systeme.

# Anyway, a mixture class should encapsulates all details concerning a mixture model:
#   - internal parameters for the pdfs
#   - can compute sufficient statistics for EM
#   - can sample a model
#   - can generate random valid parameters for a new pdf (using class method)
class MixtureModel:
    pass

class ExpMixtureModel(MixtureModel):
    """Class to model mixture of exponential pdf (eg Gaussian, exponential, Laplace, 
    etc..). This is a special case because some parts of EM are common for those
    models..."""
    pass

class GMM(ExpMixtureModel):
    """ A class to model a Gaussian Mixture Model (GMM). An instance of 
    this class is created by giving weights, mean and variances in the ctor.
    An instanciated object can be sampled, trained by EM. 
    
    The class method gen_model can be used without instanciation."""

    def init_kmean(self, data, niter = 5):
        """ Init the model with kmean."""
        k       = self.gm.k
        d       = self.gm.d
        init    = data[0:k, :]

        (code, label)   = kmean(data, init, niter)

        w   = N.ones(k) / k
        mu  = code.copy()
        if self.gm.mode == 'diag':
            va = N.zeros((k, d))
            for i in range(k):
                for j in range(d):
                    va[i,j] = N.cov(data[N.where(label==i), j], rowvar = 0)
        elif self.gm.mode == 'full':
            va  = N.zeros((k*d, d))
            for i in range(k):
                va[i*d:i*d+d,:] = \
                    N.cov(data[N.where(label==i)], rowvar = 0)
        else:
            raise GmmParamError("mode " + str(mode) + " not recognized")

        self.gm.w   = w
        self.gm.mu  = mu
        self.gm.va  = va

    def init_random(self, data):
        """ Init the model at random."""
        k   = self.gm.k
        d   = self.gm.d
        if mode == 'diag':
            w   = N.ones(k) / k
            mu  = randn(k, d)
            va  = N.fabs(randn(k, d))
        else:
            raise GmmParamError("""init_random not implemented for
                    mode %s yet""", mode)

        self.gm.w   = w
        self.gm.mu  = mu
        self.gm.va  = va

    # TODO: 
    #   - format of parameters ? For variances, list of variances matrix,
    #   keep the current format, have 3d matrices ?
    #   - To handle the different modes, we could do something "fancy" such as
    #   replacing methods, to avoid checking cases everywhere and unconsistency.
    def __init__(self, gm, init = 'kmean'):
        """ Initialize a GMM with weight w, mean mu and variances va, and initialization
        method for training init (kmean by default)"""
        self.gm = gm

        # Possible init methods
        init_methods = {'kmean': self.init_kmean, 'random' : self.init_random}

        if init not in init_methods:
            raise GmmParamError('init method %s not recognized' + str(init))

        self.init   = init_methods[init]

    def sufficient_statistics(self, data):
        """ Return normalized and non-normalized sufficient statistics
        from the model.
        
        Computes the latent variable distribution (a 
        posteriori probability) knowing the explicit data 
        for the Gaussian model (w, mu, var): gamma(t, i) = 
            P[state = i | observation = data(t); w, mu, va]

        This is basically the E step of EM for GMM."""
        n   = data.shape[0]

        # compute the gaussian pdf
        tgd	= multiple_gauss_den(data, self.gm.mu, self.gm.va)
        # multiply by the weight
        tgd	*= self.gm.w
        # Normalize to get a pdf
        gd	= tgd  / N.sum(tgd, axis=1)[:, N.newaxis]

        return gd, tgd

    def update_em(self, data, gamma):
        """Computes update of the Gaussian Mixture Model (M step)
        from the a posteriori pdf, computed by gmm_posterior
        (E step).
        """
        k       = self.gm.k
        d       = self.gm.d
        n       = data.shape[0]
        invn    = 1.0/n
        mGamma  = N.sum(gamma, axis = 0)

        if self.gm.mode == 'diag':
            mu  = N.zeros((k, d))
            va  = N.zeros((k, d))
            gamma   = gamma.transpose()
            for c in range(k):
                x   = N.dot(gamma[c:c+1,:], data)[0,:]
                xx  = N.dot(gamma[c:c+1,:], data ** 2)[0,:]

                mu[c,:] = x / mGamma[c]
                va[c,:] = xx  / mGamma[c] - mu[c,:] ** 2
            w   = invn * mGamma

        elif self.gm.mode == 'full':
            mu  = N.zeros((k, d))
            va  = N.zeros((k*d, d))

            for c in range(k):
                x   = N.sum(N.outerproduct(gamma[:, c], 
                            N.ones((1, d))) * data, axis = 0)
                xx  = N.zeros((d, d))
                
                # This should be much faster than recursing on n...
                for i in range(d):
                    for j in range(d):
                        xx[i,j] = N.sum(data[:,i] * data[:,j] * gamma[:,c], axis = 0)

                mu[c,:] = x / mGamma[c]
                va[c*d:c*d+d,:] = xx  / mGamma[c] - \
                                    N.outerproduct(mu[c,:], mu[c,:])
            w   = invn * mGamma
        else:
            raise GmmParamError("varmode not recognized")

        self.gm.set_param(w, mu, va)

class EM:
    """An EM trainer. An EM trainer
    trains from data, with a model
    
    Not really useful yet"""
    def __init__(self):
        pass

    def train(self, data, model, niter = 10):
        """
        Train a model using data, with niter iterations.

        Args:
            - data:     contains the observed features, one row is one frame, ie one 
            observation of dimension d
            - model:    object of class Mixture
            - niter:    number of iterations

        The model is trained, and its parameters updated accordingly.

        Returns:
            likelihood (one value per iteration).
        """

        # Initialize the data (may do nothing depending on the model)
        model.init(data)

        # Likelihood is kept
        like    = N.zeros(niter)

        # Em computation, with computation of the likelihood
        for i in range(niter):
            g, tgd      = model.sufficient_statistics(data)
            like[i]     = N.sum(N.log(N.sum(tgd, 1)), axis = 0)
            model.update_em(data, g)

        return like
    
class OnlineEM:
    "An online EM trainer. "
    def __init__(self):
        raise GmmError("not implemented yet")

# Misc functions
def multiple_gauss_den(data, mu, va):
    """Helper function to generate several Gaussian
    pdf (different parameters) from the same data"""
    mu  = N.atleast_2d(mu)
    va  = N.atleast_2d(va)

    K   = mu.shape[0]
    n   = data.shape[0]
    d   = data.shape[1]
    
    y   = N.zeros((n, K))
    if mu.size == va.size:
        for i in range(K):
            y[:, i] = densities.gauss_den(data, mu[i, :], va[i, :])
    else:
        for i in range(K):
            y[:, i] = densities.gauss_den(data, mu[i, :], 
                        va[d*i:d*i+d, :])

    return y

if __name__ == "__main__":
    import copy
    #=============================
    # Simple GMM with 5 components
    #=============================

    #+++++++++++++++++++++++++++++
    # Meta parameters of the model
    #   - k: Number of components
    #   - d: dimension of each Gaussian
    #   - mode: Mode of covariance matrix: full or diag
    #   - nframes: number of frames (frame = one data point = one
    #   row of d elements
    k       = 3 
    d       = 2         
    mode    = 'diag'        
    nframes = 1e3

    #+++++++++++++++++++++++++++++++++++++++++++
    # Create an artificial GMM model, samples it
    #+++++++++++++++++++++++++++++++++++++++++++
    print "Generating the mixture"
    # Generate a model with k components, d dimensions
    w, mu, va   = GM.gen_param(d, k, mode, spread = 3)
    gm          = GM(d, k, mode)
    gm.set_param(w, mu, va)

    # Sample nframes frames  from the model
    data    = gm.sample(nframes)

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
    niter   = 10
    like    = N.zeros(niter)

    print "computing..."
    for i in range(niter):
        g, tgd  = gmm.sufficient_statistics(data)
        like[i] = N.sum(N.log(N.sum(tgd, 1)), axis = 0)
        gmm.update_em(data, g)
    # # Alternative form, by using EM class: as the EM class
    # # is quite rudimentary now, it is not very useful, just save
    # # a few lines
    # em      = EM()
    # like    = em.train(data, gmm, niter)

    #+++++++++++++++
    # Draw the model
    #+++++++++++++++
    print "drawing..."
    import pylab as P
    P.subplot(2, 1, 1)
    # Draw what is happening
    P.plot(data[:, 0], data[:, 1], '.', label = '_nolegend_')

    # Real confidence ellipses
    Xre, Yre  = gm.conf_ellipses()
    P.plot(Xre[0], Yre[0], 'g', label = 'true confidence ellipsoides')
    for i in range(1,k):
        P.plot(Xre[i], Yre[i], 'g', label = '_nolegend_')

    # Initial confidence ellipses as found by kmean
    X0e, Y0e  = gm0.conf_ellipses()
    P.plot(X0e[0], Y0e[0], 'k', label = 'initial confidence ellipsoides')
    for i in range(1,k):
        P.plot(X0e[i], Y0e[i], 'k', label = '_nolegend_')

    # Values found by EM
    Xe, Ye  = lgm.conf_ellipses()
    P.plot(Xe[0], Ye[0], 'r', label = 'confidence ellipsoides found by EM')
    for i in range(1,k):
        P.plot(Xe[i], Ye[i], 'r', label = '_nolegend_')
    P.legend(loc = 0)

    #from scipy.cluster.vq import kmeans
    #code    = kmeans(data, k)[0]
    #print code
    #P.plot(code[:,0], code[:, 1], 'oy')
    #P.plot(gm0.mu[:,0], gm0.mu[:, 1], 'ok')
    P.subplot(2, 1, 2)
    P.plot(like)
    P.title('log likelihood')

    # #++++++++++++++++++
    # # Export the figure
    # #++++++++++++++++++
    # F   = P.gcf()
    # DPI = F.get_dpi()
    # DefaultSize = F.get_size_inches()
    # # the default is 100dpi for savefig:
    # F.savefig("example1.png")

    # # Now make the image twice as big, while keeping the fonts and all the
    # # same size
    # F.set_figsize_inches( (DefaultSize[0]*2, DefaultSize[1]*2) )
    # Size = F.get_size_inches()
    # print "Size in Inches", Size
    # F.savefig("example2.png")
    P.show()
