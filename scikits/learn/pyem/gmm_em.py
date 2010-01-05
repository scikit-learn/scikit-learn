# /usr/bin/python
# Last Change: Thu May 25 03:00 PM 2006 J

# TODO:
#   - interface with libem
#   - implements E and M in batch mode for mixtures, with special case for gmm 
#   ( general EM should be outside this module. This one should be GMM specific, maybe)
#   - implements E and M in online mode
#   - python has list. Does it make sense to keep the h2m format for mean, va, when
#   we can use list instead ?
#   - a Gaussian Mixture Model should be a class, really !

import numpy as N
import numpy.linalg as lin
import densities

MAX_DEV = 1e-10

# Error classes
class GmmError(Exception):
    """Base class for exceptions in this module."""
    pass

class GmmParamError(GmmError):
    """Exception raised for errors in gmm params

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message    = message
    
    def __str__(self):
        return self.message


# Function to generate a GMM, or valid parameters for GMM
def gen_rand_index(p, n):
    """Generate a N samples vector containing random index between 1 
    and length(p), each index i with probability p(i)"""
    # TODO Check args here
    
    # TODO: check each value of inverse distribution is
    # different
    invcdf  = N.cumsum(p)
    uni     = N.rand(n)
    index   = N.zeros(n)

    # This one should be a bit faster
    for k in range(len(p)-1, 0, -1):
        blop        = N.where(N.logical_and(invcdf[k-1] <= uni, 
                    uni < invcdf[k]))
        index[blop] = k
        
    return index

def gen_gmm(w, mu, va, n):
    """Generate a gaussiam mixture model with weights w, 
    mean mu and variances va. Each column of the parameters
    are one component parameter.
    """
    # Check args
    K, d, mode  = check_gmm_param(w, mu, va)

    # Generate the mixture
    S   = gen_rand_index(w, n)  # State index (ie hidden var)
    X   = N.randn(n, d)         # standard gaussian

    if mode == 'diag':
        X   = mu[S, :]  + X * N.sqrt(va[S,:])
    elif mode == 'full':
        # Naive implementation
        # cho = N.zeros(va.shape, float)
        # for k in range(K):
        #     # Using cholesky is more stable than sqrtm; sqrtm is not
        #     # available in numpy anyway, only in scipy...
        #     cho[k*d:k*d+d,:]    = lin.cholesky(va[k*d:k*d+d,:])
        # for i in range(n):
        # X[i,:]  = mu[S[i],:] + N.matrixmultiply(cho[S[i]*d:S[i]*d+d,:], X[i,:])
        
        # Faster:
        cho = N.zeros((K, va.shape[1], va.shape[1]), float)
        for k in range(K):
            # Using cholesky is more stable than sqrtm; sqrtm is not
            # available in numpy anyway, only in scipy...
            cho[k]  = lin.cholesky(va[k*d:k*d+d,:])

        for s in range(K):
            tmpind      = N.where(S == s)[0]
            X[tmpind]   = N.matrixmultiply(X[tmpind], cho[s].transpose()) + mu[s]
    else:
        raise GmmParamError('cov matrix mode not recognized')

    return X

def gen_gmm_param(d, K, varmode = 'diag', spread = 1):
    """Generate valid parameters for a gaussian mixture model.
    d is the dimension, K the number of components, and varmode
    the mode for cov matrices.

    Returns:    
        - w
        - mu
        - va
    """
    w   = abs(N.randn(K))
    w   = w / sum(w)

    mu  = spread * N.randn(K, d)
    if varmode == 'diag':
        va  = abs(N.randn(K, d))
    elif varmode == 'full':
        va  = N.randn(K * d, d)
        for k in range(K):
            va[k*d:k*d+d]   = N.matrixmultiply( va[k*d:k*d+d], 
                va[k*d:k*d+d].transpose())
    else:
        raise GmmParamError('cov matrix mode not recognized')

    return w, mu, va

def check_gmm_param(w, mu, va):
    """Check that w, mu and va are valid parameters for
    a mixture of gaussian: w should sum to 1, there should
    be the same number of component in each param, the variances
    should be positive definite, etc... 
    
    Params:
        w   = vector or list of weigths of the mixture (K elements)
        mu  = matrix: K * d
        va  = list of variances (vector K * d or square matrices Kd * d)

    returns:
        K   = number of components
        d   = dimension
        mode    = 'diag' if diagonal covariance, 'full' of full matrices
    """
        
    # Check that w is valid
    if N.fabs(N.sum(w)  - 1) > MAX_DEV:
        raise GmmParamError('weight does not sum to 1')
    
    if not len(w.shape) == 1:
        raise GmmParamError('weight is not a vector')

    # Check that mean and va have the same number of components
    K           = len(w)

    if N.ndim(mu) < 2:
        msg = "mu should be a K,d matrix, and a row vector if only 1 comp"
        raise GmmParamError(msg)
    if N.ndim(va) < 2:
        msg = """va should be a K,d / K *d, d matrix, and a row vector if
        only 1 diag comp"""
        raise GmmParamError(msg)

    (Km, d)     = mu.shape
    (Ka, da)    = va.shape

    if not K == Km:
        msg = "not same number of component in mean and weights"
        raise GmmParamError(msg)

    if not d == da:
        msg = "not same number of dimensions in mean and variances"
        raise GmmParamError(msg)

    if Km == Ka:
        mode = 'diag'
    else:
        mode = 'full'
        if not Ka == Km*d:
            msg = "not same number of dimensions in mean and variances"
            raise GmmParamError(msg)
        
    return K, d, mode
        
# For EM on GMM
def multiple_gauss_den(data, mu, va):
    """Helper function to generate several Gaussian
    pdf (different parameters) from the same data"""
    mu  = N.atleast_2d(mu)
    va  = N.atleast_2d(va)

    K   = mu.shape[0]
    n   = data.shape[0]
    d   = data.shape[1]
    
    y   = N.zeros((n, K), float)
    if mu.size == va.size:
        for i in range(K):
            y[:, i] = densities.gauss_den(data, mu[i, :], va[i, :])
    else:
        for i in range(K):
            y[:, i] = densities.gauss_den(data, mu[i, :], 
                        va[d*i:d*i+d, :])

    return y

def gmm_init_kmean(data, k, mode, init = [], niter = 10):
    """gmm_init_kmean(data, k, mode, init = [], niter = 10)
    
    Init EM for GMM with kmean from data, for k components. 
    
    Args:
        - data:     Each row of data is one frame of dimension d. 
        - k:        Number of components to look for
        - mode:     Diagonal or Full covariance matrices
        - init:     The initial centroids. The value given for k is
        ignored, and the number of row in initc is used instead. 
        If initc is not given, then the centroids are initialized 
        with the k first values of data.
        - niter:    Number of iterations in kmean.
    
    Returns:
        (w, mu, va), initial parameters for a GMM.

    Method:
        Each weight is equiprobable, each mean is one centroid returned by kmean, and
    covariances for component i is initialized with covariance of 
    data corresponding with label i. Other strategies are possible, this one
    is an easy one"""
    if len(init) == 0:
        init   = data[0:k, :]
    else:
        k       = initc.shape[0]

    if data.ndim == 1:
        d   = 1
    else:
        d   = N.shape(data)[1]

    (code, label)   = kmean(data, init, niter)

    w   = N.ones(k, float) / k
    mu  = code.copy()
    if mode == 'diag':
        va  = N.zeros((k, d), float)
        for i in range(k):
            for j in range(d):
                va[i,j] = N.cov(data[N.where(label==i), j], rowvar = 0)
    elif mode == 'full':
        va  = N.zeros((k*d, d), float)
        for i in range(k):
            va[i*d:i*d+d,:] = N.cov(data[N.where(label==i)], rowvar = 0)
    else:
        raise GmmParamError("mode " + str(mode) + " not recognized")

    return w, mu, va

def gmm_posterior(data, w, mu, va):
    """ Computes the latent variable distribution (a 
    posteriori probability) knowing the explicit data 
    for the Gaussian model (w, mu, var): gamma(t, i) = 
        P[state = i | observation = data(t); w, mu, va]

    This is basically the E step of EM for GMM.
   
    the second returned value is the non normalized version 
    of gamma, and may be needed for some computation, 
    like eg likelihood"""
    n   = data.shape[0]
    K   = len(w)

    # compute the gaussian pdf
    tgd	= multiple_gauss_den(data, mu, va)
    # multiply by the weight
    tgd	*= w
    # Normalize to get a pdf
    gd	= tgd  / N.sum(tgd, axis=1)[:, N.NewAxis]

    return gd, tgd

# This function is just calling gmm_update and gmm_posterior, with
# initialization. This is ugly, and we should have a class to model a GMM
# instead of all this code to try to guess correct values and parameters...
def gmm_em(data, niter = 10, k = 2, mode = 'diag', w = [], mu = [], va = []):
    """
    gmm_em(data, niter = 10, k = 2, mode = 'diag', w = [], mu = [], va = []):

    Compute the parameters of a Gaussian Mixture Model using EM algorithm, 
    with initial values w, mu and va (overwritten by the function).

    Args:
        - data:     contains the observed features, one row is one frame, ie one 
        observation of dimension d
        - niter:    number of iterations
        - mode:     'diag' or 'full', depending on the wanted model for cov
        matrices.
        - K:        number of components
        - w, mu, va initial parameters for the GMM. All or none must be given.
        If no initial values are given, initialized by gmm_init_kmean; if they
        are given, mode and k are ignored, and guessed from the given parameters
        instead.

    Returns:
        w, mu, va, like as found by EM, where like is the likelihood for each 
        iteration.
    """
    if len(w) == 0:
        w, mu, va   = gmm_init_kmean(data, k, mode, 5)
    else:
        k, d, mode  = check_gmm_parm(w, mu, va)
    
    like    = N.zeros(niter, float)
    for i in range(niter):
        g, tgd  = gmm_posterior(data, w, mu, va)
        like[i] = N.sum(N.log(N.sum(tgd, 1)))
        w, mu, va   = gmm_update(X, g, d, k, mode)

    return w, mu, va, like
    
def gmm_update(data, gamma, d, K, varmode):
    """Computes update of the Gaussian Mixture Model (M step)
    from the a posteriori pdf, computed by gmm_posterior
    (E step).
    """
    n       = data.shape[0]
    invn    = 1.0/n
    mGamma  = N.sum(gamma)

    if varmode == 'diag':
        mu  = N.zeros((K, d), float)
        va  = N.zeros((K, d), float)
        for k in range(K):
            x       = N.sum(N.outerproduct(gamma[:, k], 
                        N.ones((1, d))) * data)
            xx      = N.sum(N.outerproduct(gamma[:, k], 
                        N.ones((1, d))) * (data ** 2))

            mu[k,:] = x / mGamma[k]
            va[k,:] = xx  / mGamma[k] - mu[k,:] ** 2
        w   = invn * mGamma

    elif varmode == 'full':
        mu  = N.zeros((K, d), float)
        va  = N.zeros((K*d, d), float)

        # This is really slow because of the loop on each frame
        # ways to improve that: either use 3d matrices and see 
        # the memory usage increasing in n*d*d,
        # using C (or pyrex, or ...) or using the trick commented
        # which loops on the dimensions instead. If the number of dimensions
        # is high, this won't help, though...
        for k in range(K):
            x   = N.sum(N.outerproduct(gamma[:, k], 
                        N.ones((1, d), float)) * data)
            xx  = N.zeros((d, d), float)
            # for i in range(n):
            #     xx  += gamma[i, k] * N.outerproduct(data[i,:], 
            #                 data[i,:])

            # This should be much faster
            for i in range(d):
                for j in range(d):
                    xx[i,j] = N.sum(data[:,i] * data[:,j] * gamma[:,k])

            mu[k,:] = x / mGamma[k]
            va[k*d:k*d+d,:] = xx  / mGamma[k] - \
                                N.outerproduct(mu[k,:], mu[k,:])
        w   = invn * mGamma
    else:
        raise GmmParamError("varmode not recognized")

    return w, mu, va

# Misc functions
def gmm_ellipses(mu, va, c = [0, 1], npoints = 100):
    """Returns a list of ellipses describing the Gmm
    defined by mu and va. c is the dimension we are projecting
    the variances on a 2d space.
    
    Returns:
        -Xe:    a list of x coordinates for the ellipses
        -Ye:    a list of y coordinates for the ellipses

    Example:
        Suppose we have w, mu and va as parameters for a mixture, then:
        
        X       = gen_gmm(w, mu, va, 1000)
        Xe, Ye  = gmm_ellipses(mu, va)
        pylab.plot(X[:,0], X[:, 1], '.')
        for k in len(w):
            pylab.plot(Xe[k], Ye[k], 'r')
            
        Will plot samples X draw from the mixture model, and
        plot the ellipses of equi-probability from the mean with
        fixed level of confidence 0.39. 
        
    TODO: be able to modify the confidence interval to arbitrary
    value (to do in gauss_ell)"""
    K   = mu.shape[0]
    w   = N.ones(K, float) / K
    
    K, d, mode  = check_gmm_param(w, mu, va)

    # TODO: adjustable level (to do in gauss_ell). 
    # For now, a level of 0.39 means that we draw
    # ellipses for 1 standard deviation. 
    Xe  = []
    Ye  = []   
    if mode == 'diag':
        for i in range(K):
            xe, ye  = densities.gauss_ell(mu[i,:], va[i,:], dim = c, 
                    npoints = npoints)
            Xe.append(xe)
            Ye.append(ye)
    elif mode == 'full':
        for i in range(K):
            xe, ye  = densities.gauss_ell(mu[i,:], va[i*d:i*d+d,:], dim = c, 
                    npoints = npoints)
            Xe.append(xe)
            Ye.append(ye)

    return Xe, Ye

import c_gmm
def kmean(data, init, iter = 10):
    """Simple kmean implementation for EM
    
    returns a tuple (code, label), where code are the final
    centroids, and label are the class label indec for each
    frame (ie row) of data"""

    data    = N.atleast_2d(data)
    init    = N.atleast_2d(init)

    (n, d)  = data.shape
    (k, d1) = init.shape

    if not d == d1:
        msg = "data and init centers do not have same dimensions..."
        raise GmmParamError(msg)
    
    code    = N.asarray(init.copy(), float)
    for i in range(iter):
        # Compute the nearest neighbour for each obs
        # using the current code book
        label   = c_gmm._vq(data, code)
        # Update the code by computing centroids using the new code book
        for j in range(k):
            code[j,:] = N.mean(data[N.where(label==j)]) 

    return code, label

def _vq(data, code):
    """ Please do not use directly. Use kmean instead"""
    # No attempt to be efficient has been made...
    (n, d)  = data.shape
    (k, d)  = code.shape

    label   = N.zeros(n, int)
    for i in range(n):
        d           = N.sum((data[i, :] - code) ** 2, 1)
        label[i]    = N.argmin(d)

    return label
    
# Test functions usable for now
def test_kmean():
    X   = N.array([[3.0, 3], [4, 3], [4, 2],
            [9, 2], [5, 1], [6, 2], [9, 4], 
            [5, 2], [5, 4], [7, 4], [6, 5]])

    initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 

    codet1  = N.array([[3.0000, 3.0000],
            [6.2000, 4.0000], 
            [5.8000, 1.8000]])
            
    codet2  = N.array([[11.0/3, 8.0/3], 
            [6.7500, 4.2500],
            [6.2500, 1.7500]])

    code    = initc.copy()

    code1   = kmean(X, code, 1)[0]
    code2   = kmean(X, code, 2)[0]

    import numpy.testing as testing
    try:
        testing.assert_array_almost_equal(code1, codet1)
        print "kmean test 1 succeded"
    except AssertionError:
        print "kmean test 1 failed"

    try:
        testing.assert_array_almost_equal(code2, codet2)
        print "kmean test 2 succeded"
    except AssertionError:
        print "kmean test 2 failed"

# Test functions (unusable by other people for now....)
def test_gmm_den():
    """"""
    import tables
    
    filename    = 'gmm_test.h5'
    
    k   = 3
    d   = 2
    mode        = 'full'
    w, mu, va   = gen_gmm_param(d, k, mode, 2)
    X           = gen_gmm(w, mu, va, 1e3)

    h5file      = tables.openFile(filename, "w")

    h5file.createArray(h5file.root, 'mu', mu)
    h5file.createArray(h5file.root, 'va', va)
    h5file.createArray(h5file.root, 'X', X)

    Y1 =  densities.gauss_den(X, mu[0,:], va[0:2,:])
    Y2 =  densities.gauss_den(X, mu[1,:], va[2:4,:])
    Y3 =  densities.gauss_den(X, mu[2,:], va[4:6,:])

    Y  =  multiple_gauss_den(X, mu, va)

    h5file.createArray(h5file.root, 'Y', Y)
    h5file.createArray(h5file.root, 'Y2', Y2)
    h5file.createArray(h5file.root, 'Y1', Y1)
    h5file.createArray(h5file.root, 'Y3', Y3)

    h5file.close()

def test_gmm_em():
    """"""
    import numpy as N
    import tables
    import pylab as P
    filename    = 'data.h5'
    
    k               = 3
    d               = 2
    mode            = 'full'
    wr, mur, var    = gen_gmm_param(d, k, mode, 2)
    X               = gen_gmm(wr, mur, var, 1e3)

    h5file      = tables.openFile(filename, "w")

    h5file.createArray(h5file.root, 'wr', wr)
    h5file.createArray(h5file.root, 'mur', mur)
    h5file.createArray(h5file.root, 'var', var)
    h5file.createArray(h5file.root, 'X', X)

    P.plot(X[:, 0], X[:, 1], '.')
    Xe, Ye  = gmm_ellipses(mur, var, npoints = 100)
    for i in range(k):
        P.plot(Xe[i], Ye[i], 'g')

    # init EM: simply use kmean
    initc   = X[0:k, :]
    (code, label)   = kmean(X, initc, 10)

    h5file.createArray(h5file.root, 'initc', initc)

    w0      = N.ones(k, float) / k
    mu0     = code.copy()
    if mode == 'diag':
        va0     = N.zeros((k, d), float)
        for i in range(k):
            for j in range(d):
                va0[i,j] = N.cov(X[N.where(label==i), j], rowvar = 0)
    elif mode == 'full':
        va0     = N.zeros((k*d, d), float)
        for i in range(k):
            va0[i*d:i*d+d,:] = N.cov(X[N.where(label==i)], rowvar = 0)

    h5file.createArray(h5file.root, 'w0', w0)
    h5file.createArray(h5file.root, 'mu0', mu0)
    h5file.createArray(h5file.root, 'va0', va0)

    
    # # Use random values
    # w0  = N.ones(k, float) / k
    # mu0 = N.randn(k, d)
    # va0 = N.fabs(N.randn(k, d))

    w   = w0.copy()
    mu  = mu0.copy()
    va  = va0.copy()

    X0e, Y0e  = gmm_ellipses(mu, va, npoints = 100)
    for i in range(k):
        P.plot(X0e[i], Y0e[i], 'k')

    niter   = 10
    like    = N.zeros(niter, float)

    for i in range(niter):
        g, tgd  = gmm_posterior(X, w, mu, va)
        #h5file.createArray(h5file.root, 'g', g)
        #h5file.createArray(h5file.root, 'tgd', tgd)
        like[i] = N.sum(N.log(N.sum(tgd, 1)))
        print like[i]
        w, mu, va    = gmm_update(X, g, d, k, mode)

    h5file.createArray(h5file.root, 'w', w)
    h5file.createArray(h5file.root, 'mu', mu)
    h5file.createArray(h5file.root, 'va', va)

    X0e, Y0e  = gmm_ellipses(mu, va, npoints = 100)
    for i in range(k):
        P.plot(X0e[i], Y0e[i], 'r')
    P.figure()
    P.plot(like)
    h5file.close()
    P.show()

def _bench1():
    #===========================================
    # Diag GMM with 5 components of dimension 20
    #===========================================
    k       = 5
    d       = 20
    mode    = 'diag'

    print "Generating the mixture"
    # Generate a model with k components, d dimensions
    wr, mur, var    = gen_gmm_param(d, k, mode, 3)
    X               = gen_gmm(wr, mur, var, 1e4)

    print "Init the mixture"
    # Init the mixture with kmean
    w0, mu0, va0    = gmm_init_kmean(X, k, mode, niter = 5)
    
    # # Use random values instead of kmean
    # w0  = N.ones(k, float) / k
    # mu0 = N.randn(k, d)
    # va0 = N.fabs(N.randn(k, d))

    print "EM computing..."
    # Copy the initial values because we want to draw them later...
    w   = w0.copy()
    mu  = mu0.copy()
    va  = va0.copy()

    # The actual EM, with likelihood computation
    niter   = 10
    like    = N.zeros(niter, float)

    for i in range(niter):
        print "post"
        g, tgd  = gmm_posterior(X, w, mu, va)
        like[i] = N.sum(N.log(N.sum(tgd, 1)))
        print "update"
        w, mu, va   = gmm_update(X, g, d, k, mode)


def benchmark():
    import profile
    profile.run('_bench1()', 'bench1prof')

if __name__ == "__main__":
    #=============================
    # Simple GMM with 3 components
    #=============================
    import pylab as P
    k       = 4
    d       = 2
    mode    = 'diag'

    print "Generating the mixture"
    # Generate a model with k components, d dimensions
    wr, mur, var    = gen_gmm_param(d, k, mode, 3)
    X               = gen_gmm(wr, mur, var, 1e3)

    print "Init the mixture"
    # Init the mixture with kmean
    w0, mu0, va0    = gmm_init_kmean(X, k, mode, niter = 5)
    
    # # Use random values instead of kmean
    # w0  = N.ones(k, float) / k
    # mu0 = N.randn(k, d)
    # va0 = N.fabs(N.randn(k, d))

    # Copy the initial values because we want to draw them later...
    w   = w0.copy()
    mu  = mu0.copy()
    va  = va0.copy()

    # The actual EM, with likelihood computation
    niter   = 10
    like    = N.zeros(niter, float)

    print "computing..."
    for i in range(niter):
        print "post"
        g, tgd  = gmm_posterior(X, w, mu, va)
        like[i] = N.sum(N.log(N.sum(tgd, 1)))
        print "update"
        w, mu, va   = gmm_update(X, g, d, k, mode)

    print "drawing..."
    # Draw what is happening
    P.subplot(2, 1, 1)
    P.plot(X[:, 0], X[:, 1], '.', label = '_nolegend_')

    # Real confidence ellipses
    Xre, Yre  = gmm_ellipses(mur, var)
    P.plot(Xre[0], Yre[0], 'g', label = 'true confidence ellipsoides')
    for i in range(1,k):
        P.plot(Xre[i], Yre[i], 'g', label = '_nolegend_')

    # Initial confidence ellipses as found by kmean
    X0e, Y0e  = gmm_ellipses(mu0, va0)
    P.plot(X0e[0], Y0e[0], 'k', label = 'initial confidence ellipsoides')
    for i in range(1,k):
        P.plot(X0e[i], Y0e[i], 'k', label = '_nolegend_')

    # Values found by EM
    Xe, Ye  = gmm_ellipses(mu, va)
    P.plot(Xe[0], Ye[0], 'r', label = 'confidence ellipsoides found by EM')
    for i in range(1,k):
        P.plot(Xe[i], Ye[i], 'r', label = '_nolegend_')
    P.legend(loc = 0)
    P.subplot(2, 1, 2)
    P.plot(like)
    P.title('log likelihood')

    # # Export the figure
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
