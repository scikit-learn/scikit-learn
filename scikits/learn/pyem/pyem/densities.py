#! /usr/bin/python
#
# Copyrighted David Cournapeau
# Last Change: Wed Aug 30 12:00 PM 2006 J

import numpy as N
import numpy.linalg as lin
from numpy.random import randn
from scipy.stats import chi2

# Error classes
class DenError(Exception):
    """Base class for exceptions in this module.
    
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error"""
    def __init__(self, message):
        self.message    = message
    
    def __str__(self):
        return self.message

# The following function do all the fancy stuff to check that parameters
# are Ok, and call the right implementation if args are OK.
def gauss_den(x, mu, va, log = False):
    """ Compute multivariate Gaussian density at points x for 
    mean mu and variance va.
    
    Vector are row vectors, except va which can be a matrix
    (row vector variance for diagonal variance)
    
    If log is True, than the log density is returned 
    (useful for underflow ?)"""
    mu  = N.atleast_2d(mu)
    va  = N.atleast_2d(va)
    x   = N.atleast_2d(x)
    
    #=======================#
    # Checking parameters   #
    #=======================#
    if len(N.shape(mu)) != 2:
        raise DenError("mu is not rank 2")
        
    if len(N.shape(va)) != 2:
        raise DenError("va is not rank 2")
        
    if len(N.shape(x)) != 2:
        raise DenError("x is not rank 2")
        
    (n, d)      = x.shape
    (dm0, dm1)  = mu.shape
    (dv0, dv1)  = va.shape
    
    # Check x and mu same dimension
    if dm0 != 1:
        msg = "mean must be a row vector!"
        raise DenError(msg)
    if dm1 != d:
        msg = "x and mu not same dim"
        raise DenError(msg)
    # Check va and mu same size
    if dv1 != d:
        msg = "mu and va not same dim"
        raise DenError(msg)
    if dv0 != 1 and dv0 != d:
        msg = "va not square"
        raise DenError(msg)

    #===============#
    # Computation   #
    #===============#
    if d == 1:
        # scalar case
        return _scalar_gauss_den(x[:, 0], mu[0, 0], va[0, 0], log)
    elif dv0 == 1:
        # Diagonal matrix case
        return _diag_gauss_den(x, mu, va, log)
    elif dv1 == dv0:
        # full case
        return  _full_gauss_den(x, mu, va, log)
    else:
        raise DenError("variance mode not recognized, this is a bug")

# Those 3 functions do almost all the actual computation
def _scalar_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    Call gauss_den instead"""
    d       = mu.size
    inva    = 1/va
    fac     = (2*N.pi) ** (-d/2.0) * N.sqrt(inva)
    y       = ((x-mu) ** 2) * -0.5 * inva
    if not log:
        y   = fac * N.exp(y)
    else:
        y   = y + log(fac)

    return y
    
#from ctypes import cdll, c_uint, c_int, c_double, POINTER
#_gden   = cdll.LoadLibrary('src/libgden.so')
#_gden.gden_diag.restype     = c_int
#_gden.gden_diag.argtypes    = [POINTER(c_double), c_uint, c_uint,
#        POINTER(c_double), POINTER(c_double), POINTER(c_double)]

def _diag_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    Call gauss_den instead"""
    # Diagonal matrix case
    d   = mu.size
    if not log:
        inva    = 1/va[0,0]
        fac     = (2*N.pi) ** (-d/2.0) * N.sqrt(inva)
        y       =  (x[:,0] - mu[0,0]) ** 2 * inva * -0.5
        for i in range(1, d):
            inva    = 1/va[0,i]
            #fac     *= (2*N.pi) ** (-d/2.0) * N.sqrt(inva)
            fac     *= N.sqrt(inva)
            y       += (x[:,i] - mu[0,i]) ** 2 * inva * -0.5
        y   = fac * N.exp(y)
    #d   = mu.size
    #n   = x.shape[0]
    #if not log:
    #    y   = N.zeros(n)
    #    inva= 1/va
    #    _gden.gden_diag(x.ctypes.data_as(POINTER(c_double)),
    #        n, d,
    #        mu.ctypes.data_as(POINTER(c_double)),
    #        (inva).ctypes.data_as(POINTER(c_double)),
    #        y.ctypes.data_as(POINTER(c_double)))
    else:
        y   = _scalar_gauss_den(x[:,0], mu[0,0], va[0,0], log)
        for i in range(1, d):
            y    +=  _scalar_gauss_den(x[:,i], mu[0,i], va[0,i], log)
    return y

def _full_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in full matrix case. 
    
    It assumes all args are conformant, so it should 
    not be used directly Call gauss_den instead
    
    Does not check if va is definite positive (on inversible 
    for that matter), so the inverse computation and/or determinant
    would throw an exception."""
    d       = mu.size
    inva    = lin.inv(va)
    fac     = 1 / N.sqrt( (2*N.pi) ** d * N.fabs(lin.det(va)))

    # # Slow version
    # n       = N.size(x, 0)
    # y       = N.zeros(n)
    # for i in range(n):
    #     y[i] = N.dot(x[i,:],
    #              N.dot(inva, N.transpose(x[i,:])))
    # y *= -0.5

    # we are using a trick with sum to "emulate" 
    # the matrix multiplication inva * x without any explicit loop
    y   = N.dot((x-mu), inva)
    y   = -0.5 * N.sum(y * (x-mu), 1)

    if not log:
        y   = fac * N.exp(y)
    else:
        y   = y + N.log(fac)
 
    return y

# To plot a confidence ellipse from multi-variate gaussian pdf
def gauss_ell(mu, va, dim = [0, 1], npoints = 100, level = 0.39):
    """ Given a mean and covariance for multi-variate
    gaussian, returns npoints points for the ellipse
    of confidence given by level (all points will be inside
    the ellipsoides with a probability equal to level)
    
    Returns the coordinate x and y of the ellipse"""
    
    mu      = N.atleast_1d(mu)
    va      = N.atleast_1d(va)
    c       = N.array(dim)

    if mu.size == va.size:
        mode    = 'diag'
    else:
        if va.ndim == 2:
            if va.shape[0] == va.shape[1]:
                mode    = 'full'
            else:
                raise DenError("variance not square")
        else:
            raise DenError("mean and variance are not dim conformant")

    chi22d  = chi2(2)
    mahal   = N.sqrt(chi22d.ppf(level))
    
    # Generates a circle of npoints
    theta   = N.linspace(0, 2 * N.pi, npoints)
    circle  = mahal * N.array([N.cos(theta), N.sin(theta)])

    # Get the dimension which we are interested in:
    mu  = mu[dim]
    if mode == 'diag':
        va      = va[dim]
        elps    = N.outer(mu, N.ones(npoints))
        elps    += N.dot(N.diag(N.sqrt(va)), circle)
    elif mode == 'full':
        va  = va[c,:][:,c]
        # Method: compute the cholesky decomp of each cov matrix, that is
        # compute cova such as va = cova * cova' 
        # WARN: scipy is different than matlab here, as scipy computes a lower
        # triangular cholesky decomp: 
        #   - va = cova * cova' (scipy)
        #   - va = cova' * cova (matlab)
        # So take care when comparing results with matlab !
        cova    = lin.cholesky(va)
        elps    = N.outer(mu, N.ones(npoints))
        elps    += N.dot(cova, circle)
    else:
        raise DenParam("var mode not recognized")

    return elps[0, :], elps[1, :]

if __name__ == "__main__":
    import pylab

    #=========================================
    # Test plotting a simple diag 2d variance:
    #=========================================
    va  = N.array([5, 3])
    mu  = N.array([2, 3])

    # Generate a multivariate gaussian of mean mu and covariance va
    X       = randn(1e3, 2)
    Yc      = N.dot(N.diag(N.sqrt(va)), X.transpose())
    Yc      = Yc.transpose() + mu

    # Plotting
    Xe, Ye  = gauss_ell(mu, va, npoints = 100)
    pylab.figure()
    pylab.plot(Yc[:, 0], Yc[:, 1], '.')
    pylab.plot(Xe, Ye, 'r')

    #=========================================
    # Test plotting a simple full 2d variance:
    #=========================================
    va  = N.array([[0.2, 0.1],[0.1, 0.5]])
    mu  = N.array([0, 3])

    # Generate a multivariate gaussian of mean mu and covariance va
    X       = randn(1e3, 2)
    Yc      = N.dot(lin.cholesky(va), X.transpose())
    Yc      = Yc.transpose() + mu

    # Plotting
    Xe, Ye  = gauss_ell(mu, va, npoints = 100, level=0.95)
    pylab.figure()
    pylab.plot(Yc[:, 0], Yc[:, 1], '.')
    pylab.plot(Xe, Ye, 'r')
    pylab.show()
