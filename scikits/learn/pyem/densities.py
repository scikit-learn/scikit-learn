#! /usr/bin/python
#
# Copyrighted David Cournapeau
# Last Change: Thu Jul 12 04:00 PM 2007 J
"""This module implements various basic functions related to multivariate
gaussian, such as pdf estimation, confidence interval/ellipsoids, etc..."""

__docformat__ = 'restructuredtext'

import numpy as N
import numpy.linalg as lin
#from numpy.random import randn
from scipy.stats import chi2
import misc

# Error classes
class DenError(Exception):
    """Base class for exceptions in this module.
    
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error"""
    def __init__(self, message):
        self.message    = message
        Exception.__init__(self)
    
    def __str__(self):
        return self.message

# The following function do all the fancy stuff to check that parameters
# are Ok, and call the right implementation if args are OK.
def gauss_den(x, mu, va, log = False):
    """Compute multivariate Gaussian density at points x for 
    mean mu and variance va.
    
    :Parameters:
        x : ndarray
            points where to estimate the pdf.  each row of the array is one
            point of d dimension
        mu : ndarray
            mean of the pdf. Should have same dimension d than points in x.
        va : ndarray
            variance of the pdf. If va has d elements, va is interpreted as the
            diagonal elements of the actual covariance matrix. Otherwise,
            should be a dxd matrix (and positive definite).
        log : boolean
            if True, returns the log-pdf instead of the pdf.

    :Returns:
        pdf : ndarray
            Returns a rank 1 array of the pdf at points x.

    Note
    ----
        Vector are row vectors, except va which can be a matrix
        (row vector variance for diagonal variance)."""
    
    lmu  = N.atleast_2d(mu)
    lva  = N.atleast_2d(va)
    lx   = N.atleast_2d(x)
    
    #=======================#
    # Checking parameters   #
    #=======================#
    if len(N.shape(lmu)) != 2:
        raise DenError("mu is not rank 2")
        
    if len(N.shape(lva)) != 2:
        raise DenError("va is not rank 2")
        
    if len(N.shape(lx)) != 2:
        raise DenError("x is not rank 2")
        
    d = N.shape(lx)[1]
    (dm0, dm1) = N.shape(lmu)
    (dv0, dv1) = N.shape(lva)
    
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
        return _scalar_gauss_den(lx[:, 0], lmu[0, 0], lva[0, 0], log)
    elif dv0 == 1:
        # Diagonal matrix case
        return _diag_gauss_den(lx, lmu, lva, log)
    elif dv1 == dv0:
        # full case
        return  _full_gauss_den(lx, lmu, lva, log)
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
    inva    *= -0.5
    y       = ((x-mu) ** 2) * inva
    if not log:
        y   = fac * N.exp(y)
    else:
        y   += N.log(fac)

    return y
    
def _diag_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    Call gauss_den instead"""
    # Diagonal matrix case
    d   = mu.size
    #n   = x.shape[0]
    if not log:
        inva = 1/va[0]
        fac = (2*N.pi) ** (-d/2.0) * N.prod(N.sqrt(inva))
        inva *= -0.5
        x = x - mu
        x **= 2
        y = fac * N.exp(N.dot(x, inva))
    else:
        # XXX optimize log case as non log case above
        y = _scalar_gauss_den(x[:, 0], mu[0, 0], va[0, 0], log)
        for i in range(1, d):
            y +=  _scalar_gauss_den(x[:, i], mu[0, i], va[0, i], log)
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

    # we are using a trick with sum to "emulate" 
    # the matrix multiplication inva * x without any explicit loop
    #y   = -0.5 * N.sum(N.dot((x-mu), inva) * (x-mu), 1)
    y   = -0.5 * N.dot(N.dot((x-mu), inva) * (x-mu), 
                       N.ones((mu.size, 1), x.dtype))[:, 0]

    if not log:
        y   = fac * N.exp(y)
    else:
        y   = y + N.log(fac)
 
    return y

# To get coordinatea of a confidence ellipse from multi-variate gaussian pdf
def gauss_ell(mu, va, dim = misc.DEF_VIS_DIM, npoints = misc.DEF_ELL_NP, \
        level = misc.DEF_LEVEL):
    """Given a mean and covariance for multi-variate
    gaussian, returns the coordinates of the confidense ellipsoid.
    
    Compute npoints coordinates for the ellipse of confidence of given level
    (all points will be inside the ellipsoides with a probability equal to
    level).
    
    :Parameters:
        mu : ndarray
            mean of the pdf
        va : ndarray
            variance of the pdf
        dim : sequence
            sequences of two integers which represent the dimensions where to
            project the ellipsoid.
        npoints: int
            number of points to generate for the ellipse.
        level : float
            level of confidence (between 0 and 1).

    :Returns:
        Returns the coordinate x and y of the ellipse."""
    if level >= 1 or level <= 0:
        raise ValueError("level should be a scale strictly between 0 and 1.""")
    
    mu = N.atleast_1d(mu)
    va = N.atleast_1d(va)
    d = N.shape(mu)[0]
    c = N.array(dim)

    if N.any(c < 0) or N.any(c >= d):
        raise ValueError("dim elements should be >= 0 and < %d (dimension"\
                " of the variance)" % d)
    if N.size(mu) == N.size(va):
        mode    = 'diag'
    else:
        if N.ndim(va) == 2:
            if N.shape(va)[0] == N.shape(va)[1]:
                mode    = 'full'
            else:
                raise DenError("variance not square")
        else:
            raise DenError("mean and variance are not dim conformant")

    # When X is a sample from multivariante N(mu, sigma), (X-mu)Sigma^-1(X-mu)
    # follows a Chi2(d) law. Here, we only take 2 dimension, so Chi2 with 2
    # degree of freedom (See Wasserman. This is easy to see with characteristic
    # functions)
    chi22d  = chi2(2)
    mahal   = N.sqrt(chi22d.ppf(level))
    
    # Generates a circle of npoints
    theta   = N.linspace(0, 2 * N.pi, npoints)
    circle  = mahal * N.array([N.cos(theta), N.sin(theta)])

    # Get the dimension which we are interested in:
    mu  = mu[c]
    if mode == 'diag':
        va      = va[c]
        elps    = N.outer(mu, N.ones(npoints))
        elps    += N.dot(N.diag(N.sqrt(va)), circle)
    elif mode == 'full':
        va  = va[c, :][:, c]
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
        raise ValueError("var mode not recognized")

    return elps[0, :], elps[1, :]

def logsumexp(x):
    """Compute log(sum(exp(x), 1)) while avoiding underflow.
    
    :Parameters:
        x : ndarray
            data in log domain to sum"""
    axis = 1
    mc = N.max(x, axis)
    return mc + N.log(N.sum(N.exp(x-mc[:, N.newaxis]), axis))

def multiple_gauss_den(data, mu, va, log = False):
    """Helper function to generate several Gaussian
    pdf (different parameters) at the same points

    :Parameters:
        data : ndarray
            points where to estimate the pdfs (n,d).
        mu : ndarray
            mean of the pdf, of shape (k,d). One row of dimension d per
            different component, the number of rows k being the number of
            component
        va : ndarray
            variance of the pdf. One row per different component for diagonal
            covariance (k, d), or d rows per component for full matrix pdf
            (k*d,d).
        log : boolean
            if True, returns the log-pdf instead of the pdf.

    :Returns:
        Returns a (n, k) array, each column i being the pdf of the ith mean and
        ith variance."""
    mu = N.atleast_2d(mu)
    va = N.atleast_2d(va)

    k = N.shape(mu)[0]
    n = N.shape(data)[0]
    d = N.shape(mu)[1]
    
    y = N.zeros((k, n))
    if N.size(mu) == N.size(va):
        for i in range(k):
            y[i] = gauss_den(data, mu[i, :], va[i, :], log)
        return y.T
    else:
        for i in range(k):
            y[i] = gauss_den(data, mu[i, :], va[d*i:d*i+d, :], log)
        return y.T

if __name__ == "__main__":
    pass
