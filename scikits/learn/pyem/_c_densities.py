#! /usr/bin/python
#
# Copyrighted David Cournapeau
# Last Change: Thu Nov 09 05:00 PM 2006 J

# This module uses a C implementation through ctypes, for diagonal cases
# TODO:
#   - portable way to find/open the shared library
#   - full cov matrice

import numpy as N
import numpy.linalg as lin
from numpy.random import randn
from scipy.stats import chi2
import densities as D

import ctypes
from ctypes import cdll, c_uint, c_int, c_double, POINTER
from numpy.ctypeslib import ndpointer, load_library

ctypes_major    = int(ctypes.__version__.split('.')[0])
if ctypes_major < 1:
    msg =  "version of ctypes is %s, expected at least %s" \
            % (ctypes.__version__, '1.0.0')
    raise ImportError(msg)

# Requirements for diag gden
_gden   = load_library('c_gden.so', __file__)
arg1    = ndpointer(dtype=N.float64)
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype=N.float64)
arg5    = ndpointer(dtype=N.float64)
arg6    = ndpointer(dtype=N.float64)
_gden.gden_diag.argtypes    = [arg1, arg2, arg3, arg4, arg5, arg6]
_gden.gden_diag.restype     = c_int

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
    
    ** Expect centered data (ie with mean removed) **

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
    
def _diag_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    ** Expect centered data (ie with mean removed) **

    Call gauss_den instead"""
    # Diagonal matrix case
    d   = mu.size
    n   = x.shape[0]
    if not log:
        y       = N.zeros(n)
        vat     = va.copy()
        # _gden.gden_diag(N.require(x, requirements = 'C'), n, d, 
        #         N.require(mu, requirements = 'C'),
        #         N.require(inva, requirements = 'C'),
        #         N.require(y, requirements = 'C'))
        x       = N.require(x, requirements = 'C')
        mu      = N.require(mu, requirements = 'C')
        vat     = N.require(vat, requirements = 'C')
        y       = N.require(y, requirements = 'C')
        _gden.gden_diag(x, n, d, mu, vat, y)
        return y
        # _gden.gden_diag.restype     = c_int
        # _gden.gden_diag.argtypes    = [POINTER(c_double), c_uint, c_uint,
        #         POINTER(c_double), POINTER(c_double), POINTER(c_double)]

        # y   = N.zeros(n)
        # inva= 1/va
        # _gden.gden_diag(x.ctypes.data_as(POINTER(c_double)),
        #     n, d,
        #     mu.ctypes.data_as(POINTER(c_double)),
        #     inva.ctypes.data_as(POINTER(c_double)),
        #     y.ctypes.data_as(POINTER(c_double)))
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
    
    ** Expect centered data (ie with mean removed) **

    Does not check if va is definite positive (on inversible 
    for that matter), so the inverse computation and/or determinant
    would throw an exception."""
    d       = mu.size
    inva    = lin.inv(va)
    fac     = 1 / N.sqrt( (2*N.pi) ** d * N.fabs(lin.det(va)))

    # we are using a trick with sum to "emulate" 
    # the matrix multiplication inva * x without any explicit loop
    y   = N.dot((x-mu), inva)
    y   = -0.5 * N.sum(y * (x-mu), 1)

    if not log:
        y   = fac * N.exp(y)
    else:
        y   = y + N.log(fac)
 
    return y

if __name__ == "__main__":
    #=========================================
    # Test accuracy between pure and C python
    #=========================================
    mu  = N.array([2.0, 3])
    va  = N.array([5.0, 3])

    # Generate a multivariate gaussian of mean mu and covariance va
    nframes = 1e4
    X       = randn(nframes, 2)
    Yc      = N.dot(N.diag(N.sqrt(va)), X.transpose())
    Yc      = Yc.transpose() + mu

    Y   = D.gauss_den(Yc, mu, va)
    Yt  = gauss_den(Yc, mu, va)

    print "Diff is " + str(N.sqrt(N.sum((Y-Yt) ** 2))/nframes/2)
