#! /usr/bin/python
#
# Copyrighted David Cournapeau
# Last Change: Mon May 29 01:00 PM 2006 J

import numpy as N
import numpy.linalg as lin

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
    data    = x - mu

    if d == 1:
        # scalar case
        return _scalar_gauss_den(data[:, 0], mu[0, 0], va[0, 0], log)
    elif dv0 == 1:
        # Diagonal matrix case
        return _diag_gauss_den(data, mu, va, log)
    elif dv1 == dv0:
        # full case
        return  _full_gauss_den(data, mu, va, log)
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
    y       = (x ** 2) * -0.5 * inva
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
    y   = _scalar_gauss_den(x[:,0], mu[0,0], va[0,0], log)
    if not log:
        for i in range(1, d):
            y    *=  _scalar_gauss_den(x[:,i], mu[0,i], va[0,i], log)
    else:
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

    # # Slow version
    # n       = N.size(x, 0)
    # y       = N.zeros(n, float)
    # for i in range(n):
    #     y[i] = N.matrixmultiply(x[i,:],
    #              N.matrixmultiply(inva, N.transpose(x[i,:])))
    # y *= -0.5

    # we are using a trick with sum to "emulate" 
    # the matrix multiplication inva * x without any explicit loop
    y   = N.matrixmultiply(x, inva)
    y   = -0.5 * N.sum(y * x, 1)

    if not log:
        y   = fac * N.exp(y)
    else:
        y   = y + N.log(fac)
 
    return y

# To plot a confidence ellipse from multi-variate gaussian pdf
def gauss_ell(mu, va, dim = [0, 1], npoints = 100):
    """ Given a mean and covariance for multi-variate
    gaussian, returns npoints points for the ellipse
    of confidence 0.39
    
    Returns the coordinate x and y of the ellipse"""
    
    # TODO: Get a confidence interval using the Chi2 distribution
    # of points at a given mahalanobis distance...
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

    level   = 0.39
    
    # Generates a circle of npoints
    theta   = N.linspace(0, 2 * N.pi, npoints)
    circle  = N.array([N.cos(theta), N.sin(theta)])

    # Get the dimension which we are interested in:
    mu  = mu[dim]
    if mode == 'diag':
        va      = va[dim]
        elps    = N.outerproduct(mu, N.ones(npoints, float))
        elps    += N.matrixmultiply(N.diag(N.sqrt(va)), circle)
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
        elps    = N.outerproduct(mu, N.ones(npoints, float))
        elps    += N.matrixmultiply(cova, circle)
    else:
        raise DenParam("var mode not recognized")

    return elps[0, :], elps[1, :]

def test_gauss_den():
    """"""
    # import tables
    # import numpy as N
    # 
    # filename    = 'dendata.h5'

    # # # Dimension 1
    # # d   = 1
    # # mu  = 1.0
    # # va  = 2.0

    # # X   = N.randn(1e3, 1)

    # # Y   = gauss_den(X, mu, va)

    # # h5file      = tables.openFile(filename, "w")

    # # h5file.createArray(h5file.root, 'X', X)
    # # h5file.createArray(h5file.root, 'mu', mu)
    # # h5file.createArray(h5file.root, 'va', va)
    # # h5file.createArray(h5file.root, 'Y', Y)

    # # h5file.close()

    # # # Dimension 2, diag
    # # d   = 2
    # # mu  = N.array([1.0, -2.0])
    # # va  = N.array([1.0, 2.0])

    # # X   = N.randn(1e3, 2)

    # # Y   = gauss_den(X, mu, va)

    # # h5file      = tables.openFile(filename, "w")

    # # h5file.createArray(h5file.root, 'X', X)
    # # h5file.createArray(h5file.root, 'mu', mu)
    # # h5file.createArray(h5file.root, 'va', va)
    # # h5file.createArray(h5file.root, 'Y', Y)

    # # Dimension 2, full
    # d   = 2
    # mu  = N.array([[0.2, -1.0]])
    # va  = N.array([[1.2, 0.1], [0.1, 0.5]])

    # X   = N.randn(1e3, 2)

    # Y   = gauss_den(X, mu, va)

    # h5file      = tables.openFile(filename, "w")

    # h5file.createArray(h5file.root, 'X', X)
    # h5file.createArray(h5file.root, 'mu', mu)
    # h5file.createArray(h5file.root, 'va', va)
    # h5file.createArray(h5file.root, 'Y', Y)

    # h5file.close()

    import numpy.testing as testing
    #=================
    # Small test in 1d
    #=================
    va  = 2.0
    mu  = 1.0
    X   = N.linspace(-2, 2, 10)[:, N.NewAxis]

    Yt  = N.array([0.02973257230591, 0.05512079811082, 0.09257745306945, 
            0.14086453882683,
            0.19418015562214, 0.24250166773127, 0.27436665745048, 0.28122547107069,
            0.26114678964743, 0.21969564473386])

    Y   = gauss_den(X, mu, va)
    try:
        testing.assert_array_almost_equal(Y, Yt)
        print "1d test succeded"
    except AssertionError:
        print "test fails in 1d"

    #============================
    # Small test in 2d (diagonal)
    #============================
    mu  = N.atleast_2d([-1.0, 2.0])
    va  = N.atleast_2d([2.0, 3.0])
    X1  = N.linspace(-2, 2, 10)[:, N.NewAxis]
    X2  = N.linspace(-1, 3, 10)[:, N.NewAxis]
    X   = N.concatenate(([X1, X2]), 1)
    
    Yt  = N.array([0.01129091565384, 0.02025416837152, 0.03081845516786, 
            0.03977576221540, 0.04354490552910, 0.04043592581117, 
            0.03184994053539, 0.02127948225225, 0.01205937178755, 
            0.00579694938623 ])

    Y   = gauss_den(X, mu, va)
    try:
        testing.assert_array_almost_equal(Y, Yt)
        print "2d diag test succeded"
    except AssertionError:
        print "test fails in 2d diag"

    #============================
    # Small test in 2d (full mat)
    #============================
    mu  = N.array([[0.2, -1.0]])
    va  = N.array([[1.2, 0.1], [0.1, 0.5]])
    X1  = N.linspace(-2, 2, 10)[:, N.NewAxis]
    X2  = N.linspace(-3, 3, 10)[:, N.NewAxis]
    X   = N.concatenate(([X1, X2]), 1)
    
    Yt  = N.array([0.00096157109751, 0.01368908714856,
        0.07380823191162, 0.15072050533842, 
        0.11656739937861, 0.03414436965525,
        0.00378789836599, 0.00015915297541, 
        0.00000253261067, 0.00000001526368])

    Y   = gauss_den(X, mu, va)
    try:
        testing.assert_array_almost_equal(Y, Yt)
        print "2d full test succeded"
    except AssertionError:
        print "test fails in 2d full"
         

if __name__ == "__main__":
    import pylab

    #=========================================
    # Test plotting a simple diag 2d variance:
    #=========================================
    va  = N.array([5, 3])
    mu  = N.array([2, 3])

    # Generate a multivariate gaussian of mean mu and covariance va
    X       = N.randn(1e3, 2)
    Yc      = N.matrixmultiply(N.diag(N.sqrt(va)), X.transpose())
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
    X       = N.randn(1e3, 2)
    Yc      = N.matrixmultiply(lin.cholesky(va), X.transpose())
    Yc      = Yc.transpose() + mu

    # Plotting
    Xe, Ye  = gauss_ell(mu, va, npoints = 100)
    pylab.figure()
    pylab.plot(Yc[:, 0], Yc[:, 1], '.')
    pylab.plot(Xe, Ye, 'r')
    pylab.show()

    savefig('example.png')
