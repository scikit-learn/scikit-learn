import numpy as N
from numpy.random import randn
import densities as D
import tables

def bench1(mode = 'diag'):
    #===========================================
    # Diag Gaussian of dimension 20
    #===========================================
    d       = 20
    n       = 1e4
    niter   = 20
    mode    = 'diag'

    # Generate a model with k components, d dimensions
    mu  = randn(1, d)
    if mode == 'diag':
        va  = abs(randn(1, d))
    elif mode == 'full':
        va  = randn(d, d)
        va  = N.dot(va, va.transpose())

    X   = randn(n, d)
    print "Compute %d times densities, %d dimension, %d frames" % (niter, d, n)
    for i in range(niter):
        Y   = D.gauss_den(X, mu, va)
    
    # Check values
    h5file  = tables.openFile('diag.dat', "r")
    X   = h5file.root.input.read()
    mu  = h5file.root.mu.read()
    va  = h5file.root.va.read()
    Yt  = h5file.root.output.read()
    Y   = D.gauss_den(X, mu, va)

    try:
        N.testing.assert_array_almost_equal(Y, Yt) 
    except AssertionError:
        print N.sqrt(N.sum((Y-Yt) **2)) / n
        raise "Not accurate !!!!"

if __name__ == "__main__":
    import profile
    profile.run('bench1()', 'gdenprof')
    import pstats
    p = pstats.Stats('gdenprof')
    print p.sort_stats('cumulative').print_stats(20)

