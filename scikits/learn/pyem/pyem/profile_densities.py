import numpy as N
from numpy.random import randn
import densities as D
import _c_densities as DC
import tables

def bench(func, mode = 'diag'):
    #===========================================
    # Diag Gaussian of dimension 20
    #===========================================
    d       = 30
    n       = 1e5
    niter   = 100

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
        Y   = func(X, mu, va)
    
    # Check values
    h5file  = tables.openFile('diag.dat', "r")
    X   = h5file.root.input.read()
    mu  = h5file.root.mu.read()
    va  = h5file.root.va.read()
    Yt  = h5file.root.output.read()
    Y   = func(X, mu, va)

    try:
        N.testing.assert_array_almost_equal(Y, Yt) 
    except AssertionError:
        print N.sum(Y)
        print N.sqrt(N.sum((Y-Yt) **2)) / n
        raise "Not accurate !!!!"

def benchpy():
    bench(D.gauss_den)

def benchc():
    bench(DC.gauss_den)

def benchpyfull():
    bench(D.gauss_den, 'full')

def benchcfull():
    bench(DC.gauss_den, 'full')

if __name__ == "__main__":
    import profile
    import pstats
    profile.run('benchpy()', 'gdenprof')
    p = pstats.Stats('gdenprof')
    print p.sort_stats('cumulative').print_stats(20)
    profile.run('benchc()', 'gdenprof')
    p = pstats.Stats('gdenprof')
    print p.sort_stats('cumulative').print_stats(20)
