import numpy as N
from numpy.random import randn
from pyem import densities as D
from pyem import _c_densities as DC
import tables

def bench(func, mode = 'diag'):
    #===========================================
    # Diag Gaussian of dimension 20
    #===========================================
    d       = 30
    n       = 1e5
    niter   = 10

    print "Compute %d times densities, %d dimension, %d frames" % (niter, d, n)
    # Generate a model with k components, d dimensions
    mu  = randn(1, d)
    if mode == 'diag':
        va  = abs(randn(1, d))
    elif mode == 'full':
        va  = randn(d, d)
        va  = N.dot(va, va.transpose())

    X   = randn(n, d)
    for i in range(niter):
        Y   = func(X, mu, va)

def benchpy():
    bench(D.gauss_den)

def benchc():
    bench(DC.gauss_den)

def benchpyfull():
    bench(D.gauss_den, 'full')

def benchcfull():
    bench(DC.gauss_den, 'full')

if __name__ == "__main__":
    import hotshot, hotshot.stats
    profile_file    = 'gdenpy.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(benchpy)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()

    profile_file    = 'gdenc.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(benchc)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()
