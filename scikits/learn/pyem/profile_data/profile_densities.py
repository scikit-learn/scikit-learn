import numpy as N
from numpy.random import randn

from numpy.ctypeslib import load_library, ndpointer
from ctypes import cdll, c_uint, c_int, c_double, POINTER

lib = load_library("blop.so", "file")

arg1    = ndpointer(dtype=N.float64)
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype=N.float64)
arg5    = ndpointer(dtype=N.float64)

lib.compute.argtypes    = [arg1, arg2, arg3, arg4, arg5]
lib.compute.restype     = c_int
# Compare computing per component likelihood for frame per row vs frame per column
def component_likelihood(x, mu, va, log = False):
    """expect one frame to be one row (rank 2). mu and var are rank 1 array."""
    d = mu.size

    return N.exp(N.sum((x - mu) ** 2, 1))

def component_likelihood2(x, mu, va, log = False):
    """expect one frame to be one column (rank 2). mu and var are rank 1 array."""
    d = mu.size

    y = (x[0] - mu[0]) ** 2
    for i in range(1, d):
        y += (x[i] - mu[i]) ** 2

    return N.exp(y)

def component_likelihood3(x, mu, va, log = False):
    """expect one frame to be one row (rank 2). mu and var are rank 1 array."""
    d = mu.size

    y = N.empty(x.shape[0], x.dtype)
    return lib.compute(x, x.shape[0], d, mu, y)

def bench(func, mode = 'diag'):
    d       = 30
    n       = 1e5
    niter   = 10

    print "Compute %d times densities, %d dimension, %d frames" % (niter, d, n)
    mu  = randn(d)
    va  = abs(randn(d))
    
    X   = randn(n, d)
    for i in range(niter):
        Y   = func(X, mu, va)

def bench2(func, mode = 'diag'):
    d       = 30
    n       = 1e5
    niter   = 10

    print "Compute %d times densities, %d dimension, %d frames" % (niter, d, n)
    mu  = randn(d)
    va  = abs(randn(d))
    
    X   = randn(d, n)
    for i in range(niter):
        Y   = func(X, mu, va)

def benchpy():
    bench(component_likelihood)

def benchpy3():
    bench(component_likelihood3)

def benchpy2():
    bench2(component_likelihood2)

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
    prof.runcall(benchpy2)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()

    profile_file    = 'gdenc.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(benchpy3)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()
