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
    x -= mu
    x **= 2
    return N.exp(N.dot(x, N.ones((mu.size, 1), x.dtype)))

def component_likelihood3(x, mu, va, log = False):
    """expect one frame to be one row (rank 2). mu and var are rank 1 array."""
    y = N.empty(x.shape[0], x.dtype)
    lib.compute(x, x.shape[0], x.shape[1], mu, y)
    return y

def bench(func, mode = 'diag'):
    d       = 30
    n       = 1e5
    niter   = 10

    print "Compute %d times densities, %d dimension, %d frames" % (niter, d, n)
    mu  = 0.1 * randn(d)
    va  = 0.1 * abs(randn(d))
    
    X   = 0.1 * randn(n, d)
    for i in range(niter):
        Y   = func(X, mu, va)

def benchpy():
    bench(component_likelihood)

def benchpy3():
    bench(component_likelihood3)

def benchpy2():
    bench2(component_likelihood2)

if __name__ == "__main__":
    #import hotshot, hotshot.stats
    #profile_file    = 'gdenpy.prof'
    #prof    = hotshot.Profile(profile_file, lineevents=1)
    #prof.runcall(benchpy)
    #p = hotshot.stats.load(profile_file)
    #print p.sort_stats('cumulative').print_stats(20)
    #prof.close()

    #profile_file    = 'gdenc.prof'
    #prof    = hotshot.Profile(profile_file, lineevents=1)
    #prof.runcall(benchpy3)
    #p = hotshot.stats.load(profile_file)
    #print p.sort_stats('cumulative').print_stats(20)
    #prof.close()

    #import cProfile as profile
    #profile.run('benchpy()', 'fooprof')
    benchpy()
    benchpy3()
