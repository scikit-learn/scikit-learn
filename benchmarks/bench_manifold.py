# -*- coding: utf-8 -*-
"""
Benchmarks of the manifold module
"""
import gc
from time import time
import numpy as np

from scikits.learn.datasets.samples_generator import swissroll
from scikits.learn.manifold import LLE, LaplacianEigenmap
from scikits.learn.manifold import DiffusionMap, HessianMap

def bench(clf, X_train):
    gc.collect()

    # start time
    tstart = time()
    clf = clf.fit(X_train)
    delta = (time() - tstart)
    # stop time
    return delta

def compute_bench(n_samples):
    def LLEFactory():
        return LLE(n_coords=2, n_neighbors=8)

    def LaplacianEigenmapFactory():
        return LaplacianEigenmap(n_coords=2, n_neighbors=8)

    def DiffusionMapFactory():
        return DiffusionMap(n_coords=2)

    def HessianMapFactory():
        return HessianMap(n_coords=2, n_neighbors=8)

    lle_results = []
    laplacianeigenmap_results = []
    diffusionmap_results = []
    hessianmap_results = []

    n_tests = 1000
    it = 0

    for ns in n_samples:
        it += 1
        print '============'
        print 'Iteration %s' % it
        print '============'
        X, Y = swissroll(n_samples=ns)

        print "benching LLE: "
        lle_results.append(bench(LLEFactory(), X))
        print "benching Laplacian Eigenmaps: "
        laplacianeigenmap_results.append(bench(LaplacianEigenmapFactory(), X))
        print "benching Diffusion Maps: "
        diffusionmap_results.append(bench(DiffusionMapFactory(), X))
        print "benching Hessian Eigenmaps: "
        hessianmap_results.append(bench(HessianMapFactory(), X))

    return lle_results, laplacianeigenmap_results, diffusionmap_results, \
        hessianmap_results

if __name__ == '__main__':
    import pylab as pl

    alpha = 0.01 # regularization parameter

    n_features = 500
    list_n_samples = range(500, 3001, 500);
    lle_results, laplacianeigenmap_results, diffusionmap_results, \
        hessianmap_results = compute_bench(list_n_samples)

    pl.close('all')
    pl.title('Manifold benchmark')
    pl.plot(list_n_samples, lle_results, 'b-', label='LLE')
    pl.plot(list_n_samples, laplacianeigenmap_results,'r-', label='Laplacian Eigenmaps')
    pl.plot(list_n_samples, diffusionmap_results, 'g-', label='Diffusion Maps')
    pl.plot(list_n_samples, hessianmap_results,'y-', label='Hessian Eigenmaps')
    pl.legend()
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.show()

