"""Benchmarks of Lasso regularization path computation using LARS and CD

The input data is mostly low rank but is a fat infinite tail.
"""
import gc
from time import time
import sys

from collections import defaultdict

import numpy as np
from numpy import random as nr

from scikits.learn.cluster.k_means_ import KMeans, MiniBatchKMeans

def compute_bench(samples_range, features_range):

    it = 0
    iterations = 200
    results = []

    max_it = len(samples_range) * len(features_range)
    for n_samples in samples_range:
        for n_features in features_range:
            it += 1
            print '=============================='
            print 'Iteration %03d of %03d' %(it, max_it)
            print '=============================='
            bkmeans = MiniBatchKMeans(init='k-means++',
                                  k=10)
            data = nr.rand(iterations, n_samples, n_features)
            tstart = time()
            for i in xrange(iterations):
                bkmeans.partial_fit(data[i])
            delta = time() - tstart
            print "%0.3fs" % delta
            results.append(delta)
    return results


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d # register the 3d projection
    import matplotlib.pyplot as plt

    samples_range = np.linspace(10, 500, 15).astype(np.int)
    features_range = np.linspace(10, 500, 15).astype(np.int)
    results = compute_bench(samples_range, features_range)

    max_time = max(results)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    i = 1
    X, Y = np.meshgrid(samples_range, features_range)
    Z = np.asarray(results).reshape(samples_range.shape[0],
                                    features_range.shape[0])
    ax.plot_surface(X, Y, Z.T, cstride=1, rstride=1, color='b', alpha=0.8)


#    for timings in results:
#        ax= fig.add_subplot(2, 2, i, projection='3d')
#        X, Y = np.meshgrid(samples_range, features_range)
#        Z = np.asarray(timings).reshape(samples_range.shape[0],
#                                        features_range.shape[0])
#
#        # plot the actual surface
#        ax.plot_surface(X, Y, Z.T, cstride=1, rstride=1, color=c, alpha=0.8)
#
#        # dummy point plot to stick the legend to since surface plot do not
#        # support legends (yet?)
#        #ax.plot([1], [1], [1], color=c, label=label)
#
#        ax.set_xlabel('n_samples')
#        ax.set_ylabel('n_features')
#        ax.set_zlabel('time (s)')
#        ax.set_zlim3d(0.0, max_time * 1.1)
#        ax.set_title(label)
#        #ax.legend()
#        i += 1
    plt.show()
