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
    results = defaultdict(lambda: [])
    chunk = 100

    max_it = len(samples_range) * len(features_range)
    for n_samples in samples_range:
        for n_features in features_range:
            it += 1
            print '=============================='
            print 'Iteration %03d of %03d' %(it, max_it)
            print '=============================='
            print ''
            data = nr.random_integers(-50, 50, (n_samples, n_features))

            print 'K-Means'
            tstart = time()
            kmeans = KMeans(init='k-means++',
                            k=10).fit(data)

            delta = time() - tstart
            print "Speed: %0.3fs" % delta
            print "Inertia: %0.5f" % kmeans.inertia_

            results['kmeans_speed'].append(delta)
            results['kmeans_quality'].append(kmeans.inertia_ / \
                                        (n_samples * n_features))

            print 'Fast K-Means'
            # let's prepare the data in small chunks
            chunks_data = []
            tstart = time()
            mbkmeans = MiniBatchKMeans(init='k-means++',
                                      k=10)

            for i in xrange(iterations):
                j = i * chunk % len(data)
                chunks_data.append(data[j:j+chunk])
            tstart = time()
            for i in xrange(iterations):
                mbkmeans.partial_fit(chunks_data[i])
            mbkmeans.partial_fit(data)
            delta = time() - tstart
            print "Speed: %0.3fs" % delta
            print "Inertia: %f" % mbkmeans.inertia_

            results['minibatchkmeans_speed'].append(delta)
            results['minibatchkmeans_quality'].append(mbkmeans.inertia_ / \
                                        (n_samples * n_features))

    return results

def compute_bench_2(chunks, iterations):
    results = defaultdict(lambda: [])
    samples_range = 3
    features_range = 50000
    data = nr.random_integers(-50, 50, (n_samples, n_features))
    max_it = len(iterations) * len(chunks)
    for iteration in iterations:
        for chunk in chunks:
            it += 1
            print '=============================='
            print 'Iteration %03d of %03d' %(it, max_it)
            print '=============================='
            print ''

            print 'Fast K-Means'
            # let's prepare the data in small chunks
            chunks_data = []
            tstart = time()
            mbkmeans = MiniBatchKMeans(init='k-means++',
                                       k=10)

            for i in xrange(iteration):
                j = i * chunk % len(data)
                chunks_data.append(data[j:j+chunk])
            tstart = time()
            for i in xrange(iterations):
                mbkmeans.partial_fit(chunks_data[i])
            mbkmeans.partial_fit(data)
            delta = time() - tstart
            print "Speed: %0.3fs" % delta
            print "Inertia: %0.3fs" % mbkmeans.inertia_

            results['minibatchkmeans_speed'].append(delta)
            results['minibatchkmeans_quality'].append(mbkmeans.inertia_ / \
                                        (n_samples * n_features))

    return results



if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d # register the 3d projection
    import matplotlib.pyplot as plt

    samples_range = np.linspace(15, 150, 5).astype(np.int)
    features_range = np.linspace(150, 50000, 5).astype(np.int)
    chunks = np.linspace(50, 1250, 5).astype(np.int)
    iterations = np.linspace(100, 500, 5).astype(np.int)

    results = compute_bench(samples_range, features_range)
    results_2 = compute_bench(chunks, iterations)

    fig = plt.figure()
    for c, (label, timings) in zip('brcy',
                                    sorted(results.iteritems())):
        if 'speed' in label:
            ax = fig.add_subplot(2, 2, 1, projection='3d')
        else:
            ax = fig.add_subplot(2, 2, 2, projection='3d')

        X, Y = np.meshgrid(samples_range, features_range)
        Z = np.asarray(timings).reshape(samples_range.shape[0],
                                        features_range.shape[0])
        ax.plot_surface(X, Y, Z.T, cstride=1, rstride=1, color=c, alpha=0.5)

    i = 0
    for c, (label, timings) in zip('br',
                                   sorted(results_2.iteritems())):
        i += 1
        ax = fig.add_subplot(2, 2, i + 2, projection='3d')
        X, Y = np.meshgrid(samples_range, features_range)
        Z = np.asarray(timings).reshape(samples_range.shape[0],
                                        features_range.shape[0])
        ax.plot_surface(X, Y, Z.T, cstride=1, rstride=1, color=c, alpha=0.8)

    plt.show()
