"""
Plot the scaling of the nearest neighbors algorithms with k, D, and N
"""
from time import time

import numpy as np
import pylab as pl

from sklearn import neighbors, datasets


def get_data(N, D, dataset='dense'):
    if dataset == 'dense':
        np.random.seed(0)
        return np.random.random((N, D))
    elif dataset == 'digits':
        X = datasets.load_digits().data
        i = np.argsort(X[0])[::-1]
        X = X[:, i]
        return X[:N, :D]
    else:
        raise ValueError("invalid dataset: %s" % dataset)


def plot_neighbors_vs_N(Nrange=2 ** np.arange(11),
                        D=32,
                        k=5,
                        leaf_size=30,
                        dataset='digits',
                        ax=None):
    print '------------------------------------------------------------'
    print ' plot %s neighbors as a function of N:' % dataset
    print ' ', Nrange
    print '------------------------------------------------------------'

    results_build = {'ball_tree': np.zeros(len(Nrange)),
                     'kd_tree': np.zeros(len(Nrange)),
                     'brute': np.zeros(len(Nrange))}
    results_query = {'ball_tree': np.zeros(len(Nrange)),
                     'kd_tree': np.zeros(len(Nrange)),
                     'brute': np.zeros(len(Nrange))}

    for i, N in enumerate(Nrange):
        print "N = %i (%i out of %i)" % (N, i + 1, len(Nrange))
        X = get_data(N, D, dataset)
        for algorithm in results_build:
            nbrs = neighbors.NearestNeighbors(n_neighbors=min(k, N),
                                              algorithm=algorithm,
                                              leaf_size=leaf_size)
            t0 = time()
            nbrs.fit(X)
            t1 = time()
            nbrs.kneighbors(X)
            t2 = time()

            results_build[algorithm][i] = (t1 - t0)
            results_query[algorithm][i] = (t2 - t1)

    if ax is None:
        pl.figure()
        ax = pl.subplot(111)

    color_dict = {}

    for alg in results_build:
        l = ax.loglog(Nrange, results_build[alg] + results_query[alg],
                      label=alg)
        color_dict[alg] = l[0].get_color()

    for alg in results_build:
        ax.plot(Nrange, results_build[alg], ls='--', c=color_dict[alg])
        ax.plot(Nrange, results_query[alg], ls=':', c=color_dict[alg])

    pl.legend(loc=4)
    pl.xlabel('N')
    pl.ylabel('time (s)')
    pl.title('Time vs N for %s (D = %i, k = %i)' % (dataset, D, k))


def plot_neighbors_vs_D(Drange=2 ** np.arange(7),
                        N=1000,
                        k=5,
                        leaf_size=30,
                        dataset='digits',
                        ax=None):
    print '------------------------------------------------------------'
    print ' plot %s neighbors as a function of D:' % dataset
    print ' ', Drange
    print '------------------------------------------------------------'

    results_build = {'ball_tree': np.zeros(len(Drange)),
                     'kd_tree': np.zeros(len(Drange)),
                     'brute': np.zeros(len(Drange))}
    results_query = {'ball_tree': np.zeros(len(Drange)),
                     'kd_tree': np.zeros(len(Drange)),
                     'brute': np.zeros(len(Drange))}

    for i, D in enumerate(Drange):
        print "D = %i (%i out of %i)" % (D, i + 1, len(Drange))
        X = get_data(N, D, dataset)
        for algorithm in results_build:
            nbrs = neighbors.NearestNeighbors(n_neighbors=k,
                                              algorithm=algorithm,
                                              leaf_size=leaf_size)
            t0 = time()
            nbrs.fit(X)
            t1 = time()
            nbrs.kneighbors(X)
            t2 = time()

            results_build[algorithm][i] = (t1 - t0)
            results_query[algorithm][i] = (t2 - t1)

    if ax is None:
        pl.figure()
        ax = pl.subplot(111)

    color_dict = {}

    for alg in results_build:
        l = ax.loglog(Drange, results_build[alg] + results_query[alg],
                      label=alg)
        color_dict[alg] = l[0].get_color()

    for alg in results_build:
        ax.plot(Drange, results_build[alg], ls='--', c=color_dict[alg])
        ax.plot(Drange, results_query[alg], ls=':', c=color_dict[alg])

    pl.legend(loc=4)
    pl.xlabel('D')
    pl.ylabel('time (s)')
    pl.title('Time vs D for %s (N = %i, k = %i)' % (dataset, N, k))


def plot_neighbors_vs_k(krange=2 ** np.arange(10),
                        N=1000,
                        D=20,
                        leaf_size=30,
                        dataset='digits',
                        ax=None):
    print '------------------------------------------------------------'
    print ' plot %s neighbors as a function of k:' % dataset
    print ' ', krange
    print '------------------------------------------------------------'

    results_build = {'ball_tree': np.zeros(len(krange)),
                     'kd_tree': np.zeros(len(krange)),
                     'brute': np.zeros(len(krange))}
    results_query = {'ball_tree': np.zeros(len(krange)),
                     'kd_tree': np.zeros(len(krange)),
                     'brute': np.zeros(len(krange))}

    X = get_data(N, D, dataset)

    for i, k in enumerate(krange):
        print "k = %i (%i out of %i)" % (k, i + 1, len(krange))
        for algorithm in results_build:
            nbrs = neighbors.NearestNeighbors(n_neighbors=k,
                                              algorithm=algorithm,
                                              leaf_size=leaf_size)
            t0 = time()
            nbrs.fit(X)
            t1 = time()
            nbrs.kneighbors(X)
            t2 = time()

            results_build[algorithm][i] = (t1 - t0)
            results_query[algorithm][i] = (t2 - t1)

    if ax is None:
        pl.figure()
        ax = pl.subplot(111)

    color_dict = {}

    for alg in results_build:
        l = ax.loglog(krange, results_build[alg] + results_query[alg],
                      label=alg)
        color_dict[alg] = l[0].get_color()

    for alg in results_build:
        ax.plot(krange, results_build[alg], ls='--', c=color_dict[alg])
        ax.plot(krange, results_query[alg], ls=':', c=color_dict[alg])

    pl.legend(loc=4)
    pl.xlabel('k')
    pl.ylabel('time (s)')
    pl.title('Time vs k for %s (N = %i, D = %i)' % (dataset, N, D))

if __name__ == '__main__':
    for plot_func in [plot_neighbors_vs_N,
                      plot_neighbors_vs_D,
                      plot_neighbors_vs_k]:
        pl.figure(figsize=(8, 10))
        for dataset, plt in [('dense', 211), ('digits', 212)]:
            plot_func(dataset=dataset, ax=pl.subplot(plt))
            pl.grid(True)
            if plt == 211:
                pl.xlabel('')

    pl.show()
