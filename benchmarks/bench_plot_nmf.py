"""
Benchmarks of non-negative matrix factorization (NMF).
"""

from __future__ import print_function

from collections import defaultdict
from time import time

import numpy as np

from sklearn.decomposition.nmf import MultiplicativeNMF, ProjectedGradientNMF
from sklearn.datasets.samples_generator import make_low_rank_matrix


def benchmark(samples_range, features_range, rank=50, tolerance=1e-5):
    timeset = defaultdict(lambda: [])
    err = defaultdict(lambda: [])

    def record(model, name, time):
        loss = model.reconstruction_err_

        timeset[name].append(time)
        err[name].append(loss)

        print("Frobenius loss: %.5f" % loss)
        print("Elapsed time: %.2fs" % time)
        print()

    for n_samples in samples_range:
        for n_features in features_range:
            print("%2d samples, %2d features" % (n_samples, n_features))
            print('=======================')
            X = np.abs(make_low_rank_matrix(n_samples, n_features,
                                            effective_rank=rank,
                                            tail_strength=0.2))

            print("benchmarking nndsvd-nmf: ")
            tstart = time()
            m = ProjectedGradientNMF(n_components=30, tol=tolerance,
                                     init='nndsvd')
            m.fit(X)
            record(m, 'nndsvd-nmf', time() - tstart)
            del m

            print("benchmarking nndsvda-nmf: ")
            tstart = time()
            m = ProjectedGradientNMF(n_components=30, tol=tolerance,
                                     init='nndsvda')
            m.fit(X)
            record(m, 'nndsvda-nmf', time() - tstart)
            del m

            print("benchmarking nndsvdar-nmf: ")
            tstart = time()
            m = ProjectedGradientNMF(n_components=30, tol=tolerance,
                                     init='nndsvdar')
            m.fit(X)
            record(m, 'nndsvdar-nmf', time() - tstart)
            del m

            print("benchmarking random-nmf")
            tstart = time()
            m = ProjectedGradientNMF(n_components=30, tol=tolerance,
                                     init="random", random_state=31,
                                     max_iter=1000)
            m.fit(X)
            record(m, 'random-nmf', time() - tstart)
            del m

            print("benchmarking alt-random-nmf")
            tstart = time()
            m = MultiplicativeNMF(n_components=30, tol=tolerance,
                                  init="random")
            m.fit(X)
            record(m, 'alt-random-nmf', time() - tstart)
            del m

    return timeset, err


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d  # register the 3d projection
    axes3d
    import matplotlib.pyplot as plt

    samples_range = np.linspace(50, 500, 3).astype(np.int)
    features_range = np.linspace(50, 500, 3).astype(np.int)
    timeset, err = benchmark(samples_range, features_range)

    for i, results in enumerate((timeset, err)):
        fig = plt.figure('Non-negative matrix factorization benchmark')
        ax = fig.gca(projection='3d')
        for c, (label, timings) in zip('rbgcm', sorted(results.iteritems())):
            X, Y = np.meshgrid(samples_range, features_range)
            Z = np.asarray(timings).reshape(samples_range.shape[0],
                                            features_range.shape[0])
            # plot the actual surface
            ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3,
                            color=c)
            # dummy point plot to stick the legend to since surface plot do not
            # support legends (yet?)
            ax.plot([1], [1], [1], color=c, label=label)

        ax.set_xlabel('n_samples')
        ax.set_ylabel('n_features')
        zlabel = 'Time (s)' if i == 0 else 'reconstruction error'
        ax.set_zlabel(zlabel)
        ax.legend()
        plt.show()
