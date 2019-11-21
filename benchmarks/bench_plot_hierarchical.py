from collections import defaultdict
from time import time

import numpy as np
from numpy import random as nr

from sklearn.cluster import AgglomerativeClustering


def compute_bench(samples_range, features_range):

    it = 0
    results = defaultdict(lambda: [])

    max_it = len(samples_range) * len(features_range)
    for n_samples in samples_range:
        for n_features in features_range:
            it += 1
            print('==============================')
            print('Iteration %03d of %03d' % (it, max_it))
            print('n_samples %05d; n_features %02d' % (n_samples, n_features))
            print('==============================')
            print()
            data = nr.randint(-50, 51, (n_samples, n_features))

            for linkage in ("single", "average", "complete", "ward"):
                print(linkage.capitalize())
                tstart = time()
                AgglomerativeClustering(
                    linkage=linkage,
                    n_clusters=10
                ).fit(data)

                delta = time() - tstart
                print("Speed: %0.3fs" % delta)
                print()

                results[linkage].append(delta)

    return results


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d  # noqa: F401
    import matplotlib.pyplot as plt

    samples_range = np.linspace(1000, 15000, 5).astype(np.int)
    features_range = np.linspace(2, 50, 5).astype(np.int)

    results = compute_bench(samples_range, features_range)

    max_time = max([max(i) for i in [t for (label, t) in results.items()]])

    fig = plt.figure('scikit-learn agglomerative clustering benchmark results')
    n = 1
    for c, (label, timings) in zip('brcy',
                                   sorted(results.items())):
        ax = fig.add_subplot(2, 2, n, projection='3d')
        ax.set_zlim3d(0.0, max_time * 1.1)

        X, Y = np.meshgrid(samples_range, features_range)
        Z = np.asarray(timings).reshape(samples_range.shape[0],
                                        features_range.shape[0])
        ax.plot_surface(X, Y, Z.T, cstride=1, rstride=1, color=c, alpha=0.5)
        ax.set_xlabel('n_samples')
        ax.set_ylabel('n_features')
        ax.set_zlabel('time (s)')
        ax.set_title(label.capitalize() + " linkage")

        n += 1

    plt.show()
