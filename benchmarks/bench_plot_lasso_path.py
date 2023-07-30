"""Benchmarks of Lasso regularization path computation using Lars and CD

The input data is mostly low rank but is a fat infinite tail.
"""
import gc
import sys
from collections import defaultdict
from time import time

import numpy as np

from sklearn.datasets import make_regression
from sklearn.linear_model import lars_path, lars_path_gram, lasso_path


def compute_bench(samples_range, features_range):
    it = 0

    results = defaultdict(lambda: [])

    max_it = len(samples_range) * len(features_range)
    for n_samples in samples_range:
        for n_features in features_range:
            it += 1
            print("====================")
            print("Iteration %03d of %03d" % (it, max_it))
            print("====================")
            dataset_kwargs = {
                "n_samples": n_samples,
                "n_features": n_features,
                "n_informative": n_features // 10,
                "effective_rank": min(n_samples, n_features) / 10,
                # 'effective_rank': None,
                "bias": 0.0,
            }
            print("n_samples: %d" % n_samples)
            print("n_features: %d" % n_features)
            X, y = make_regression(**dataset_kwargs)

            gc.collect()
            print("benchmarking lars_path (with Gram):", end="")
            sys.stdout.flush()
            tstart = time()
            G = np.dot(X.T, X)  # precomputed Gram matrix
            Xy = np.dot(X.T, y)
            lars_path_gram(Xy=Xy, Gram=G, n_samples=y.size, method="lasso")
            delta = time() - tstart
            print("%0.3fs" % delta)
            results["lars_path (with Gram)"].append(delta)

            gc.collect()
            print("benchmarking lars_path (without Gram):", end="")
            sys.stdout.flush()
            tstart = time()
            lars_path(X, y, method="lasso")
            delta = time() - tstart
            print("%0.3fs" % delta)
            results["lars_path (without Gram)"].append(delta)

            gc.collect()
            print("benchmarking lasso_path (with Gram):", end="")
            sys.stdout.flush()
            tstart = time()
            lasso_path(X, y, precompute=True)
            delta = time() - tstart
            print("%0.3fs" % delta)
            results["lasso_path (with Gram)"].append(delta)

            gc.collect()
            print("benchmarking lasso_path (without Gram):", end="")
            sys.stdout.flush()
            tstart = time()
            lasso_path(X, y, precompute=False)
            delta = time() - tstart
            print("%0.3fs" % delta)
            results["lasso_path (without Gram)"].append(delta)

    return results


if __name__ == "__main__":
    from mpl_toolkits.mplot3d import axes3d  # noqa register the 3d projection
    import matplotlib.pyplot as plt

    samples_range = np.linspace(10, 2000, 5).astype(int)
    features_range = np.linspace(10, 2000, 5).astype(int)
    results = compute_bench(samples_range, features_range)

    max_time = max(max(t) for t in results.values())

    fig = plt.figure("scikit-learn Lasso path benchmark results")
    i = 1
    for c, (label, timings) in zip("bcry", sorted(results.items())):
        ax = fig.add_subplot(2, 2, i, projection="3d")
        X, Y = np.meshgrid(samples_range, features_range)
        Z = np.asarray(timings).reshape(samples_range.shape[0], features_range.shape[0])

        # plot the actual surface
        ax.plot_surface(X, Y, Z.T, cstride=1, rstride=1, color=c, alpha=0.8)

        # dummy point plot to stick the legend to since surface plot do not
        # support legends (yet?)
        # ax.plot([1], [1], [1], color=c, label=label)

        ax.set_xlabel("n_samples")
        ax.set_ylabel("n_features")
        ax.set_zlabel("Time (s)")
        ax.set_zlim3d(0.0, max_time * 1.1)
        ax.set_title(label)
        # ax.legend()
        i += 1
    plt.show()
