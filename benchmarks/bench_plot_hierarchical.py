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
            print("==============================")
            print("Iteration %03d of %03d" % (it, max_it))
            print("n_samples %05d; n_features %02d" % (n_samples, n_features))
            print("==============================")
            print()
            data = nr.randint(-50, 51, (n_samples, n_features))

            for linkage in ("single", "average", "complete", "ward"):
                print(linkage.capitalize())
                tstart = time()
                AgglomerativeClustering(linkage=linkage, n_clusters=10).fit(data)

                delta = time() - tstart
                print("Speed: %0.3fs" % delta)
                print()

                results[linkage].append(delta)

    return results


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    samples_range = np.linspace(1000, 15000, 8).astype(int)
    features_range = np.array([2, 10, 20, 50])

    results = compute_bench(samples_range, features_range)

    max_time = max([max(i) for i in [t for (label, t) in results.items()]])

    colors = plt.get_cmap("tab10")(np.linspace(0, 1, 10))[:4]
    lines = {linkage: None for linkage in results.keys()}
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
    fig.suptitle("Scikit-learn agglomerative clustering benchmark results", fontsize=16)
    for c, (label, timings) in zip(colors, sorted(results.items())):
        timing_by_samples = np.asarray(timings).reshape(
            samples_range.shape[0], features_range.shape[0]
        )

        for n in range(timing_by_samples.shape[1]):
            ax = axs.flatten()[n]
            (lines[label],) = ax.plot(
                samples_range, timing_by_samples[:, n], color=c, label=label
            )
            ax.set_title("n_features = %d" % features_range[n])
            if n >= 2:
                ax.set_xlabel("n_samples")
            if n % 2 == 0:
                ax.set_ylabel("time (s)")

    fig.subplots_adjust(right=0.8)
    fig.legend(
        [lines[link] for link in sorted(results.keys())],
        sorted(results.keys()),
        loc="center right",
        fontsize=8,
    )

    plt.show()
