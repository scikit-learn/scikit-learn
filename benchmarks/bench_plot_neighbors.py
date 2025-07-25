"""
Plot the scaling of the nearest neighbors algorithms with k, D, and N
"""

from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

from sklearn import datasets, neighbors


def get_data(N, D, dataset="dense"):
    if dataset == "dense":
        np.random.seed(0)
        return np.random.random((N, D))
    elif dataset == "digits":
        X, _ = datasets.load_digits(return_X_y=True)
        i = np.argsort(X[0])[::-1]
        X = X[:, i]
        return X[:N, :D]
    else:
        raise ValueError("invalid dataset: %s" % dataset)


def barplot_neighbors(
    Nrange=2 ** np.arange(1, 11),
    Drange=2 ** np.arange(7),
    krange=2 ** np.arange(10),
    N=1000,
    D=64,
    k=5,
    leaf_size=30,
    dataset="digits",
):
    algorithms = ("kd_tree", "brute", "ball_tree")
    fiducial_values = {"N": N, "D": D, "k": k}

    # ------------------------------------------------------------
    # varying N
    N_results_build = {alg: np.zeros(len(Nrange)) for alg in algorithms}
    N_results_query = {alg: np.zeros(len(Nrange)) for alg in algorithms}

    for i, NN in enumerate(Nrange):
        print("N = %i (%i out of %i)" % (NN, i + 1, len(Nrange)))
        X = get_data(NN, D, dataset)
        for algorithm in algorithms:
            nbrs = neighbors.NearestNeighbors(
                n_neighbors=min(NN, k), algorithm=algorithm, leaf_size=leaf_size
            )
            t0 = time()
            nbrs.fit(X)
            t1 = time()
            nbrs.kneighbors(X)
            t2 = time()

            N_results_build[algorithm][i] = t1 - t0
            N_results_query[algorithm][i] = t2 - t1

    # ------------------------------------------------------------
    # varying D
    D_results_build = {alg: np.zeros(len(Drange)) for alg in algorithms}
    D_results_query = {alg: np.zeros(len(Drange)) for alg in algorithms}

    for i, DD in enumerate(Drange):
        print("D = %i (%i out of %i)" % (DD, i + 1, len(Drange)))
        X = get_data(N, DD, dataset)
        for algorithm in algorithms:
            nbrs = neighbors.NearestNeighbors(
                n_neighbors=k, algorithm=algorithm, leaf_size=leaf_size
            )
            t0 = time()
            nbrs.fit(X)
            t1 = time()
            nbrs.kneighbors(X)
            t2 = time()

            D_results_build[algorithm][i] = t1 - t0
            D_results_query[algorithm][i] = t2 - t1

    # ------------------------------------------------------------
    # varying k
    k_results_build = {alg: np.zeros(len(krange)) for alg in algorithms}
    k_results_query = {alg: np.zeros(len(krange)) for alg in algorithms}

    X = get_data(N, DD, dataset)

    for i, kk in enumerate(krange):
        print("k = %i (%i out of %i)" % (kk, i + 1, len(krange)))
        for algorithm in algorithms:
            nbrs = neighbors.NearestNeighbors(
                n_neighbors=kk, algorithm=algorithm, leaf_size=leaf_size
            )
            t0 = time()
            nbrs.fit(X)
            t1 = time()
            nbrs.kneighbors(X)
            t2 = time()

            k_results_build[algorithm][i] = t1 - t0
            k_results_query[algorithm][i] = t2 - t1

    plt.figure(figsize=(8, 11))

    for sbplt, vals, quantity, build_time, query_time in [
        (311, Nrange, "N", N_results_build, N_results_query),
        (312, Drange, "D", D_results_build, D_results_query),
        (313, krange, "k", k_results_build, k_results_query),
    ]:
        ax = plt.subplot(sbplt, yscale="log")
        plt.grid(True)

        tick_vals = []
        tick_labels = []

        bottom = 10 ** np.min(
            [min(np.floor(np.log10(build_time[alg]))) for alg in algorithms]
        )

        for i, alg in enumerate(algorithms):
            xvals = 0.1 + i * (1 + len(vals)) + np.arange(len(vals))
            width = 0.8

            c_bar = plt.bar(xvals, build_time[alg] - bottom, width, bottom, color="r")
            q_bar = plt.bar(xvals, query_time[alg], width, build_time[alg], color="b")

            tick_vals += list(xvals + 0.5 * width)
            tick_labels += ["%i" % val for val in vals]

            plt.text(
                (i + 0.02) / len(algorithms),
                0.98,
                alg,
                transform=ax.transAxes,
                ha="left",
                va="top",
                bbox=dict(facecolor="w", edgecolor="w", alpha=0.5),
            )

            plt.ylabel("Time (s)")

        ax.xaxis.set_major_locator(ticker.FixedLocator(tick_vals))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(tick_labels))

        for label in ax.get_xticklabels():
            label.set_rotation(-90)
            label.set_fontsize(10)

        title_string = "Varying %s" % quantity

        descr_string = ""

        for s in "NDk":
            if s == quantity:
                pass
            else:
                descr_string += "%s = %i, " % (s, fiducial_values[s])

        descr_string = descr_string[:-2]

        plt.text(
            1.01,
            0.5,
            title_string,
            transform=ax.transAxes,
            rotation=-90,
            ha="left",
            va="center",
            fontsize=20,
        )

        plt.text(
            0.99,
            0.5,
            descr_string,
            transform=ax.transAxes,
            rotation=-90,
            ha="right",
            va="center",
        )

        plt.gcf().suptitle("%s data set" % dataset.capitalize(), fontsize=16)

    plt.figlegend((c_bar, q_bar), ("construction", "N-point query"), "upper right")


if __name__ == "__main__":
    barplot_neighbors(dataset="digits")
    barplot_neighbors(dataset="dense")
    plt.show()
