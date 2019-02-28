# coding: utf-8
"""
Benchmarks of polynomial features for dense matrices
implemented in 0.20.2 against implementation from PR #13290.
"""
# Authors: Xavier Dupré (benchmark)
# License: MIT

from time import time
from itertools import combinations, chain
from itertools import combinations_with_replacement as combinations_w_r

import numpy as np
from numpy.random import rand
from numpy.testing import assert_almost_equal
import matplotlib.pyplot as plt
import pandas

from sklearn.preprocessing import PolynomialFeatures


##############################
# Implementations to benchmark.
##############################

def _combinations(n_features, degree, interaction_only, include_bias):
    comb = (combinations if interaction_only else combinations_w_r)
    start = int(not include_bias)
    return chain.from_iterable(comb(range(n_features), i)
                               for i in range(start, degree + 1))


def fct_polynomial_features_0_20_2(X, degree, interaction_only, order):
    "PolynomialFeature 0.20.2 for dense matrices (no bias)."
    comb = _combinations(X.shape[1], degree, interaction_only,
                         include_bias=False)
    comb = list(comb)

    def compute_feat(X, comb=comb):
        XP = np.empty((X.shape[0], len(comb)), dtype=X.dtype, order=order)
        for i, comb in enumerate(comb):
            XP[:, i] = X[:, comb].prod(1)
        return XP

    return compute_feat


def fct_polynomial_features(X, degree, interaction_only, order):
    "Current implementation of PolynomialFeatures."
    model = PolynomialFeatures(degree=degree, include_bias=False,
                               interaction_only=interaction_only,
                               order=order)
    model.fit(X)

    def compute_feat(X, model=model):
        return model.transform(X)

    return compute_feat


##############################
# Benchmarks
##############################


def bench(n_obs, n_features, degrees, interactions_only, orders,
          repeat=10, verbose=False):
    res = []
    for n in n_obs:
        for nfeat in n_features:
            for order in orders:
                X = np.empty((n, nfeat), order=order)
                X[:, :] = rand(n, nfeat)[:, :]
                for degree in degrees:
                    for interaction_only in interactions_only:

                        obs = dict(n=n, nfeat=nfeat, order=order,
                                   interaction_only=interaction_only,
                                   degree=degree)
                        fct1 = fct_polynomial_features_0_20_2(
                            X, degree, interaction_only, order)
                        fct2 = fct_polynomial_features(
                            X, degree, interaction_only, order)

                        # creates different inputs to avoid caching in any ways
                        Xs = []
                        for r in range(repeat):
                            x = np.empty((n, nfeat), order=order)
                            x[:, :] = rand(n, nfeat)[:, :]
                            Xs.append(x)

                        # measures the baseline
                        st = time()
                        r = 0
                        for X in Xs:
                            p1 = fct1(X)
                            r += 1
                            if time() - st >= 1:
                                break  # stops if longer than a second
                        end = time()
                        obs["time_0_20_2"] = (end - st) / r

                        # measures the new implementation
                        st = time()
                        r2 = 0
                        for X in Xs:
                            p2 = fct2(X)
                            r2 += 1
                            if r2 >= r:
                                break
                        end = time()
                        obs["time_current"] = (end - st) / r
                        res.append(obs)
                        if verbose and (len(res) % 10 == 0 or n >= 10000):
                            print("bench", len(res), ":", obs)

                        # checks that both produce the same outputs
                        assert_almost_equal(p1, p2)
    return res


##############################
# Plots.
##############################

def plot_results(df, verbose=False):
    nrows = len(set(df.degree))
    fig, ax = plt.subplots(nrows, 4, figsize=(16, 12))
    pos = 0

    for di, degree in enumerate(sorted(set(df.degree))):
        pos = 0
        for order in sorted(set(df.order)):
            for interaction_only in sorted(set(df.interaction_only)):
                a = ax[di, pos]
                if di == ax.shape[0] - 1:
                    a.set_xlabel("N observations", fontsize='x-small')
                if pos == 0:
                    a.set_ylabel("Time (s) degree={}".format(degree),
                                 fontsize='x-small')

                for color, nfeat in zip('brgy', sorted(set(df.nfeat))):
                    subset = df[(df.degree == degree) & (df.nfeat == nfeat) &
                                (df.interaction_only == interaction_only) &
                                (df.order == order)]
                    subset = subset.sort_values("n")
                    if verbose:
                        print(subset)
                    label = "nf={} l=0.20.2".format(nfeat)
                    subset.plot(x="n", y="time_0_20_2", label=label, ax=a,
                                logx=True, logy=True, c=color, style='--')
                    label = "nf={} l=now".format(nfeat)
                    subset.plot(x="n", y="time_current", label=label, ax=a,
                                logx=True, logy=True, c=color)

                a.legend(loc=0, fontsize='x-small')
                if di == 0:
                    a.set_title("order={} interaction_only={}".format(
                        order, interaction_only), fontsize='x-small')
                pos += 1

    plt.suptitle("Benchmark for PolynomialFeatures", fontsize=16)


def run_bench(repeat=100, verbose=False):
    n_obs = [1, 10, 100, 1000, 10000]
    n_features = [10, 20]
    degrees = [2, 3]
    interactions_only = [False, True]
    orders = ['C', 'F']

    start = time()
    results = bench(n_obs, n_features, degrees, interactions_only, orders,
                    repeat=repeat, verbose=verbose)
    end = time()

    results_df = pandas.DataFrame(results)
    print("Total time = %0.3f sec\n" % (end - start))

    # plot the results
    plot_results(results_df, verbose=verbose)
    return results_df


if __name__ == '__main__':
    df = run_bench(verbose=True)
    plt.savefig("bench_polynomial_features.png")
    df.to_csv("bench_polynomial_features.csv", index=False)
    plt.show()
