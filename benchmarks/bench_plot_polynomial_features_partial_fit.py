# coding: utf-8
"""
Benchmarks of polynomial features for dense matrices
implemented in 0.20.2 against implementation from PR #13290.
*PolynomialFeatures* is used in a pipeline.
"""
# Authors: Xavier Dupré (benchmark)
# License: MIT
from time import perf_counter as time
from itertools import combinations, chain
from itertools import combinations_with_replacement as combinations_w_r

import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
import pandas
from scipy import sparse
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import SGDClassifier
from sklearn.utils.testing import ignore_warnings
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted


##############################
# Implementations to benchmark.
##############################


def _combinations_poly(n_features, degree, interaction_only, include_bias):
    "Computes all polynomial features combinations."
    comb = (combinations if interaction_only else combinations_w_r)
    start = int(not include_bias)
    return chain.from_iterable(comb(range(n_features), i)
                               for i in range(start, degree + 1))


class CustomPolynomialFeatures(BaseEstimator, TransformerMixin):
    "Custom implemententations of PolynomialFeatures."
    def __init__(self, kind='poly-fast', degree=2, interaction_only=False,
                 include_bias=True, order='C'):
        BaseEstimator.__init__(self)
        TransformerMixin.__init__(self)
        self.kind = kind
        self.degree = degree
        self.include_bias = include_bias
        self.interaction_only = interaction_only
        self.order = order

    def fit(self, X, y=None):
        self.n_input_features_ = X.shape[1]
        poly = PolynomialFeatures(degree=self.degree,
                                  interaction_only=self.interaction_only,
                                  include_bias=self.include_bias,
                                  order=self.order)
        poly.fit(X)
        self.n_output_features_ = len(poly.get_feature_names())
        check_array(X, accept_sparse=False)
        return self

    def transform(self, X):
        check_is_fitted(self, ['n_input_features_', 'n_output_features_'])
        n_features = X.shape[1]
        if n_features != self.n_input_features_:
            raise ValueError("X shape does not match training shape")
        if self.kind == 'poly-fast':
            return self._transform_poly_fast(X)
        elif self.kind == 'poly-slow':
            return self._transform_poly_slow(X)
        else:
            raise ValueError(
                "Unknown extended features '{}'.".format(self.kind))

    def _transform_poly_fast(self, X):
        """
        Implementation in PR #13290.
        """
        if sparse.isspmatrix(X):
            raise NotImplementedError("Not implemented for sparse matrices.")

        n_samples, n_features = X.shape
        XP = np.empty((n_samples, self.n_output_features_),
                      dtype=X.dtype, order=self.order)

        if self.include_bias:
            XP[:, 0] = 1

        current_col = 1 if self.include_bias else 0
        for d in range(0, self.degree):
            if d == 0:
                XP[:, current_col:current_col + n_features] = X
                index = list(range(current_col,
                                   current_col + n_features))
                current_col += n_features
                index.append(current_col)
            else:
                new_index = []
                end = index[-1]
                for feature_idx in range(0, n_features):
                    a = index[feature_idx]
                    new_index.append(current_col)
                    start = a
                    if self.interaction_only:
                        start += index[feature_idx + 1] - \
                                 index[feature_idx]
                    next_col = current_col + end - start
                    if next_col <= current_col:
                        break
                    np.multiply(XP[:, start:end],
                                X[:, feature_idx:feature_idx + 1],
                                out=XP[:, current_col:next_col],
                                where=True, casting='no')
                    current_col = next_col

                new_index.append(current_col)
                index = new_index
        return XP

    def _transform_poly_slow(self, X):
        """
        Implemented in 0.20.2.
        """
        if sparse.isspmatrix(X):
            raise NotImplementedError("Not implemented for sparse matrices.")
        else:
            comb = _combinations_poly(X.shape[1], self.degree,
                                      self.interaction_only,
                                      include_bias=self.include_bias)
            XP = np.empty((X.shape[0], self.n_output_features_),
                          dtype=X.dtype, order=self.order)
            for i, comb in enumerate(comb):
                XP[:, i] = X[:, comb].prod(1)
            return XP


def measure_sgd_only(X_train, y_train, Xs, ys):
    model = SGDClassifier()
    model.fit(PolynomialFeatures().fit_transform(X_train), y_train)
    Xs = [PolynomialFeatures().fit_transform(x) for x in Xs]
    st = time()
    for X, y in zip(Xs, ys):
        model.partial_fit(X, y)
    end = time()
    return end - st


def measure_sgd_poly_skl(X_train, y_train, Xs, ys):
    model = make_pipeline(PolynomialFeatures(), SGDClassifier())
    model.fit(X_train, y_train)
    st = time()
    for X, y in zip(Xs, ys):
        X2 = model.steps[0][1].transform(X)
        model.steps[1][1].partial_fit(X2, y)
    end = time()
    return end - st


def measure_sgd_poly_fast(X_train, y_train, Xs, ys):
    model = make_pipeline(CustomPolynomialFeatures(kind='poly-fast'),
                          SGDClassifier())
    model.fit(X_train, y_train)
    st = time()
    for X, y in zip(Xs, ys):
        X2 = model.steps[0][1].transform(X)
        model.steps[1][1].partial_fit(X2, y)
    end = time()
    return end - st


def measure_sgd_poly_slow(X_train, y_train, Xs, ys):
    model = make_pipeline(CustomPolynomialFeatures(kind='poly-slow'),
                          SGDClassifier())
    model.fit(X_train, y_train)
    st = time()
    for X, y in zip(Xs, ys):
        X2 = model.steps[0][1].transform(X)
        model.steps[1][1].partial_fit(X2, y)
    end = time()
    return end - st


##############################
# Benchmarks
##############################


def build_x_y(ntrain, nfeat):
    """Creates a random classification problem with ntrain
    rows and nfeat features."""
    X_train = np.empty((ntrain, nfeat))
    X_train[:, :] = rand(ntrain, nfeat)[:, :]
    X_trainsum = X_train.sum(axis=1)
    eps = rand(ntrain) - 0.5
    X_trainsum_ = X_trainsum + eps
    y_train = (X_trainsum_ >= X_trainsum).ravel().astype(int)
    return X_train, y_train


@ignore_warnings(category=FutureWarning)
def bench(n_obs, n_features, repeat=1000, verbose=False, profiles=None):
    res = []
    for n in n_obs:
        for nfeat in n_features:

            X_train, y_train = build_x_y(1000, nfeat)
            obs = dict(n_obs=n, nfeat=nfeat)

            # creates different inputs to avoid caching in any ways
            Xs = []
            ys = []
            for r in range(repeat):
                X, y = build_x_y(n, nfeat)
                Xs.append(X)
                ys.append(y)

            t1 = measure_sgd_only(X_train, y_train, Xs, ys)
            t2 = measure_sgd_poly_skl(X_train, y_train, Xs, ys)
            t3 = measure_sgd_poly_fast(X_train, y_train, Xs, ys)
            t4 = measure_sgd_poly_slow(X_train, y_train, Xs, ys)

            obs["time_sgd"] = t1 / r
            obs["time_pipe_skl"] = t2 / r
            obs["time_pipe_fast"] = t3 / r
            obs["time_pipe_slow"] = t4 / r
            res.append(obs)

            if verbose and (len(res) % 1 == 0 or n >= 10000):
                print("bench", len(res), ":", obs)

    return res


##############################
# Plots.
##############################


def plot_results(df, verbose=False):
    nrows = max(len(set(df.nfeat)), 2)
    ncols = max(1, 2)
    fig, ax = plt.subplots(nrows, ncols,
                           figsize=(nrows * 4, ncols * 4))
    colors = "gbry"
    row = 0
    for nfeat in sorted(set(df.nfeat)):
        pos = 0
        for _ in range(1):
            a = ax[row, pos]
            if row == ax.shape[0] - 1:
                a.set_xlabel("N observations", fontsize='x-small')
            if pos == 0:
                a.set_ylabel("Time (s) nfeat={}".format(nfeat),
                             fontsize='x-small')

            subset = df[df.nfeat == nfeat]
            if subset.shape[0] == 0:
                continue
            subset = subset.sort_values("n_obs")
            if verbose:
                print(subset)

            label = "SGD-ONLY"
            subset.plot(x="n_obs", y="time_sgd", label=label, ax=a,
                        logx=True, logy=True, c=colors[0], style='--')
            label = "SGD-SKL"
            subset.plot(x="n_obs", y="time_pipe_skl", label=label, ax=a,
                        logx=True, logy=True, c=colors[1], style='--')
            label = "SGD-FAST"
            subset.plot(x="n_obs", y="time_pipe_fast", label=label, ax=a,
                        logx=True, logy=True, c=colors[2])
            label = "SGD-SLOW"
            subset.plot(x="n_obs", y="time_pipe_slow", label=label, ax=a,
                        logx=True, logy=True, c=colors[3])

            a.legend(loc=0, fontsize='x-small')
            if row == 0:
                a.set_title("--", fontsize='x-small')
            pos += 1
        row += 1

    plt.suptitle("Benchmark for PolynomialFeatures with SGDClassifier",
                 fontsize=16)


def run_bench(repeat=100, verbose=False):
    n_obs = [10, 100, 1000]
    n_features = [5, 10, 50]

    start = time()
    results = bench(n_obs, n_features, repeat=repeat, verbose=verbose,
                    profiles=[(100, 10)])
    end = time()

    results_df = pandas.DataFrame(results)
    print("Total time = %0.3f sec\n" % (end - start))

    # plot the results
    plot_results(results_df, verbose=verbose)
    return results_df


if __name__ == '__main__':
    df = run_bench(verbose=True)
    plt.savefig("bench_plot_polynomial_features_partial_fit.png")
    df.to_csv("bench_plot_polynomial_features_partial_fit.csv", index=False)
    plt.show()
