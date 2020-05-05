"""
Benchmarks on the power iterations phase in randomized SVD.

We test on various synthetic and real datasets the effect of increasing
the number of power iterations in terms of quality of approximation
and running time. A number greater than 0 should help with noisy matrices,
which are characterized by a slow spectral decay.

We test several policy for normalizing the power iterations. Normalization
is crucial to avoid numerical issues.

The quality of the approximation is measured by the spectral norm discrepancy
between the original input matrix and the reconstructed one (by multiplying
the randomized_svd's outputs). The spectral norm is always equivalent to the
largest singular value of a matrix. (3) justifies this choice. However, one can
notice in these experiments that Frobenius and spectral norms behave
very similarly in a qualitative sense. Therefore, we suggest to run these
benchmarks with `enable_spectral_norm = False`, as Frobenius' is MUCH faster to
compute.

The benchmarks follow.

(a) plot: time vs norm, varying number of power iterations
    data: many datasets
    goal: compare normalization policies and study how the number of power
    iterations affect time and norm

(b) plot: n_iter vs norm, varying rank of data and number of components for
    randomized_SVD
    data: low-rank matrices on which we control the rank
    goal: study whether the rank of the matrix and the number of components
    extracted by randomized SVD affect "the optimal" number of power iterations

(c) plot: time vs norm, varying datasets
    data: many datasets
    goal: compare default configurations

We compare the following algorithms:
-   randomized_svd(..., power_iteration_normalizer='none')
-   randomized_svd(..., power_iteration_normalizer='LU')
-   randomized_svd(..., power_iteration_normalizer='QR')
-   randomized_svd(..., power_iteration_normalizer='auto')
-   fbpca.pca() from https://github.com/facebook/fbpca (if installed)

Conclusion
----------
- n_iter=2 appears to be a good default value
- power_iteration_normalizer='none' is OK if n_iter is small, otherwise LU
  gives similar errors to QR but is cheaper. That's what 'auto' implements.

References
----------
(1) Finding structure with randomness: Stochastic algorithms for constructing
    approximate matrix decompositions
    Halko, et al., 2009 https://arxiv.org/abs/0909.4061

(2) A randomized algorithm for the decomposition of matrices
    Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert

(3) An implementation of a randomized algorithm for principal component
    analysis
    A. Szlam et al. 2014
"""

# Author: Giorgio Patrini

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import gc
import pickle
from time import time
from collections import defaultdict
import os.path

from sklearn.utils import gen_batches
from sklearn.utils.validation import check_random_state
from sklearn.utils.extmath import randomized_svd
from sklearn.datasets import make_low_rank_matrix, make_sparse_uncorrelated
from sklearn.datasets import (fetch_lfw_people,
                              fetch_openml,
                              fetch_20newsgroups_vectorized,
                              fetch_olivetti_faces,
                              fetch_rcv1)

try:
    import fbpca
    fbpca_available = True
except ImportError:
    fbpca_available = False

# If this is enabled, tests are much slower and will crash with the large data
enable_spectral_norm = False

# TODO: compute approximate spectral norms with the power method as in
# Estimating the largest eigenvalues by the power and Lanczos methods with
# a random start, Jacek Kuczynski and Henryk Wozniakowski, SIAM Journal on
# Matrix Analysis and Applications, 13 (4): 1094-1122, 1992.
# This approximation is a very fast estimate of the spectral norm, but depends
# on starting random vectors.

# Determine when to switch to batch computation for matrix norms,
# in case the reconstructed (dense) matrix is too large
MAX_MEMORY = np.int(2e9)

# The following datasets can be downloaded manually from:
# CIFAR 10: https://www.cs.toronto.edu/~kriz/cifar-10-python.tar.gz
# SVHN: http://ufldl.stanford.edu/housenumbers/train_32x32.mat
CIFAR_FOLDER = "./cifar-10-batches-py/"
SVHN_FOLDER = "./SVHN/"

datasets = ['low rank matrix', 'lfw_people', 'olivetti_faces', '20newsgroups',
            'mnist_784', 'CIFAR', 'a3a', 'SVHN', 'uncorrelated matrix']

big_sparse_datasets = ['big sparse matrix', 'rcv1']


def unpickle(file_name):
    with open(file_name, 'rb') as fo:
        return pickle.load(fo, encoding='latin1')["data"]


def handle_missing_dataset(file_folder):
    if not os.path.isdir(file_folder):
        print("%s file folder not found. Test skipped." % file_folder)
        return 0


def get_data(dataset_name):
    print("Getting dataset: %s" % dataset_name)

    if dataset_name == 'lfw_people':
        X = fetch_lfw_people().data
    elif dataset_name == '20newsgroups':
        X = fetch_20newsgroups_vectorized().data[:, :100000]
    elif dataset_name == 'olivetti_faces':
        X = fetch_olivetti_faces().data
    elif dataset_name == 'rcv1':
        X = fetch_rcv1().data
    elif dataset_name == 'CIFAR':
        if handle_missing_dataset(CIFAR_FOLDER) == "skip":
            return
        X1 = [unpickle("%sdata_batch_%d" % (CIFAR_FOLDER, i + 1))
              for i in range(5)]
        X = np.vstack(X1)
        del X1
    elif dataset_name == 'SVHN':
        if handle_missing_dataset(SVHN_FOLDER) == 0:
            return
        X1 = sp.io.loadmat("%strain_32x32.mat" % SVHN_FOLDER)['X']
        X2 = [X1[:, :, :, i].reshape(32 * 32 * 3) for i in range(X1.shape[3])]
        X = np.vstack(X2)
        del X1
        del X2
    elif dataset_name == 'low rank matrix':
        X = make_low_rank_matrix(n_samples=500, n_features=np.int(1e4),
                                 effective_rank=100, tail_strength=.5,
                                 random_state=random_state)
    elif dataset_name == 'uncorrelated matrix':
        X, _ = make_sparse_uncorrelated(n_samples=500, n_features=10000,
                                        random_state=random_state)
    elif dataset_name == 'big sparse matrix':
        sparsity = np.int(1e6)
        size = np.int(1e6)
        small_size = np.int(1e4)
        data = np.random.normal(0, 1, np.int(sparsity/10))
        data = np.repeat(data, 10)
        row = np.random.uniform(0, small_size, sparsity)
        col = np.random.uniform(0, small_size, sparsity)
        X = sp.sparse.csr_matrix((data, (row, col)), shape=(size, small_size))
        del data
        del row
        del col
    else:
        X = fetch_openml(dataset_name).data
    return X


def plot_time_vs_s(time, norm, point_labels, title):
    plt.figure()
    colors = ['g', 'b', 'y']
    for i, l in enumerate(sorted(norm.keys())):
        if l != "fbpca":
            plt.plot(time[l], norm[l], label=l, marker='o', c=colors.pop())
        else:
            plt.plot(time[l], norm[l], label=l, marker='^', c='red')

        for label, x, y in zip(point_labels, list(time[l]), list(norm[l])):
            plt.annotate(label, xy=(x, y), xytext=(0, -20),
                         textcoords='offset points', ha='right', va='bottom')
    plt.legend(loc="upper right")
    plt.suptitle(title)
    plt.ylabel("norm discrepancy")
    plt.xlabel("running time [s]")


def scatter_time_vs_s(time, norm, point_labels, title):
    plt.figure()
    size = 100
    for i, l in enumerate(sorted(norm.keys())):
        if l != "fbpca":
            plt.scatter(time[l], norm[l], label=l, marker='o', c='b', s=size)
            for label, x, y in zip(point_labels, list(time[l]), list(norm[l])):
                plt.annotate(label, xy=(x, y), xytext=(0, -80),
                             textcoords='offset points', ha='right',
                             arrowprops=dict(arrowstyle="->",
                                             connectionstyle="arc3"),
                             va='bottom', size=11, rotation=90)
        else:
            plt.scatter(time[l], norm[l], label=l, marker='^', c='red', s=size)
            for label, x, y in zip(point_labels, list(time[l]), list(norm[l])):
                plt.annotate(label, xy=(x, y), xytext=(0, 30),
                             textcoords='offset points', ha='right',
                             arrowprops=dict(arrowstyle="->",
                                             connectionstyle="arc3"),
                             va='bottom', size=11, rotation=90)

    plt.legend(loc="best")
    plt.suptitle(title)
    plt.ylabel("norm discrepancy")
    plt.xlabel("running time [s]")


def plot_power_iter_vs_s(power_iter, s, title):
    plt.figure()
    for l in sorted(s.keys()):
        plt.plot(power_iter, s[l], label=l, marker='o')
    plt.legend(loc="lower right", prop={'size': 10})
    plt.suptitle(title)
    plt.ylabel("norm discrepancy")
    plt.xlabel("n_iter")


def svd_timing(X, n_comps, n_iter, n_oversamples,
               power_iteration_normalizer='auto', method=None):
    """
    Measure time for decomposition
    """
    print("... running SVD ...")
    if method is not 'fbpca':
        gc.collect()
        t0 = time()
        U, mu, V = randomized_svd(X, n_comps, n_oversamples, n_iter,
                                  power_iteration_normalizer,
                                  random_state=random_state, transpose=False)
        call_time = time() - t0
    else:
        gc.collect()
        t0 = time()
        # There is a different convention for l here
        U, mu, V = fbpca.pca(X, n_comps, raw=True, n_iter=n_iter,
                             l=n_oversamples+n_comps)
        call_time = time() - t0

    return U, mu, V, call_time


def norm_diff(A, norm=2, msg=True):
    """
    Compute the norm diff with the original matrix, when randomized
    SVD is called with *params.

    norm: 2 => spectral; 'fro' => Frobenius
    """

    if msg:
        print("... computing %s norm ..." % norm)
    if norm == 2:
        # s = sp.linalg.norm(A, ord=2)  # slow
        value = sp.sparse.linalg.svds(A, k=1, return_singular_vectors=False)
    else:
        if sp.sparse.issparse(A):
            value = sp.sparse.linalg.norm(A, ord=norm)
        else:
            value = sp.linalg.norm(A, ord=norm)
    return value


def scalable_frobenius_norm_discrepancy(X, U, s, V):
    # if the input is not too big, just call scipy
    if X.shape[0] * X.shape[1] < MAX_MEMORY:
        A = X - U.dot(np.diag(s).dot(V))
        return norm_diff(A, norm='fro')

    print("... computing fro norm by batches...")
    batch_size = 1000
    Vhat = np.diag(s).dot(V)
    cum_norm = .0
    for batch in gen_batches(X.shape[0], batch_size):
        M = X[batch, :] - U[batch, :].dot(Vhat)
        cum_norm += norm_diff(M, norm='fro', msg=False)
    return np.sqrt(cum_norm)


def bench_a(X, dataset_name, power_iter, n_oversamples, n_comps):

    all_time = defaultdict(list)
    if enable_spectral_norm:
        all_spectral = defaultdict(list)
        X_spectral_norm = norm_diff(X, norm=2, msg=False)
    all_frobenius = defaultdict(list)
    X_fro_norm = norm_diff(X, norm='fro', msg=False)

    for pi in power_iter:
        for pm in ['none', 'LU', 'QR']:
            print("n_iter = %d on sklearn - %s" % (pi, pm))
            U, s, V, time = svd_timing(X, n_comps, n_iter=pi,
                                       power_iteration_normalizer=pm,
                                       n_oversamples=n_oversamples)
            label = "sklearn - %s" % pm
            all_time[label].append(time)
            if enable_spectral_norm:
                A = U.dot(np.diag(s).dot(V))
                all_spectral[label].append(norm_diff(X - A, norm=2) /
                                           X_spectral_norm)
            f = scalable_frobenius_norm_discrepancy(X, U, s, V)
            all_frobenius[label].append(f / X_fro_norm)

        if fbpca_available:
            print("n_iter = %d on fbca" % (pi))
            U, s, V, time = svd_timing(X, n_comps, n_iter=pi,
                                       power_iteration_normalizer=pm,
                                       n_oversamples=n_oversamples,
                                       method='fbpca')
            label = "fbpca"
            all_time[label].append(time)
            if enable_spectral_norm:
                A = U.dot(np.diag(s).dot(V))
                all_spectral[label].append(norm_diff(X - A, norm=2) /
                                           X_spectral_norm)
            f = scalable_frobenius_norm_discrepancy(X, U, s, V)
            all_frobenius[label].append(f / X_fro_norm)

    if enable_spectral_norm:
        title = "%s: spectral norm diff vs running time" % (dataset_name)
        plot_time_vs_s(all_time, all_spectral, power_iter, title)
    title = "%s: Frobenius norm diff vs running time" % (dataset_name)
    plot_time_vs_s(all_time, all_frobenius, power_iter, title)


def bench_b(power_list):

    n_samples, n_features = 1000, 10000
    data_params = {'n_samples': n_samples, 'n_features': n_features,
                   'tail_strength': .7, 'random_state': random_state}
    dataset_name = "low rank matrix %d x %d" % (n_samples, n_features)
    ranks = [10, 50, 100]

    if enable_spectral_norm:
        all_spectral = defaultdict(list)
    all_frobenius = defaultdict(list)
    for rank in ranks:
        X = make_low_rank_matrix(effective_rank=rank, **data_params)
        if enable_spectral_norm:
            X_spectral_norm = norm_diff(X, norm=2, msg=False)
        X_fro_norm = norm_diff(X, norm='fro', msg=False)

        for n_comp in [np.int(rank/2), rank, rank*2]:
            label = "rank=%d, n_comp=%d" % (rank, n_comp)
            print(label)
            for pi in power_list:
                U, s, V, _ = svd_timing(X, n_comp, n_iter=pi, n_oversamples=2,
                                        power_iteration_normalizer='LU')
                if enable_spectral_norm:
                    A = U.dot(np.diag(s).dot(V))
                    all_spectral[label].append(norm_diff(X - A, norm=2) /
                                               X_spectral_norm)
                f = scalable_frobenius_norm_discrepancy(X, U, s, V)
                all_frobenius[label].append(f / X_fro_norm)

    if enable_spectral_norm:
        title = "%s: spectral norm diff vs n power iteration" % (dataset_name)
        plot_power_iter_vs_s(power_iter, all_spectral, title)
    title = "%s: Frobenius norm diff vs n power iteration" % (dataset_name)
    plot_power_iter_vs_s(power_iter, all_frobenius, title)


def bench_c(datasets, n_comps):
    all_time = defaultdict(list)
    if enable_spectral_norm:
        all_spectral = defaultdict(list)
    all_frobenius = defaultdict(list)

    for dataset_name in datasets:
        X = get_data(dataset_name)
        if X is None:
            continue

        if enable_spectral_norm:
            X_spectral_norm = norm_diff(X, norm=2, msg=False)
        X_fro_norm = norm_diff(X, norm='fro', msg=False)
        n_comps = np.minimum(n_comps, np.min(X.shape))

        label = "sklearn"
        print("%s %d x %d - %s" %
              (dataset_name, X.shape[0], X.shape[1], label))
        U, s, V, time = svd_timing(X, n_comps, n_iter=2, n_oversamples=10,
                                   method=label)

        all_time[label].append(time)
        if enable_spectral_norm:
            A = U.dot(np.diag(s).dot(V))
            all_spectral[label].append(norm_diff(X - A, norm=2) /
                                       X_spectral_norm)
        f = scalable_frobenius_norm_discrepancy(X, U, s, V)
        all_frobenius[label].append(f / X_fro_norm)

        if fbpca_available:
            label = "fbpca"
            print("%s %d x %d - %s" %
                  (dataset_name, X.shape[0], X.shape[1], label))
            U, s, V, time = svd_timing(X, n_comps, n_iter=2, n_oversamples=2,
                                       method=label)
            all_time[label].append(time)
            if enable_spectral_norm:
                A = U.dot(np.diag(s).dot(V))
                all_spectral[label].append(norm_diff(X - A, norm=2) /
                                           X_spectral_norm)
            f = scalable_frobenius_norm_discrepancy(X, U, s, V)
            all_frobenius[label].append(f / X_fro_norm)

    if len(all_time) == 0:
        raise ValueError("No tests ran. Aborting.")

    if enable_spectral_norm:
        title = "normalized spectral norm diff vs running time"
        scatter_time_vs_s(all_time, all_spectral, datasets, title)
    title = "normalized Frobenius norm diff vs running time"
    scatter_time_vs_s(all_time, all_frobenius, datasets, title)


if __name__ == '__main__':
    random_state = check_random_state(1234)

    power_iter = np.linspace(0, 6, 7, dtype=int)
    n_comps = 50

    for dataset_name in datasets:
        X = get_data(dataset_name)
        if X is None:
            continue
        print(" >>>>>> Benching sklearn and fbpca on %s %d x %d" %
              (dataset_name, X.shape[0], X.shape[1]))
        bench_a(X, dataset_name, power_iter, n_oversamples=2,
                n_comps=np.minimum(n_comps, np.min(X.shape)))

    print(" >>>>>> Benching on simulated low rank matrix with variable rank")
    bench_b(power_iter)

    print(" >>>>>> Benching sklearn and fbpca default configurations")
    bench_c(datasets + big_sparse_datasets, n_comps)

    plt.show()
