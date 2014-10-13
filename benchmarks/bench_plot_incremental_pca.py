"""
========================
IncrementalPCA benchmark
========================

Benchmarks for IncrementalPCA

"""

import numpy as np
import gc
from time import time
from collections import defaultdict
import matplotlib.pyplot as plt
from sklearn.datasets import fetch_lfw_people
from sklearn.decomposition import IncrementalPCA, RandomizedPCA, PCA


def plot_results(X, y, label):
    plt.plot(X, y, label=label, marker='o')


def benchmark(estimator, data):
    gc.collect()
    print("Benching %s" % estimator)
    t0 = time()
    estimator.fit(data)
    training_time = time() - t0
    data_t = estimator.transform(data)
    data_r = estimator.inverse_transform(data_t)
    reconstruction_error = np.mean(np.abs(data - data_r))
    return {'time': training_time, 'error': reconstruction_error}


def plot_feature_times(all_times, batch_size, all_components, data):
    plt.figure()
    plot_results(all_components, all_times['pca'], label="PCA")
    plot_results(all_components, all_times['ipca'],
                 label="IncrementalPCA, bsize=%i" % batch_size)
    plot_results(all_components, all_times['rpca'], label="RandomizedPCA")
    plt.legend(loc="upper left")
    plt.suptitle("Algorithm runtime vs. n_components\n \
                 LFW, size %i x %i" % data.shape)
    plt.xlabel("Number of components (out of max %i)" % data.shape[1])
    plt.ylabel("Time (seconds)")


def plot_feature_errors(all_errors, batch_size, all_components, data):
    plt.figure()
    plot_results(all_components, all_errors['pca'], label="PCA")
    plot_results(all_components, all_errors['ipca'],
                 label="IncrementalPCA, bsize=%i" % batch_size)
    plot_results(all_components, all_errors['rpca'], label="RandomizedPCA")
    plt.legend(loc="lower left")
    plt.suptitle("Algorithm error vs. n_components\n"
                 "LFW, size %i x %i" % data.shape)
    plt.xlabel("Number of components (out of max %i)" % data.shape[1])
    plt.ylabel("Mean absolute error")


def plot_batch_times(all_times, n_features, all_batch_sizes, data):
    plt.figure()
    plot_results(all_batch_sizes, all_times['pca'], label="PCA")
    plot_results(all_batch_sizes, all_times['rpca'], label="RandomizedPCA")
    plot_results(all_batch_sizes, all_times['ipca'], label="IncrementalPCA")
    plt.legend(loc="lower left")
    plt.suptitle("Algorithm runtime vs. batch_size for n_components %i\n \
                 LFW, size %i x %i" % (
                 n_features, data.shape[0], data.shape[1]))
    plt.xlabel("Batch size")
    plt.ylabel("Time (seconds)")


def plot_batch_errors(all_errors, n_features, all_batch_sizes, data):
    plt.figure()
    plot_results(all_batch_sizes, all_errors['pca'], label="PCA")
    plot_results(all_batch_sizes, all_errors['ipca'], label="IncrementalPCA")
    plt.legend(loc="lower left")
    plt.suptitle("Algorithm error vs. batch_size for n_components %i\n \
                 LFW, size %i x %i" % (
                 n_features, data.shape[0], data.shape[1]))
    plt.xlabel("Batch size")
    plt.ylabel("Mean absolute error")


def fixed_batch_size_comparison(data):
    all_features = [i.astype(int) for i in np.linspace(data.shape[1] // 10,
                                                       data.shape[1], num=5)]
    batch_size = 1000
    # Compare runtimes and error for fixed batch size
    all_times = defaultdict(list)
    all_errors = defaultdict(list)
    for n_components in all_features:
        pca = PCA(n_components=n_components)
        rpca = RandomizedPCA(n_components=n_components, random_state=1999)
        ipca = IncrementalPCA(n_components=n_components, batch_size=batch_size)
        results_dict = {k: benchmark(est, data) for k, est in [('pca', pca),
                                                               ('ipca', ipca),
                                                               ('rpca', rpca)]}

        for k in sorted(results_dict.keys()):
            all_times[k].append(results_dict[k]['time'])
            all_errors[k].append(results_dict[k]['error'])

    plot_feature_times(all_times, batch_size, all_features, data)
    plot_feature_errors(all_errors, batch_size, all_features, data)


def variable_batch_size_comparison(data):
    batch_sizes = [i.astype(int) for i in np.linspace(data.shape[0] // 10,
                                                      data.shape[0], num=10)]

    for n_components in [i.astype(int) for i in
                         np.linspace(data.shape[1] // 10,
                                     data.shape[1], num=4)]:
        all_times = defaultdict(list)
        all_errors = defaultdict(list)
        pca = PCA(n_components=n_components)
        rpca = RandomizedPCA(n_components=n_components, random_state=1999)
        results_dict = {k: benchmark(est, data) for k, est in [('pca', pca),
                                                               ('rpca', rpca)]}

        # Create flat baselines to compare the variation over batch size
        all_times['pca'].extend([results_dict['pca']['time']] *
                                len(batch_sizes))
        all_errors['pca'].extend([results_dict['pca']['error']] *
                                 len(batch_sizes))
        all_times['rpca'].extend([results_dict['rpca']['time']] *
                                 len(batch_sizes))
        all_errors['rpca'].extend([results_dict['rpca']['error']] *
                                  len(batch_sizes))
        for batch_size in batch_sizes:
            ipca = IncrementalPCA(n_components=n_components,
                                  batch_size=batch_size)
            results_dict = {k: benchmark(est, data) for k, est in [('ipca',
                                                                   ipca)]}
            all_times['ipca'].append(results_dict['ipca']['time'])
            all_errors['ipca'].append(results_dict['ipca']['error'])

        plot_batch_times(all_times, n_components, batch_sizes, data)
        # RandomizedPCA error is always worse (approx 100x) than other PCA
        # tests
        plot_batch_errors(all_errors, n_components, batch_sizes, data)

faces = fetch_lfw_people(resize=.2, min_faces_per_person=5)
# limit dataset to 5000 people (don't care who they are!)
X = faces.data[:5000]
n_samples, h, w = faces.images.shape
n_features = X.shape[1]

X -= X.mean(axis=0)
X /= X.std(axis=0)

fixed_batch_size_comparison(X)
variable_batch_size_comparison(X)
plt.show()
