# Author Gordon Walsh (github: g-walsh) gordon.p.walsh@gmail.com
# Test of implementation of new gmm initialisation.
# Generate the responsibilities and means from existing 'kmeans' and 'random'
# also using a new 'rand_data' which selects a random sample of data points of
# size n_components.

# data generation code from https://jakevdp.github.io/

import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture


# Generate some data
from sklearn.datasets.samples_generator import make_blobs
X, y_true = make_blobs(n_samples=400, centers=4,
                       cluster_std=0.60, random_state=0)
X = X[:, ::-1]  # flip axes for better plotting

n_samples = 400
n_components = 4

# Plot the data with K Means Labels
# kmeans = KMeans(4, random_state=0)
# labels = kmeans.fit(X).predict(X)
# plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')
# plt.show()


def genGmm(init_params, seed=None):
    r = np.random.RandomState(seed)
    if init_params == 'kmeans':
        init_means = calculateMeans(kmeansMean(r))
    elif init_params == 'random':
        init_means = calculateMeans(randMean(r))
    elif init_params == 'rand_data':
        init_means = calculateMeans(randPointMean(r))
    else:
        raise ValueError("Unimplemented initialisation method '%s'"
                         % init_params)

    gmm = GaussianMixture(n_components=4, means_init=init_means, tol=1e-9,
                          max_iter=2000, random_state=r).fit(X)
    labels = gmm.predict(X)
    # plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')
    # plt.scatter(init_means[:, 0], init_means[:, 1], s=100,
    #             marker='x', c='black')
    return labels, init_means, seed, init_params
    # plt.show()


def kmeansMean(r):
    # Calculate the responsibilities by kmeans
    resp_km = np.zeros((n_samples, n_components))
    label = KMeans(n_clusters=n_components,
                   n_init=1, random_state=r).fit(X).labels_
    resp_km[np.arange(n_samples), label] = 1
    return resp_km
    # This will label all data points with one of the components absolutely.


def randPointMean(r):
    # Generate responsibilities to pick random points from the data.
    resp_select_point = np.zeros((n_samples, n_components))
    points = r.choice(range(n_samples), n_components, replace=False)
    for n, i in enumerate(points):
        resp_select_point[i, n] = 1
    return resp_select_point
    # This will label one random data point for each component. All others 0.


def randMean(r):
    # Generate random responsibilities for all points.
    resp_random_orig = r.rand(n_samples, n_components)
    resp_random_orig /= resp_random_orig.sum(axis=1)[:, np.newaxis]
    return resp_random_orig
    # This will label all points with random weighting.
    # Sum of responsibilities across a given point is 1.


def calculateMeans(resp):
    # Generate the means of the components. These are the initial parameters.
    nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
    means = np.dot(resp.T, X) / nk[:, np.newaxis]
    return means


# Tests and plots


def testSeed(seed):
    plt.figure(figsize=(9, 3))
    plt.subplots_adjust(bottom=.1, top=0.9, hspace=.15, wspace=.05,
                        left=.05, right=.95)
    methods = ['random', 'kmeans', 'rand_data']

    for n, i in enumerate(methods):
        plt.subplot(1, 3, n+1)
        labels, ini, seed, params = genGmm(i, seed)
        plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis', lw=0.5,
                    edgecolors='black')
        plt.scatter(ini[:, 0], ini[:, 1], s=150, marker='P', c='orange', lw=1.5,
                    edgecolors='black')
        plt.title(i)

    plt.show()
