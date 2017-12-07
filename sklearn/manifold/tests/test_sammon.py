import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn.metrics import euclidean_distances
from sklearn.manifold import Sammon
from sklearn.utils.testing import assert_raises


def test_sammon_precomputed():
    dists = 1.0 - np.eye(4)
    fitter = Sammon(5, dissimilarity="precomputed", verbose=1, eps=0)
    X = fitter.fit_transform(dists)
    sampled_dists = euclidean_distances(X)
    assert_array_almost_equal(sampled_dists, dists, decimal=3)


def test_sammon_euclidean():
    points = 1.0 - np.eye(3)
    fitter = Sammon(2, dissimilarity="euclidean")
    X = fitter.fit_transform(points)

    sampled_dists = euclidean_distances(X)
    true_dists = euclidean_distances(points)
    assert_array_almost_equal(sampled_dists, true_dists, decimal=3)


def test_sammon_other_metric():
    points = 1.0 - np.eye(3)
    fitter = Sammon(2, dissimilarity="baobabs")
    assert_raises(ValueError, lambda: fitter.fit_transform(points))


def test_sammon_init():
    points = 1.0 - np.eye(3)
    init = 1.0 - np.eye(4)
    fitter = Sammon(2, dissimilarity="baobabs")
    assert_raises(
        ValueError,
        lambda: fitter.fit_transform(points, init=init))


def test_sammon_init_2():
    points = 1.0 - np.eye(3)
    init = 1.0 - np.eye(3)
    fitter = Sammon(4, dissimilarity="euclidean", eps=0)
    X = fitter.fit_transform(points, init=init)

    sampled_dists = euclidean_distances(X)
    true_dists = euclidean_distances(points)
    assert_array_almost_equal(sampled_dists, true_dists, decimal=3)
