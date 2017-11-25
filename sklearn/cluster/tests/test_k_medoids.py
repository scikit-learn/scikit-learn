"""Testing for K-Medoids"""
import numpy as np
from scipy.sparse import csc_matrix

from sklearn.cluster import KMedoids, KMeans
from sklearn.datasets import load_iris
from sklearn.exceptions import NotFittedError
from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal, assert_raises, \
    assert_raise_message
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_allclose

rng = np.random.RandomState(0)
X = rng.rand(100, 5)


def test_kmedoids_input_validation_and_fit_check():
    # Invalid parameters
    assert_raise_message(ValueError, "n_clusters should be a nonnegative "
                                     "integer. 0 was given",
                         KMedoids(n_clusters=0).fit, X)

    assert_raise_message(ValueError, "n_clusters should be a nonnegative "
                                     "integer. None was given",
                         KMedoids(n_clusters=None).fit, X)

    assert_raise_message(ValueError, "init needs to be one of the following: "
                                     "['random', 'heuristic']",
                         KMedoids(init=None).fit, X)

    # Trying to fit 3 samples to 8 clusters
    Xsmall = rng.rand(5, 2)
    assert_raise_message(ValueError, "The number of medoids 8 must be less "
                                     "than the number of samples 5.",
                         KMedoids(n_clusters=8).fit, Xsmall)

    # Test if NotFittedError is raised appropriately
    assert_raises(NotFittedError, KMedoids().transform, X)
    assert_raises(NotFittedError, KMedoids().predict, X)


def test_kmedoids_fit_naive_with_all_pairwise_distance_functions():
    for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.keys():

        if distance_metric == 'precomputed':
            continue

        model = KMedoids(n_clusters=3, distance_metric=distance_metric)
        Xnaive = np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        model.fit(Xnaive)

        assert_array_equal(model.cluster_centers_,
                           [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        assert_array_equal(model.labels_, [0, 1, 2])
        assert_equal(model.inertia_, 0.0)


def test_kmedoids_iris_with_all_pairwise_distance_functions():
    X_iris = load_iris()['data']

    ref_model = KMeans(n_clusters=3)

    ref_model.fit(X_iris)

    avg_dist_to_closest_centroid = np.sum(np.min(
        euclidean_distances(X_iris, Y=ref_model.cluster_centers_), axis=1)
    ) / X_iris.shape[0]

    for init in ['random', 'heuristic']:
        for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.keys():
            if distance_metric == 'precomputed':
                continue

            model = KMedoids(n_clusters=3,
                             distance_metric=distance_metric,
                             init=init,
                             random_state=rng)
            model.fit(X_iris)

            distances = PAIRWISE_DISTANCE_FUNCTIONS[distance_metric](X_iris)
            avg_dist_to_random_medoid = np.mean(distances.ravel())
            avg_dist_to_closest_medoid = model.inertia_ / X_iris.shape[0]
            # We want distance-to-closest-medoid to be reduced from average
            # distance by more than 50%
            assert_greater(0.5 * avg_dist_to_random_medoid,
                           avg_dist_to_closest_medoid)
            # When K-Medoids is using Euclidean distance,
            # we can compare its performance to
            # K-Means. We want to average distance to cluster centers
            # be similar between K-Means and K-Medoids
            if distance_metric == "euclidean":
                assert_allclose(avg_dist_to_closest_medoid,
                                avg_dist_to_closest_centroid, rtol=0.1)


def test_kmedoids_fit_predict_transform():
    model = KMedoids(random_state=rng)

    labels1 = model.fit_predict(X)
    assert_equal(len(labels1), 100)
    assert_array_equal(labels1, model.labels_)

    labels2 = model.predict(X)
    assert_array_equal(labels1, labels2)

    Xt1 = model.fit_transform(X)
    assert_array_equal(Xt1.shape, (100, model.n_clusters))

    Xt2 = model.transform(X)
    assert_array_equal(Xt1, Xt2)


def test_outlier_robustness():
    kmeans = KMeans(n_clusters=2, random_state=rng)
    kmedoids = KMedoids(n_clusters=2, random_state=rng)

    X = [[-11, 0], [-10, 0], [-9, 0],
         [0, 0], [1, 0], [2, 0], [1000, 0]]

    kmeans.fit(X)
    kmedoids.fit(X)

    assert_array_equal(kmeans.labels_, [0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(kmedoids.labels_, [0, 0, 0, 1, 1, 1, 1])


def test_kmedoids_on_sparse_input():
    model = KMedoids(n_clusters=2, random_state=rng)
    row = np.array([1, 0])
    col = np.array([0, 4])
    data = np.array([1, 1])
    X = csc_matrix((data, (row, col)), shape=(2, 5))
    labels = model.fit_predict(X)
    assert_equal(len(labels), 2)
    assert_array_equal(labels, model.labels_)
