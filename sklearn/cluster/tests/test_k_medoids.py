"""Testing for K-Medoids"""
import numpy as np

from sklearn.exceptions import NotFittedError
from sklearn.utils.testing import assert_equal, assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import raises

from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.datasets import load_iris

from sklearn.cluster import KMedoids, KMeans


rng = np.random.RandomState(0)
X = rng.rand(100, 5)

@raises(ValueError)
def test_kmedoids_fit_fails_n_clusters_is_zero():
    # n_clusters is 0
    KMedoids(n_clusters=0).fit(X)


@raises(ValueError)
def test_kmedoids_fit_fails_n_clusters_is_none():
    # n_clusters is None
    KMedoids(n_clusters=None).fit(X)


@raises(ValueError)
def test_kmedoids_fit_fails_bad_clustering_method():
    # Bad clustering
    KMedoids(n_clusters=5, clustering_method='foo').fit(X)


@raises(ValueError)
def test_kmedoids_fit_fails_init_is_none():
    KMedoids(init=None).fit(X)


def test_kmedoids_fit_succeeds():
    model = KMedoids().fit(X)

    assert_true(model is not None)

    model = KMedoids(n_clusters=5).fit(X)

    assert_true(model is not None)

    model = KMedoids(n_clusters=5, clustering_method='pam').fit(X)

    assert_true(model is not None)


def test_kmedoids_fit_fails_too_few_samples_vs_clusters():
    model = KMedoids(n_clusters=8)

    Xsmall = rng.rand(5, 2)

    # Trying to fit 3 samples to 8 clusters -> Wrong!
    assert_raises(ValueError, model.fit, Xsmall)


def test_kmedoids_fails_predict_before_fit():
    model = KMedoids()

    assert_raises(NotFittedError, model.predict, X)


def test_kmedoids_fails_transform_before_fit():
    model = KMedoids()

    assert_raises(NotFittedError, model.transform, X)


def test_kmedoids_fit_naive_with_all_pairwise_distance_functions():
    for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.values():
        if distance_metric is None:
            continue

        model = KMedoids(n_clusters=3, distance_metric=distance_metric)

        Xnaive = np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        model.fit(Xnaive)

        assert_array_equal(model.labels_, [0, 1, 2])


def test_kmedoids_iris_with_all_pairwise_distance_functions():
    np.random.seed(0)

    Xiris = load_iris()['data']

    refModel = KMeans(n_clusters=3)

    refModel.fit(Xiris)

    avgDistToClosestCentroid = \
        np.sum(np.min(euclidean_distances(Xiris, Y=refModel.cluster_centers_),
                      axis=1)) / Xiris.shape[0]

    for init in ['random', 'heuristic']:

        for (distance_metric_name,
             distance_metric) in PAIRWISE_DISTANCE_FUNCTIONS.items():

            if distance_metric is None:
                continue

            model = KMedoids(n_clusters=3,
                             distance_metric=distance_metric,
                             init=init)

            D = distance_metric(Xiris)

            avgDistToRandomMedoid = np.mean(D.ravel())

            model.fit(Xiris)

            avgDistToClosestMedoid = model.inertia(Xiris) / Xiris.shape[0]

            # We want distance-to-closest-medoid to be reduced from average
            # distance by more than 50%
            assert_greater(0.5 * avgDistToRandomMedoid,
                           avgDistToClosestMedoid)

            # When K-Medoids is using Euclidean distance,
            # we can compare its performance to
            # K-Means. We want to average distance to cluster centers
            # be similar between K-Means and K-Medoids
            if distance_metric_name == "euclidean":
                assert_greater(0.1,
                               np.abs(avgDistToClosestMedoid -
                                      avgDistToClosestCentroid))


def test_kmedoids_fit_predict():
    model = KMedoids()

    labels1 = model.fit_predict(X)

    assert_equal(len(labels1), 100)

    assert_array_equal(labels1, model.labels_)

    labels2 = model.predict(X)

    assert_array_equal(labels1, labels2)


def test_kmedoids_fit_transform():
    model = KMedoids()

    Xt1 = model.fit_transform(X)

    assert_array_equal(Xt1.shape, (100, model.n_clusters))

    Xt2 = model.transform(X)

    assert_array_equal(Xt1, Xt2)
