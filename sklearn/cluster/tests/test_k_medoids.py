"""Testing for K-Medoids"""
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import raises

from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.datasets import load_iris

from sklearn.cluster import KMedoids, KMeans

from sklearn.utils.validation import NotFittedError


def get_dummy_X():

    np.random.seed(0)

    return np.random.rand(100, 5)


@raises(ValueError)
def test_kmedoids_fit_fails_n_clusters_is_zero():

    X = get_dummy_X()

    # n_clusters is 0
    KMedoids(n_clusters=0).fit(X)


@raises(ValueError)
def test_kmedoids_fit_fails_n_clusters_is_none():

    X = get_dummy_X()

    # n_clusters is None
    KMedoids(n_clusters=None).fit(X)


@raises(ValueError)
def test_kmedoids_fit_fails_bad_clustering_method():

    X = get_dummy_X()

    # Bad clustering
    KMedoids(n_clusters=5, clustering_method='foo').fit(X)


@raises(ValueError)
def test_kmedoids_fit_fails_init_is_none():

    X = get_dummy_X()

    KMedoids(init=None).fit(X)


def test_kmedoids_fit_succeeds():

    X = get_dummy_X()

    model = KMedoids().fit(X)

    assert_true(model is not None)

    model = KMedoids(n_clusters=5).fit(X)

    assert_true(model is not None)

    model = KMedoids(n_clusters=5, clustering_method='pam').fit(X)

    assert_true(model is not None)


@raises(ValueError)
def test_kmedoids_fit_fails_too_few_samples_vs_clusters():

    model = KMedoids(n_clusters=8)

    np.random.seed(0)

    X = np.random.rand(5, 2)

    # Trying to fit 3 samples to 8 clusters -> Wrong!
    model.fit(X)


@raises(NotFittedError)
def test_kmedoids_fails_predict_before_fit():

    X = get_dummy_X()

    model = KMedoids()

    model.predict(X)


@raises(NotFittedError)
def test_kmedoids_fails_transform_before_fit():

    X = get_dummy_X()

    model = KMedoids()

    model.transform(X)


def test_kmedoids_fit_naive():

    model = KMedoids(n_clusters=3)

    X = np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    model.fit(X)

    assert_array_equal(model.labels_, [0, 1, 2])


def test_kmedoids_fit_naive_with_all_pairwise_distance_functions():

    for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.keys():

        model = KMedoids(n_clusters=3, distance_metric=distance_metric)

        X = np.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        model.fit(X)

        assert_array_equal(model.labels_, [0, 1, 2])


def test_kmedoids_iris_with_all_pairwise_distance_functions():

    np.random.seed(0)

    X = load_iris()['data']

    refModel = KMeans(n_clusters=3)

    refModel.fit(X)

    avgDistToClosestCentroid = \
        np.sum(np.min(euclidean_distances(X, Y=refModel.cluster_centers_),
                      axis=1)) / X.shape[0]

    for init in ['random', 'heuristic']:

        for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.keys():

            model = KMedoids(n_clusters=3,
                             distance_metric=distance_metric,
                             init=init)

            D = PAIRWISE_DISTANCE_FUNCTIONS[distance_metric](X)

            avgDistToRandomMedoid = np.mean(D.ravel())

            model.fit(X)

            avgDistToClosestMedoid = model.inertia(X) / X.shape[0]

            # We want distance-to-closest-medoid to be reduced from average
            # distance by more than 50%
            assert_greater(0.5 * avgDistToRandomMedoid,
                           avgDistToClosestMedoid)

            # When K-Medoids is using Euclidean distance,
            # we can compare its performance to
            # K-Means. We want to average distance to cluster centers
            # be similar between K-Means and K-Medoids
            if distance_metric == "euclidean":
                assert_greater(0.1,
                               np.abs(avgDistToClosestMedoid -
                                      avgDistToClosestCentroid))


def test_kmedoids_fit_predict():

    model = KMedoids()

    X = np.random.rand(100, 5)

    labels1 = model.fit_predict(X)

    assert_equal(len(labels1), 100)

    assert_array_equal(labels1, model.labels_)

    labels2 = model.predict(X)

    assert_array_equal(labels1, labels2)


def test_kmedoids_fit_transform():

    model = KMedoids()

    X = np.random.rand(100, 5)

    Xt1 = model.fit_transform(X)

    assert_array_equal(Xt1.shape, (100, model.n_clusters))

    Xt2 = model.transform(X)

    assert_array_equal(Xt1, Xt2)
