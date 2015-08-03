"""Testing for K-Medoids"""
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import raises

from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.datasets import load_iris

from exceptions import ValueError

from sklearn.cluster import KMedoids, KMeans

@raises(ValueError)
def test_kmedoids_constructor_fails_n_clusters_is_zero():

    # n_clusters is 0
    model = KMedoids(n_clusters=0)


@raises(ValueError)
def test_kmedoids_constructor_fails_n_clusters_is_none():

    # n_clusters is None
    model = KMedoids(n_clusters=None)


@raises(ValueError)
def test_kmedoids_constructor_fails_bad_clustering_method():

    # Bad clustering
    model = KMedoids(n_clusters=5, clustering_method='foo')


@raises(ValueError)
def test_kmedoids_constructor_fails_init_is_none():

    model = KMedoids(init=None)


def test_kmedoids_constructor_succeeds():

    model = KMedoids()

    assert_true(model is not None)

    model = KMedoids(n_clusters=5)

    assert_true(model is not None)

    model = KMedoids(n_clusters=5, clustering_method='pam')

    assert_true(model is not None)


@raises(ValueError)
def test_kmedoids_fit_fails_too_few_samples_vs_clusters():

    model = KMedoids(n_clusters=8)

    X = np.random.rand(5,2)

    # Trying to fit 3 samples to 8 clusters -> Wrong!
    model.fit(X)


def test_kmedoids_fit_naive():

    model = KMedoids(n_clusters=3)

    X = np.asarray([[1,0,0],[0,1,0],[0,0,1]])

    model.fit(X)

    assert_array_equal(model.labels_,[0,1,2])


def test_kmedoids_fit_naive_with_all_pairwise_distance_functions():

    for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.keys():

        model = KMedoids(n_clusters=3, distance_metric=distance_metric)

        X = np.asarray([[1,0,0],[0,1,0],[0,0,1]])

        model.fit(X)

        assert_array_equal(model.labels_,[0,1,2])

def test_kmedoids_iris_with_all_pairwise_distance_functions():

    X = load_iris()['data']

    refModel = KMeans(n_clusters=3)

    refModel.fit(X)
    
    avgDistanceToClosestCentroid = \
        np.sum(np.min(euclidean_distances(X, Y=refModel.cluster_centers_),
                      axis=1)) / X.shape[0]

    for distance_metric in PAIRWISE_DISTANCE_FUNCTIONS.keys():

        model = KMedoids(n_clusters=3, distance_metric=distance_metric)

        D = PAIRWISE_DISTANCE_FUNCTIONS[distance_metric](X)

        avgDistanceToRandomMedoid = np.mean(D.ravel())

        model.fit(X)

        avgDistanceToClosestMedoid = model.inertia(X) / X.shape[0]

        # We want distance-to-closest-medoid to be reduced from average 
        # distance by more than 50%
        assert_true(avgDistanceToClosestMedoid < 0.5*avgDistanceToRandomMedoid)

        # When K-Medoids is using Euclidean distance, 
        # we can compare its performance to 
        # K-Means. We want to average distance to cluster centers
        # be similar between K-Means and K-Medoids
        if distance_metric == "euclidean":
            assert_true(np.abs(avgDistanceToClosestMedoid -
                               avgDistanceToClosestCentroid) < 0.1)
        

def test_kmedoids_fit_predict():

    model = KMedoids()

    X = np.random.rand(100,5)

    labels1 = model.fit_predict(X)

    assert_equal(len(labels1),100)
    
    assert_array_equal(labels1,model.labels_)
    
    labels2 = model.predict(X)

    assert_array_equal(labels1,labels2)


def test_kmedoids_fit_transform():

    model = KMedoids()

    X = np.random.rand(100, 5)

    Xt1 = model.fit_transform(X)

    assert_array_equal(Xt1.shape, (100, model.n_clusters))

    Xt2 = model.transform(X)

    assert_array_equal(Xt1, Xt2)



    
