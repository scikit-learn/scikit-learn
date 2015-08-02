"""Testing for K-Medoids"""
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import raises

from exceptions import ValueError

from sklearn.cluster import KMedoids

@raises(ValueError)
def test_kmedoids_constructor_fails_n_clusters_is_zero():

    # n_clusters is 0
    model = KMedoids(n_clusters=0)


@raises(ValueError)
def test_kmedoids_constructor_fails_n_clusters_is_none():

    # n_clusters is None
    model = KMedoids(n_clusters=None)


@raises(ValueError)
def test_kmedoids_cunstructor_fails_bad_clustering():

    # Bad clustering
    model = KMedoids(n_clusters=5, clustering='foo')


def test_kmedoids_constructor_succeeds():

    model = KMedoids()

    assert_true(model is not None)

    model = KMedoids(n_clusters=5)

    assert_true(model is not None)

    model = KMedoids(n_clusters=5, clustering='pam')

    assert_true(model is not None)


def test_kmedoids_fit():


    model = KMedoids()

    X = np.random.rand(100,5)

    model.fit(X)
        

def test_kmedoids_fit_predict():

    model = KMedoids()

    X = np.random.rand(100,5)

    labels1 = model.fit_predict(X)
    
    assert_array_equal(labels1,model.labels_)
    
    labels2 = model.predict(X)

    assert_array_equal(labels1,labels2)


def test_kmedoids_fit_transform():

    model = KMedoids()

    X = np.random.rand(100,5)

    Xt1 = model.fit_transform(X)

    Xt2 = model.transform(X)

    assert_array_equal(Xt1,Xt2)
