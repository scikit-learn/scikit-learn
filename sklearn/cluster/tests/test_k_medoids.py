"""Testing for K-Medoids"""
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
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
