"""
Tests for HDBSCAN clustering algorithm
Shamelessly based on (i.e. ripped off from) the DBSCAN test code
"""
#import pickle
import numpy as np
from scipy.spatial import distance
from scipy import sparse
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import assert_not_in
from sklearn.cluster.hdbscan_ import HDBSCAN
from sklearn.cluster.hdbscan_ import hdbscan
from sklearn.cluster.tests.common import generate_clustered_data
from sklearn.metrics.pairwise import pairwise_distances

n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)


def test_hdbscan_distance_matrix():
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)
    
    labels = hdbscan(D, metric='precomputed')
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels) # ignore noise 
    assert_equal(n_clusters_1, n_clusters)

    labels = HDBSCAN(metric="precomputed").fit(D).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)
    
def test_hdbscan_feature_vector():    
    labels = hdbscan(X)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, n_clusters)

    labels = HDBSCAN().fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)

def test_hdbscan_no_clusters():
    labels = hdbscan(X, min_cluster_size=len(X)+1)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, 0)
    
    labels = HDBSCAN(min_cluster_size=len(X)+1).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, 0)
    
def test_hdbscan_callable_metric():
    # metric is the function reference, not the string key.
    metric = distance.euclidean

    labels = hdbscan(X, metric=metric)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, n_clusters)

    labels = HDBSCAN(metric=metric).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)

def test_hdbscan_input_lists():
    X = [[1., 2.], [3., 4.]]
    HDBSCAN().fit(X)  # must not raise exception

def test_hdbscan_badargs():
    assert_raises(ValueError,
                  hdbscan,
                  X='fail')
    assert_raises(ValueError,
                  hdbscan,
                  X=None)
    assert_raises(ValueError,
                  hdbscan,
                  X, min_cluster_size='fail')
    assert_raises(ValueError,
                  hdbscan,
                  X, min_samples='fail')
    assert_raises(ValueError,
                  hdbscan,
                  X, min_samples=-1)
    assert_raises(ValueError,
                  hdbscan,
                  X, metric='imperial')
    assert_raises(ValueError,
                  hdbscan,
                  X, metric=None) 
    assert_raises(TypeError,
                  hdbscan,
                  X, p=None)
    assert_raises(ValueError,
                  hdbscan,
                  X, p=-1)
    
    
### Probably not applicable now #########################
#def test_dbscan_sparse():
#def test_dbscan_balltree():
#def test_pickle():
#def test_dbscan_core_samples_toy():
#def test_boundaries():

