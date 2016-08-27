from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster.optics import OPTICS
from sklearn.utils.testing import assert_equal, assert_greater_equal
from .common import generate_clustered_data
import numpy as np

def test_optics():
    '''
    Tests the optics clustering method and all functions inside it
    'auto' mode
    '''


    ##########################################################################

    n_clusters = 3
    X = generate_clustered_data(n_clusters=n_clusters)
    # Parameters chosen specifically for this task.
    # Compute OPTICS
    clust = OPTICS(eps=6.0, min_samples=4, metric='euclidean')
    clust.fit(X)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(clust.labels_)) - int(-1 in clust.labels_)
    assert_equal(n_clusters_1, n_clusters)

def test_filter():
    '''
    Tests the filter function.
    '''

    n_clusters = 3
    X = generate_clustered_data(n_clusters=n_clusters)
    # Parameters chosen specifically for this task.
    clust = OPTICS(eps=6.0, min_samples=4, metric='euclidean')
    # Run filter (before computing OPTICS)
    bool_memb = clust.filter(X, 0.5)
    idx_memb = clust.filter(X, 0.5, index_type='idx')
    # Test for equivalence between 'idx' and 'bool' extraction
    assert_equal(sum(bool_memb), len(idx_memb))
    # Compute OPTICS
    clust.fit(X)
    clust.extract(0.5,clustering='dbscan')
    # core points from filter and extract should be the same within 1 point,
    # with extract occasionally underestimating due to start point of the 
    # OPTICS algorithm. Here we test for at least 95% similarity in 
    # classification of core/not core
    agree = sum(clust._is_core == bool_memb)
    assert_greater_equal(float(agree)/len(X),0.95)

def test_optics2():
    '''
    Tests the optics clustering method and all functions inside it
    'dbscan' mode
    '''


    ##########################################################################
    # Compute OPTICS
    X = [[1,1]]
    clust = OPTICS(eps=0.3, min_samples=10)

    # Run the fit

    clust2 = clust.fit(X)

    assert clust2 == None
    # samples, labels = clust2.extract(0.4, 'dbscan')

    # assert samples.size == 0
    # assert labels[0] == -1


def test_empty_extract():
    '''
    Test extract where fit() has not yet been run.
    '''
    clust = OPTICS(eps=0.3, min_samples=10)
    assert clust.extract(0.01, clustering='auto') == None

def test_bad_extract():
    '''
    Test an extraction of eps too close to original eps
    '''
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers, 
                                cluster_std=0.4, random_state=0)

    ##########################################################################
    # Compute OPTICS

    clust = OPTICS(eps=0.003, min_samples=10)
    clust2 = clust.fit(X)
    assert clust2.extract(0.3) == None

def test_close_extract():
    '''
    Test extract where extraction eps is close to scaled epsPrime
    '''
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers, 
                                cluster_std=0.4, random_state=0)

    # Compute OPTICS

    clust = OPTICS(eps=0.2, min_samples=10)
    clust3 = clust.fit(X)
    clust3.extract(0.3, clustering='dbscan')
    assert max(clust3.labels_) == 3

def test_auto_extract_hier():
    # Generate sample data

    np.random.seed(0)
    n_points_per_cluster = 250

    X = np.empty((0, 2))
    X = np.r_[X, [-5,-2] + .8 * np.random.randn(n_points_per_cluster, 2)]
    X = np.r_[X, [4,-1] + .1 * np.random.randn(n_points_per_cluster, 2)]
    X = np.r_[X, [1,-2] + .2 * np.random.randn(n_points_per_cluster, 2)]
    X = np.r_[X, [-2,3] + .3 * np.random.randn(n_points_per_cluster, 2)]
    X = np.r_[X, [3,-2] + 1.6 * np.random.randn(n_points_per_cluster, 2)]
    X = np.r_[X, [5,6] + 2 * np.random.randn(n_points_per_cluster, 2)]
    # Compute OPTICS

    clust = OPTICS(eps=30.3, min_samples=9)

    # Run the fit
    clust.fit(X)

    # Extract the result
    clust.extract(0.0, 'auto') # eps not used for 'auto' extract

    assert_equal(len(set(clust.labels_)), 6)
