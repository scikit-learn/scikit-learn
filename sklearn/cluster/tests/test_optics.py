from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster.optics import OPTICS
from sklearn.utils.testing import assert_equal
from .common import generate_clustered_data

def test_optics():
    '''
    Tests the optics clustering method and all functions inside it
    'auto' mode
    '''


    ##############################################################################

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
    # Compute OPTICS
    clust = OPTICS(eps=6.0, min_samples=4, metric='euclidean')
    clust.fit(X)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(clust.labels_)) - int(-1 in clust.labels_)
    assert_equal(n_clusters_1, n_clusters)



def test_optics2():
    '''
    Tests the optics clustering method and all functions inside it
    'dbscan' mode
    '''


    ##############################################################################

    ##############################################################################
    # Compute OPTICS
    # Note the large eps; seeding problems when eps is close to
    # desired epsPrime 
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
    assert clust.extract(0.01) == None

def test_bad_extract():
    '''
    Test an extraction of an eps too close to original eps
    '''
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers, 
                                cluster_std=0.4, random_state=0)

    ##############################################################################

    ##############################################################################
    # Compute OPTICS

    clust = OPTICS(eps=0.003, min_samples=10)
    clust2 = clust.fit(X)
    assert clust2.extract(0.3) == None


