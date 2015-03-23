"""Testing for Spectral Clustering methods"""

from sklearn.externals.six.moves import cPickle

dumps, loads = cPickle.dumps, cPickle.loads

import numpy as np
from scipy import sparse

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_warns_message

from sklearn.cluster import SpectralClustering, spectral_clustering
from sklearn.cluster.spectral import spectral_embedding
from sklearn.cluster.spectral import discretize
from sklearn.metrics import pairwise_distances
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.pairwise import kernel_metrics, rbf_kernel
from sklearn.datasets.samples_generator import make_blobs


def test_spectral_clustering():
    S = np.array([[1.0, 1.0, 1.0, 0.2, 0.0, 0.0, 0.0],
                  [1.0, 1.0, 1.0, 0.2, 0.0, 0.0, 0.0],
                  [1.0, 1.0, 1.0, 0.2, 0.0, 0.0, 0.0],
                  [0.2, 0.2, 0.2, 1.0, 1.0, 1.0, 1.0],
                  [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                  [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                  [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]])

    for eigen_solver in ('arpack', 'lobpcg'):
        for assign_labels in ('kmeans', 'discretize'):
            for mat in (S, sparse.csr_matrix(S)):
                model = SpectralClustering(random_state=0, n_clusters=2,
                                           affinity='precomputed',
                                           eigen_solver=eigen_solver,
                                           assign_labels=assign_labels
                                          ).fit(mat)
                labels = model.labels_
                if labels[0] == 0:
                    labels = 1 - labels

                assert_array_equal(labels, [1, 1, 1, 0, 0, 0, 0])

                model_copy = loads(dumps(model))
                assert_equal(model_copy.n_clusters, model.n_clusters)
                assert_equal(model_copy.eigen_solver, model.eigen_solver)
                assert_array_equal(model_copy.labels_, model.labels_)


def test_spectral_amg_mode():
    # Test the amg mode of SpectralClustering
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = pairwise_distances(X)  # Distance matrix
    S = np.max(D) - D  # Similarity matrix
    S = sparse.coo_matrix(S)
    try:
        from pyamg import smoothed_aggregation_solver

        amg_loaded = True
    except ImportError:
        amg_loaded = False
    if amg_loaded:
        labels = spectral_clustering(S, n_clusters=len(centers),
                                     random_state=0, eigen_solver="amg")
        # We don't care too much that it's good, just that it *worked*.
        # There does have to be some lower limit on the performance though.
        assert_greater(np.mean(labels == true_labels), .3)
    else:
        assert_raises(ValueError, spectral_embedding, S,
                      n_components=len(centers),
                      random_state=0, eigen_solver="amg")


def test_spectral_unknown_mode():
    # Test that SpectralClustering fails with an unknown mode set.
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = pairwise_distances(X)  # Distance matrix
    S = np.max(D) - D  # Similarity matrix
    S = sparse.coo_matrix(S)
    assert_raises(ValueError, spectral_clustering, S, n_clusters=2,
                  random_state=0, eigen_solver="<unknown>")


def test_spectral_unknown_assign_labels():
    # Test that SpectralClustering fails with an unknown assign_labels set.
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = pairwise_distances(X)  # Distance matrix
    S = np.max(D) - D  # Similarity matrix
    S = sparse.coo_matrix(S)
    assert_raises(ValueError, spectral_clustering, S, n_clusters=2,
                  random_state=0, assign_labels="<unknown>")


def test_spectral_clustering_sparse():
    X, y = make_blobs(n_samples=20, random_state=0,
                      centers=[[1, 1], [-1, -1]], cluster_std=0.01)

    S = rbf_kernel(X, gamma=1)
    S = np.maximum(S - 1e-4, 0)
    S = sparse.coo_matrix(S)

    labels = SpectralClustering(random_state=0, n_clusters=2,
                                affinity='precomputed').fit(S).labels_
    assert_equal(adjusted_rand_score(y, labels), 1)


def test_affinities():
    # Note: in the following, random_state has been selected to have
    # a dataset that yields a stable eigen decomposition both when built
    # on OSX and Linux
    X, y = make_blobs(n_samples=20, random_state=0,
                      centers=[[1, 1], [-1, -1]], cluster_std=0.01
                     )
    # nearest neighbors affinity
    sp = SpectralClustering(n_clusters=2, affinity='nearest_neighbors',
                            random_state=0)
    assert_warns_message(UserWarning, 'not fully connected', sp.fit, X)
    assert_equal(adjusted_rand_score(y, sp.labels_), 1)

    sp = SpectralClustering(n_clusters=2, gamma=2, random_state=0)
    labels = sp.fit(X).labels_
    assert_equal(adjusted_rand_score(y, labels), 1)

    X = check_random_state(10).rand(10, 5) * 10

    kernels_available = kernel_metrics()
    for kern in kernels_available:
        # Additive chi^2 gives a negative similarity matrix which
        # doesn't make sense for spectral clustering
        if kern != 'additive_chi2':
            sp = SpectralClustering(n_clusters=2, affinity=kern,
                                    random_state=0)
            labels = sp.fit(X).labels_
            assert_equal((X.shape[0],), labels.shape)

    sp = SpectralClustering(n_clusters=2, affinity=lambda x, y: 1,
                            random_state=0)
    labels = sp.fit(X).labels_
    assert_equal((X.shape[0],), labels.shape)

    def histogram(x, y, **kwargs):
        # Histogram kernel implemented as a callable.
        assert_equal(kwargs, {})    # no kernel_params that we didn't ask for
        return np.minimum(x, y).sum()

    sp = SpectralClustering(n_clusters=2, affinity=histogram, random_state=0)
    labels = sp.fit(X).labels_
    assert_equal((X.shape[0],), labels.shape)

    # raise error on unknown affinity
    sp = SpectralClustering(n_clusters=2, affinity='<unknown>')
    assert_raises(ValueError, sp.fit, X)


def test_discretize(seed=8):
    # Test the discretize using a noise assignment matrix
    random_state = np.random.RandomState(seed)
    for n_samples in [50, 100, 150, 500]:
        for n_class in range(2, 10):
            # random class labels
            y_true = random_state.random_integers(0, n_class, n_samples)
            y_true = np.array(y_true, np.float)
            # noise class assignment matrix
            y_indicator = sparse.coo_matrix((np.ones(n_samples),
                                             (np.arange(n_samples),
                                              y_true)),
                                            shape=(n_samples,
                                                   n_class + 1))
            y_true_noisy = (y_indicator.toarray()
                            + 0.1 * random_state.randn(n_samples,
                                                       n_class + 1))
            y_pred = discretize(y_true_noisy, random_state)
            assert_greater(adjusted_rand_score(y_true, y_pred), 0.8)
