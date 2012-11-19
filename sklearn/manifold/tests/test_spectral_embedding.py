from nose.tools import assert_true
from nose.tools import assert_equal

from scipy import sparse
from scipy.sparse import csr_matrix
import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from nose.tools import assert_raises
from nose.plugins.skip import SkipTest

from sklearn.manifold.spectral_embedding import SpectralEmbedding
from sklearn.manifold.spectral_embedding import _graph_is_connected
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.pipeline import Pipeline
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.datasets.samples_generator import make_blobs


# non centered, sparse centers to check the
centers = np.array([
    [0.0, 5.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 4.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 5.0, 1.0],
])
n_samples = 1000
n_clusters, n_features = centers.shape
S, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                            cluster_std=1., random_state=42)


def test_spectral_embedding_two_components(seed=36):
    """Test spectral embedding with two components"""
    random_state = np.random.RandomState(seed)
    n_sample = 10
    affinity = np.zeros(shape=[n_sample * 2,
                               n_sample * 2])
    # first component
    affinity[0:n_sample,
             0:n_sample] = np.abs(random_state.randn(n_sample, n_sample)) + 2
    # second component
    affinity[n_sample::,
             n_sample::] = np.abs(random_state.randn(n_sample, n_sample)) + 2
    # connection
    affinity[0, n_sample + 1] = 1
    affinity[n_sample + 1, 0] = 1
    affinity.flat[::2 * n_sample + 1] = 0
    affinity = 0.5 * (affinity + affinity.T)

    true_label = np.zeros(shape=2 * n_sample)
    true_label[0:n_sample] = 1

    se_precomp = SpectralEmbedding(n_components=1, affinity="precomputed",
                                   random_state=np.random.RandomState(seed))
    embedded_corrdinate = np.squeeze(se_precomp.fit_transform(affinity))
    # thresholding on the first components using 0.
    label_ = np.array(embedded_corrdinate < 0, dtype="float")
    assert_equal(normalized_mutual_info_score(true_label, label_), 1.0)


def test_spectral_embedding_precomputed_affinity(seed=36):
    """Test spectral embedding with precomputed kernel"""
    gamma = 1.0
    se_precomp = SpectralEmbedding(n_components=3, affinity="precomputed",
                                   random_state=np.random.RandomState(seed))
    se_rbf = SpectralEmbedding(n_components=3, affinity="rbf",
                               gamma=gamma,
                               random_state=np.random.RandomState(seed))
    embed_precomp = se_precomp.fit_transform(rbf_kernel(S, gamma=gamma))
    embed_rbf = se_rbf.fit_transform(S)
    assert_array_almost_equal(
        se_precomp.affinity_matrix_, se_rbf.affinity_matrix_)
    assert_array_almost_equal(np.abs(embed_precomp), np.abs(embed_rbf), 0)


def test_spectral_embedding_callable_affinity(seed=36):
    """Test spectral embedding with callable affinity"""
    gamma = 0.9
    kern = rbf_kernel(S, gamma=gamma)
    se_callable = SpectralEmbedding(n_components=3,
                                    affinity=(
                                        lambda x: rbf_kernel(x, gamma=gamma)),
                                    gamma=gamma,
                                    random_state=np.random.RandomState(seed))
    se_rbf = SpectralEmbedding(n_components=3, affinity="rbf",
                               gamma=gamma,
                               random_state=np.random.RandomState(seed))
    embed_rbf = se_rbf.fit_transform(S)
    embed_callable = se_callable.fit_transform(S)
    embed_rbf = se_rbf.fit_transform(S)
    embed_callable = se_callable.fit_transform(S)
    assert_array_almost_equal(
        se_callable.affinity_matrix_, se_rbf.affinity_matrix_)
    assert_array_almost_equal(np.abs(embed_rbf), np.abs(embed_callable), 2)


def test_spectral_embedding_amg_solver(seed=36):
    """Test spectral embedding with amg solver"""
    try:
        from pyamg import smoothed_aggregation_solver
    except ImportError:
        raise SkipTest

    gamma = 0.9
    se_amg = SpectralEmbedding(n_components=3, affinity="rbf",
                               gamma=gamma, eig_solver="amg",
                               random_state=np.random.RandomState(seed))
    se_arpack = SpectralEmbedding(n_components=3, affinity="rbf",
                                  gamma=gamma, eig_solver="arpack",
                                  random_state=np.random.RandomState(seed))
    embed_amg = se_amg.fit_transform(S)
    embed_arpack = se_arpack.fit_transform(S)
    assert_array_almost_equal(
        se_amg.affinity_matrix_, se_arpack.affinity_matrix_)
    assert_array_almost_equal(np.abs(embed_amg), np.abs(embed_arpack), 2)


def test_pipline_spectral_clustering(seed=36):
    """Test using pipline to do spectral clustering"""
    random_state = np.random.RandomState(seed)
    se_rbf = SpectralEmbedding(n_components=n_clusters,
                               affinity="rbf",
                               random_state=random_state)
    se_knn = SpectralEmbedding(n_components=n_clusters,
                               affinity="nearest_neighbors",
                               n_neighbors=5,
                               random_state=random_state)
    for se in [se_rbf, se_knn]:
        km = KMeans(n_clusters=n_clusters, random_state=random_state)
        km.fit(se.fit_transform(S))
        assert_array_almost_equal(
            normalized_mutual_info_score(
                km.labels_,
                true_labels), 1.0, 2)


def test_spectral_embedding_unknown_eigensolver(seed=36):
    """Test that SpectralClustering fails with an unknown eigensolver"""
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = rbf_kernel(X)  # Distance matrix

    se_precomp = SpectralEmbedding(n_components=1, affinity="precomputed",
                                   random_state=np.random.RandomState(seed),
                                   eig_solver="<unknown>")
    assert_raises(ValueError, se_precomp.fit, S)


def test_spectral_embedding_unknown_affinity(seed=36):
    """Test that SpectralClustering fails with an unknown affinity type"""
    centers = np.array([
        [0., 0., 0.],
        [10., 10., 10.],
        [20., 20., 20.],
    ])
    X, true_labels = make_blobs(n_samples=100, centers=centers,
                                cluster_std=1., random_state=42)
    D = rbf_kernel(X)  # Distance matrix

    se_precomp = SpectralEmbedding(n_components=1, affinity="<unknown>",
                                   random_state=np.random.RandomState(seed))
    assert_raises(ValueError, se_precomp.fit, S)


def test_connectivity(seed=36):
    """Test that graph connectivity test works as expected"""
    graph = np.array([[0, 0, 0, 0, 0], 
                      [0, 0, 1, 0, 0],
                      [0, 1, 0, 1, 0],
                      [0, 0, 1, 0, 1],
                      [0, 0, 0, 1, 0]])
    assert_equal(_graph_is_connected(graph), False)
    graph = np.array([[0, 1, 0, 0, 0], 
                      [1, 0, 1, 0, 0],
                      [0, 1, 0, 1, 0],
                      [0, 0, 1, 0, 1],
                      [0, 0, 0, 1, 0]])
    assert_equal(_graph_is_connected(graph), True)
