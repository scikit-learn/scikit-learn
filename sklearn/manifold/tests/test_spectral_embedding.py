from nose.tools import assert_true
from nose.tools import assert_equal

from scipy.sparse import csr_matrix
import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from sklearn.manifold.spectral_embedding import SpectralEmbedding
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.pipeline import Pipeline
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import KMeans, SpectralClustering

S = np.array([[1, 5, 2, 1, 0, 0, 0],
              [5, 1, 3, 1, 0, 0, 0],
              [2, 3, 1, 1, 0, 0, 0],
              [1, 1, 1, 1, 2, 1, 1],
              [0, 0, 0, 2, 2, 3, 2],
              [0, 0, 0, 1, 3, 1, 4],
              [0, 0, 0, 1, 2, 4, 1],
              ])


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
    """Test spectral embedding with knn graph"""
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
    assert_array_almost_equal(np.abs(embed_rbf), np.abs(embed_callable), 0)


def test_pipline_spectral_clustering():
    """Test using pipline to do spectral clustering"""
    SE = SpectralEmbedding()
    # dummy transform for SpectralEmbedding
    SE.transform = lambda x: x
    spectral_clustering = Pipeline([
        ('se', SE),
        ('km', KMeans()),
    ])

    for n_cluster in range(1, 5):
        n_cluster = 3
        spectral_clustering.set_params(km__n_clusters=n_cluster)
        spectral_clustering.set_params(se__n_components=n_cluster)
        spectral_clustering.set_params(se__gamma=1.0)
        spectral_clustering.fit(S)
        SC = SpectralClustering(n_clusters=n_cluster)
        SC.fit(S)
        assert_array_almost_equal(
            normalized_mutual_info_score(
                spectral_clustering.steps[1][1].labels_,
                SC.labels_), 0.0, 0)
