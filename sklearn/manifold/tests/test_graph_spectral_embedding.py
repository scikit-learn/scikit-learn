import pytest
import numpy as np

from scipy.sparse import csr_matrix

from sklearn.cluster import KMeans
from sklearn.manifold import GraphSpectralEmbedding
from sklearn.metrics import adjusted_rand_score
from sklearn.datasets import make_erdos_reyni_graph, make_sbm_graph


def _kmeans_comparison(data, labels, n_clusters):
    """
    Function for comparing the ARIs of kmeans clustering for arbitrary number of
    data/labels

    Parameters
    ----------
        data: list-like
            each element in the list is a dataset to perform k-means on
        labels: list-like
            each element in the list is a set of lables with the same number of points
            as the corresponding data
        n_clusters: int
            the number of clusters to use for k-means

    Returns
    -------
        aris: list, length the same as data/labels
            the i-th element in the list is an ARI (Adjusted Rand Index) corresponding
            to the result of k-means clustering on the i-th data/labels
    """

    if len(data) != len(labels):
        raise ValueError("Must have same number of labels and data")

    aris = []
    for i in range(0, len(data)):
        kmeans_prediction = KMeans(n_clusters=n_clusters).fit_predict(data[i])
        aris.append(adjusted_rand_score(labels[i], kmeans_prediction))

    return aris


def _test_inputs(X, error_type, **kws):
    with pytest.raises(error_type):
        gse = GraphSpectralEmbedding(**kws)
        gse.fit(X)


@pytest.mark.parametrize("method", ["LSE", "ASE"])
def test_output_dim_directed(method):
    n_components = 4
    embed = GraphSpectralEmbedding(
        n_components=n_components, concat=True, algorithm=method
    )
    n = 10
    M = 20
    A = make_erdos_reyni_graph(n, M, directed=True) + 5
    assert embed.fit_transform(A).shape == (n, 8)
    assert embed.latent_left_.shape == (n, 4)
    assert embed.latent_right_.shape == (n, 4)


@pytest.mark.parametrize("method", ["LSE", "ASE"])
def test_output_dim(method, sparse=False, *args, **kwargs):
    n_components = 4
    embed = GraphSpectralEmbedding(n_components=n_components, algorithm=method)
    n = 10
    M = 20
    A = make_erdos_reyni_graph(n, M) + 5
    if sparse:
        A = csr_matrix(A)
    embed.fit(A)
    assert embed.latent_left_.shape == (n, 4)
    assert embed.latent_right_ is None


@pytest.mark.parametrize("method", ["LSE", "ASE"])
@pytest.mark.parametrize("P", [np.array([[0.8, 0.2], [0.2, 0.8]])])
@pytest.mark.parametrize("directed", [True, False])
def test_sbm_er_binary(method, P, directed, sparse=True, *args, **kwargs):
    np.random.seed(8888)

    num_sims = 50
    verts = 200
    communities = 2

    verts_per_community = [100, 100]

    sbm_wins = 0
    er_wins = 0
    for sim in range(0, num_sims):
        sbm_sample = make_sbm_graph(verts_per_community, P, directed=directed)
        er = make_sbm_graph(np.asarray([verts]), np.asarray([[0.5]]), directed=directed)
        if sparse:
            sbm_sample = csr_matrix(sbm_sample)
            er = csr_matrix(er)
        embed_sbm = GraphSpectralEmbedding(
            n_components=2, algorithm=method, concat=directed
        )
        embed_er = GraphSpectralEmbedding(
            n_components=2, algorithm=method, concat=directed
        )

        labels_sbm = np.zeros((verts), dtype=np.int8)
        labels_er = np.zeros((verts), dtype=np.int8)
        labels_sbm[100:] = 1
        labels_er[100:] = 1

        X_sbm = embed_sbm.fit_transform(sbm_sample)
        X_er = embed_er.fit_transform(er)

        if directed:
            assert X_sbm.shape == (verts, 2 * communities)
            assert X_er.shape == (verts, 2 * communities)
        else:
            assert X_sbm.shape == (verts, communities)
            assert X_er.shape == (verts, communities)

        aris = _kmeans_comparison((X_sbm, X_er), (labels_sbm, labels_er), communities)
        sbm_wins = sbm_wins + (aris[0] > aris[1])
        er_wins = er_wins + (aris[0] < aris[1])

    assert sbm_wins > er_wins


def test_input_params():
    X = make_erdos_reyni_graph(10, 20)

    # value error, check n_components type int
    _test_inputs(X, ValueError, n_components="not_int")

    # n_components must be <= min(X.shape)
    _test_inputs(X, ValueError, n_components=15)

    # value error, check n_elbow type int
    _test_inputs(X, ValueError, n_elbows="not_int")

    # value error, check n_elbow > 1
    _test_inputs(X, ValueError, n_elbows=0)

    # check algorithm string
    _test_inputs(X, ValueError, algorithm="wrong")
    _test_inputs(X, TypeError, algorithm=1)

    # svd_solver string
    _test_inputs(X, ValueError, svd_solver="wrong")

    # regularizer not int, float, or bool
    _test_inputs(X, TypeError, algorithm="lse", regularizer="wrong")
    # regularizer not greater than 0
    _test_inputs(X, ValueError, algorithm="lse", regularizer=-1)


def test_unconnected_warning():
    A = make_erdos_reyni_graph(100, 10)
    with pytest.warns(UserWarning):
        ase = GraphSpectralEmbedding(algorithm="ASE")
        ase.fit(A)
    A = make_erdos_reyni_graph(100, 10, directed=True)
    with pytest.warns(UserWarning):
        ase = GraphSpectralEmbedding(algorithm="ASE")
        ase.fit(A)
