import numpy as np
import pytest

from scipy.sparse import csc_matrix
from sklearn.decomposition import graph_regularized_nmf
from sklearn.decomposition import _gdnmf as gdnmf
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from sklearn.utils._testing import assert_array_almost_equal


def test_convergence_warning():
    convergence_warning = (
        "Maximum number of iterations 1 reached. Increase it to improve convergence."
    )
    X = np.ones((2, 3))
    y = np.ones(2)
    nn = NearestNeighbors(n_neighbors=1)
    with pytest.warns(ConvergenceWarning, match=convergence_warning):
        graph_regularized_nmf(X, y, max_iter=1, knn=nn)


def test_initialize_w_h_a_output():
    # Test that initialization does not return negative values.
    random_state = np.random.RandomState(42)
    data = np.abs(random_state.randn(10, 10))
    W, H, A = gdnmf._init_w_h_a(data, 10, n_classes=3, random_state=random_state)
    assert not ((W < 0).any() or (H < 0).any() or (A < 0).any())


def test_initialize_w_h_a_stable():
    # Test that initialization returns close values for
    # repetitive runs with the same random state.
    random_state = np.random.RandomState(42)
    data = np.abs(random_state.randn(10, 10))
    W, H, A = gdnmf._init_w_h_a(data, 10, n_classes=3, random_state=0)
    assert not ((W < 0).any() or (H < 0).any() or (A < 0).any())
    W_2, H_2, A_2 = gdnmf._init_w_h_a(data, 10, n_classes=3, random_state=0)
    assert not ((W_2 < 0).any() or (H_2 < 0).any() or (A_2 < 0).any())
    assert_array_almost_equal(W_2, W, decimal=2)
    assert_array_almost_equal(H_2, H, decimal=2)
    assert_array_almost_equal(A_2, A, decimal=2)


def test_parameter_checking():
    X = np.ones((2, 3))
    y = np.ones(2)

    msg = "Negative values in data passed to"
    with pytest.raises(ValueError, match=msg):
        graph_regularized_nmf(-X, y)
    with pytest.raises(ValueError, match=msg):
        gdnmf._init_w_h_a(-X, 1, 2)
    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, n_components=-1)
    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, n_classes=-1)
    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, max_iter=-1)
    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, tol=-1)
    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, graph_coeff=-1)
    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, label_coeff=-1)

    class DummyNearestNeighbor:
        pass

    with pytest.raises(ValueError):
        graph_regularized_nmf(X, y, knn=DummyNearestNeighbor())


@pytest.mark.parametrize("graph_coeff", (0.0, 1.0))
@pytest.mark.parametrize("label_coeff", (0.0, 1.0))
def test_gdnmf_fit_nn_output(graph_coeff, label_coeff):
    # Test that the decomposition does not contain negative values.
    X = np.c_[5.0 - np.arange(1, 6), 5.0 + np.arange(1, 6)]
    y = np.arange(X.shape[0])
    W, H, _ = graph_regularized_nmf(
        X,
        y,
        n_components=2,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=0,
    )
    assert not ((W < 0).any() or (H < 0).any())


@pytest.mark.parametrize("graph_coeff", (0.0, 0.05, 0.1))
@pytest.mark.parametrize("label_coeff", (0.0, 0.05, 0.1))
def test_gdnmf_fit_close(graph_coeff, label_coeff):
    random_state = np.random.RandomState(42)
    # Test that the fit is not too far away from the original matrix.
    X = np.abs(random_state.randn(6, 5))
    y = random_state.randint(0, 3, X.shape[0])
    W, H, _ = graph_regularized_nmf(
        X,
        y,
        n_components=5,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        max_iter=1000,
    )
    assert np.linalg.norm(X - W @ H) < 0.05 * np.linalg.norm(X)


def test_n_components_greater_n_features():
    # Smoke test for the case of more components than features.
    random_state = np.random.RandomState(42)
    X = np.abs(random_state.randn(30, 10))
    y = random_state.randint(0, 5, X.shape[0])
    graph_regularized_nmf(
        X, y, n_components=15, n_classes=5, random_state=random_state, tol=1e-2
    )


@pytest.mark.parametrize("graph_coeff", (0.0, 0.5, 1.0))
@pytest.mark.parametrize("label_coeff", (0.0, 0.5, 1.0))
def test_gdnmf_sparse_input(graph_coeff, label_coeff):
    # Test that sparse matrices are accepted as input.
    random_state = np.random.RandomState(42)
    X = np.abs(random_state.randn(10, 10))
    y = random_state.randint(0, 3, X.shape[0])
    X[:, 2 * np.arange(5)] = 0
    X_sparse = csc_matrix(X)

    W_1, H_1, _ = graph_regularized_nmf(
        X,
        y,
        n_components=5,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=0,
        tol=1e-2,
    )
    W_2, H_2, _ = graph_regularized_nmf(
        X_sparse,
        y,
        n_components=5,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=0,
        tol=1e-2,
    )

    assert_array_almost_equal(W_1, W_2)
    assert_array_almost_equal(H_1, H_2)


@pytest.mark.parametrize("graph_coeff", (0.0, 0.5, 1.0))
@pytest.mark.parametrize("label_coeff", (10, 20))
def test_same_label_samples_closer(graph_coeff, label_coeff):
    # Test the transformed samples with the same labels
    # will be closer together.
    random_state = np.random.RandomState(42)
    N_CLASSES = 3
    X = np.abs(random_state.randn(10, 10))
    y = random_state.randint(0, N_CLASSES, X.shape[0])

    W, *_ = graph_regularized_nmf(
        X,
        y,
        n_components=4,
        n_classes=N_CLASSES,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        tol=1e-5,
    )

    def intra_group_avg_dist(A):
        n = A.shape[0]
        return np.triu(pairwise_distances(A)).sum() / ((n - 1) * n / 2)

    avg_pair_dist = intra_group_avg_dist(W)
    for cls in range(N_CLASSES):
        cls_samples = W[y == cls, :]
        cls_pair_dist = intra_group_avg_dist(cls_samples)
        assert cls_pair_dist < avg_pair_dist


@pytest.mark.parametrize("graph_coeff", (50, 100))
@pytest.mark.parametrize("label_coeff", (0.0, 0.5, 1.0))
def test_preservation_of_knn(graph_coeff, label_coeff):
    # Test nearest-neighbour structure is preserved
    # during the tranformation.
    random_state = np.random.RandomState(42)
    N_CLASSES = 3
    N_NEIGHBORS = 3
    X = np.abs(random_state.randn(10, 10))
    y = random_state.randint(0, N_CLASSES, X.shape[0])
    knn = NearestNeighbors(n_neighbors=N_NEIGHBORS)
    knn.fit(X)
    C_1 = knn.kneighbors_graph(X, mode="connectivity").toarray()

    W, *_ = graph_regularized_nmf(
        X,
        y,
        n_components=4,
        n_classes=N_CLASSES,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        tol=1e-5,
        knn=NearestNeighbors(n_neighbors=N_NEIGHBORS),
    )

    knn = NearestNeighbors(n_neighbors=N_NEIGHBORS)
    knn.fit(W)
    C_2 = knn.kneighbors_graph(W, mode="connectivity").toarray()

    assert np.sum(np.abs(C_1 - C_2)) <= 10
