import numpy as np
import pytest

from scipy.sparse import csc_matrix
from sklearn.base import clone
from sklearn.decomposition import GDNMF
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
        GDNMF(max_iter=1, knn=nn).fit(X, y)


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
    nn = NearestNeighbors(n_neighbors=1)

    msg = "Negative values in data passed to"
    with pytest.raises(ValueError, match=msg):
        GDNMF(knn=nn).fit(-X, y)
    with pytest.raises(ValueError, match=msg):
        gdnmf._init_w_h_a(-X, 1, 2)

    clf = GDNMF(1, tol=0.1, knn=nn).fit(X, y)
    with pytest.raises(ValueError, match=msg):
        clf.transform(-X)


@pytest.mark.parametrize("graph_coeff", (0.0, 1.0))
@pytest.mark.parametrize("label_coeff", (0.0, 1.0))
def test_gdnmf_fit_nn_output(graph_coeff, label_coeff):
    # Test that the decomposition does not contain negative values.
    X = np.c_[5.0 - np.arange(1, 6), 5.0 + np.arange(1, 6)]
    y = np.arange(X.shape[0])
    model = GDNMF(
        n_components=2,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=0,
    )
    transf = model.fit_transform(X, y)
    assert not ((model.components_ < 0).any() or (transf < 0).any())


@pytest.mark.parametrize("graph_coeff", (0.0, 0.05, 0.1))
@pytest.mark.parametrize("label_coeff", (0.0, 0.05, 0.1))
def test_gdnmf_fit_close(graph_coeff, label_coeff):
    random_state = np.random.RandomState(42)
    # Test that the fit is not too far away from the original matrix.
    gdnmf = GDNMF(
        n_components=5,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        max_iter=1000,
    )
    X = np.abs(random_state.randn(6, 5))
    y = random_state.randint(0, 3, X.shape[0])
    assert gdnmf.fit(X, y).reconstruction_err_ < 0.05 * np.linalg.norm(X)


@pytest.mark.parametrize("graph_coeff", (0.0, 0.05, 0.1))
@pytest.mark.parametrize("label_coeff", (0.0, 0.05, 0.1))
def test_gdnmf_transform(graph_coeff, label_coeff):
    # Test that GDNMF.transform returns a value close to fit_transform.
    random_state = np.random.RandomState(42)
    X = np.abs(random_state.randn(6, 5))
    y = random_state.randint(0, 3, X.shape[0])
    nmf = GDNMF(
        n_components=3,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=0,
        tol=1e-5,
    )
    W_1 = nmf.fit_transform(X, y)
    W_2 = nmf.transform(X)
    assert_array_almost_equal(W_1, W_2)


@pytest.mark.parametrize("graph_coeff", (0.0, 0.05, 0.1))
@pytest.mark.parametrize("label_coeff", (0.0, 0.05, 0.1))
def test_gdnmf_inverse_transform(graph_coeff, label_coeff):
    # Test that GDNMF.inverse_transform returns a matrix close
    # to the original matrix.
    random_state = np.random.RandomState(0)
    X = np.abs(random_state.randn(6, 4))
    y = random_state.randint(0, 3, X.shape[0])
    m = GDNMF(
        n_components=4,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        max_iter=1000,
    )
    ft = m.fit_transform(X, y)
    X_new = m.inverse_transform(ft)
    assert np.linalg.norm(X_new - X) < 0.05 * np.linalg.norm(X)


def test_n_components_greater_n_features():
    # Smoke test for the case of more components than features.
    random_state = np.random.RandomState(42)
    X = np.abs(random_state.randn(30, 10))
    y = random_state.randint(0, 5, X.shape[0])
    GDNMF(n_components=15, n_classes=5, random_state=random_state, tol=1e-2).fit(X, y)


@pytest.mark.parametrize("graph_coeff", (0.0, 0.5, 1.0))
@pytest.mark.parametrize("label_coeff", (0.0, 0.5, 1.0))
def test_gdnmf_sparse_input(graph_coeff, label_coeff):
    # Test that sparse matrices are accepted as input.
    random_state = np.random.RandomState(42)
    X = np.abs(random_state.randn(10, 10))
    y = random_state.randint(0, 3, X.shape[0])
    X[:, 2 * np.arange(5)] = 0
    X_sparse = csc_matrix(X)

    est1 = GDNMF(
        n_components=5,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        tol=1e-2,
    )
    est2 = clone(est1)

    W_1 = est1.fit_transform(X, y)
    W_2 = est2.fit_transform(X_sparse, y)
    H_1 = est1.components_
    H_2 = est2.components_

    assert_array_almost_equal(W_1, W_2)
    assert_array_almost_equal(H_1, H_2)


@pytest.mark.parametrize("graph_coeff", (0.0, 0.01, 0.02))
@pytest.mark.parametrize("label_coeff", (0.0, 0.05, 0.1))
def test_gdnmf_sparse_transform(graph_coeff, label_coeff):
    # Test that transform works on sparse data.
    random_state = np.random.RandomState(42)
    X = np.abs(random_state.randn(10, 10))
    y = random_state.randint(0, 3, X.shape[0])
    X[:, 2 * np.arange(5)] = 0
    X = csc_matrix(X)

    model = GDNMF(
        n_components=3,
        n_classes=3,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=0,
        tol=1e-5,
    )
    W_1 = model.fit_transform(X, y)
    W_2 = model.transform(X)
    assert_array_almost_equal(W_1, W_2)


@pytest.mark.parametrize("graph_coeff", (0.0, 0.5, 1.0))
@pytest.mark.parametrize("label_coeff", (10, 20))
def test_same_label_samples_closer(graph_coeff, label_coeff):
    # Test the transformed samples with the same labels
    # will be closer together.
    random_state = np.random.RandomState(42)
    N_CLASSES = 3
    X = np.abs(random_state.randn(10, 10))
    y = random_state.randint(0, N_CLASSES, X.shape[0])
    model = GDNMF(
        n_components=4,
        n_classes=N_CLASSES,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        tol=1e-5,
    )

    W = model.fit_transform(X, y)

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
    model = GDNMF(
        n_components=4,
        n_classes=N_CLASSES,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
        random_state=random_state,
        tol=1e-5,
        knn=NearestNeighbors(n_neighbors=N_NEIGHBORS),
    )

    W = model.fit_transform(X, y)

    knn = NearestNeighbors(n_neighbors=N_NEIGHBORS)
    knn.fit(W)
    C_2 = knn.kneighbors_graph(W, mode="connectivity").toarray()

    assert np.sum(np.abs(C_1 - C_2)) <= 10
