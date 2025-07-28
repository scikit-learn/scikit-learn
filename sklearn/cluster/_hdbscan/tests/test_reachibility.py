import numpy as np
import pytest

from sklearn.cluster._hdbscan._reachability import mutual_reachability_graph
from sklearn.utils._testing import assert_allclose
from sklearn.utils.fixes import CSC_CONTAINERS, CSR_CONTAINERS


@pytest.mark.parametrize("csc_container", CSC_CONTAINERS)
def test_mutual_reachability_graph_error_sparse_format(csc_container):
    """Check that we raise an error if the sparse format is not CSR."""
    rng = np.random.RandomState(0)
    X = rng.randn(10, 10)
    X = X.T @ X
    np.fill_diagonal(X, 0.0)
    X = csc_container(X)

    err_msg = "Only sparse CSR matrices are supported"
    with pytest.raises(ValueError, match=err_msg):
        mutual_reachability_graph(X)


@pytest.mark.parametrize("csr_container", [np.asarray] + CSR_CONTAINERS)
def test_mutual_reachability_graph_inplace(csr_container):
    """Check that the operation is happening inplace."""
    rng = np.random.RandomState(0)
    X = rng.randn(10, 10)
    X = X.T @ X
    np.fill_diagonal(X, 0.0)
    X = csr_container(X)

    mr_graph = mutual_reachability_graph(X)

    assert id(mr_graph) == id(X)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_mutual_reachability_graph_equivalence_dense_sparse(csr_container):
    """Check that we get the same results for dense and sparse implementation."""
    rng = np.random.RandomState(0)
    X = rng.randn(5, 5)
    X_dense = X.T @ X
    X_sparse = csr_container(X_dense)

    mr_graph_dense = mutual_reachability_graph(X_dense, min_samples=3)
    mr_graph_sparse = mutual_reachability_graph(X_sparse, min_samples=3)

    assert_allclose(mr_graph_dense, mr_graph_sparse.toarray())


@pytest.mark.parametrize("csr_container", [np.asarray] + CSR_CONTAINERS)
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_mutual_reachability_graph_preserves_dtype(csr_container, dtype):
    """Check that the computation preserve dtype thanks to fused types."""
    rng = np.random.RandomState(0)
    X = rng.randn(10, 10)
    X = (X.T @ X).astype(dtype)
    np.fill_diagonal(X, 0.0)
    X = csr_container(X, dtype=dtype)

    assert X.dtype == dtype
    mr_graph = mutual_reachability_graph(X)
    assert mr_graph.dtype == dtype
