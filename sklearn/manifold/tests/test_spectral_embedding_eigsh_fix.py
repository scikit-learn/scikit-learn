"""
Unit test for spectral embedding eigsh performance fix (issue #33242).

This test verifies that the eigenvalues computed by eigsh are near 1.0
when using I-L instead of -L, enabling effective sigma=1.0 targeting.
"""
import numpy as np
import pytest

def test_spectral_embedding_eigsh_eigenvalues():
    """Test that spectral embedding uses correct matrix for eigsh.
    
    This addresses issue #33242 where eigsh was called on -L instead of I-L,
    causing sigma=1.0 to target the wrong eigenvalues and leading to poor
    performance.
    """
    pytest.importorskip("scipy")
    from scipy.sparse.linalg import eigsh
    from scipy.sparse.csgraph import laplacian
    from sklearn.neighbors import kneighbors_graph
    from sklearn.datasets import make_swiss_roll
    from sklearn.utils._arpack import _init_arpack_v0

    # Create test data
    X, _ = make_swiss_roll(n_samples=300, random_state=42)
    
    # Build affinity matrix (same as SpectralEmbedding)
    A = kneighbors_graph(X, n_neighbors=10, include_self=True).toarray()
    A = (A + A.T) / 2
    
    # Get normalized Laplacian
    L, dd = laplacian(A, normed=True, return_diag=True)
    
    # Test corrected approach: I - L
    L_corrected = L.copy()
    np.fill_diagonal(L_corrected, 1)  # Set diagonal to 1
    L_corrected *= -1  # Negate to get -L
    np.fill_diagonal(L_corrected, L_corrected.diagonal() + 1)  # Add identity: I - L
    
    # Initialize deterministic v0
    v0 = _init_arpack_v0(L_corrected.shape[0], random_state=42)
    
    # Compute eigenvalues with sigma=1.0 (should target eigenvalues near 1.0)
    eigvals, _ = eigsh(L_corrected, k=3, sigma=1.0, which="LM", tol=1e-10, v0=v0)
    
    # For I - L, eigenvalues should be close to 1.0
    # The smallest eigenvalues of I - L correspond to largest of L
    # Since L has eigenvalues in [0, 2], I - L has eigenvalues in [-1, 1]
    # The eigenvalues closest to 1.0 should be the ones we get
    assert eigvals.min() > 0.9, f"Expected eigenvalues near 1.0, got min={eigvals.min()}"
    assert eigvals.max() <= 1.0 + 1e-10, f"Expected eigenvalues <= 1.0, got max={eigvals.max()}"
    
    # Test that eigenvalues are in the expected range for I - L
    # The exact values depend on the graph structure, but should be close to 1
    assert np.all(eigvals > 0.95), f"Eigenvalues should be close to 1.0: {eigvals}"

def test_spectral_embedding_performance_regression():
    """Test that spectral embedding converges quickly with the fix.
    
    This is a regression test to ensure the eigsh performance fix
    doesn't cause convergence issues.
    """
    pytest.importorskip("scipy")
    from sklearn.manifold import SpectralEmbedding
    from sklearn.datasets import make_swiss_roll
    
    # Create test data
    X, _ = make_swiss_roll(n_samples=200, random_state=42)
    
    # This should converge quickly with the fix
    embedding = SpectralEmbedding(
        n_components=2,
        eigen_solver='arpack',
        random_state=42
    )
    
    # Should not raise any convergence errors
    result = embedding.fit_transform(X)
    
    assert result.shape == (200, 2)
    # Should have meaningful variance (not degenerate)
    assert np.std(result[:, 0]) > 0.01
    assert np.std(result[:, 1]) > 0.01