import numpy as np
from sklearn.cluster import SphericalKMeans

def test_spherical_kmeans_basic():
    rng = np.random.RandomState(42)
    X = rng.randn(100, 10)
    X = X / np.linalg.norm(X, axis=1, keepdims=True)  # Make unit vectors

    model = SphericalKMeans(n_clusters=3, random_state=42)
    model.fit(X)
    assert model.labels_.shape == (100,)
    assert model.cluster_centers_.shape == (3, 10)
    # Centroids should also be unit norm
    norms = np.linalg.norm(model.cluster_centers_, axis=1)
    np.testing.assert_allclose(norms, 1, atol=1e-7)

def test_predict():
    X = np.eye(5)
    model = SphericalKMeans(n_clusters=2, random_state=0)
    model.fit(X)
    labels = model.predict(X)
    assert set(labels) <= {0, 1}
