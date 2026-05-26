import numpy as np
import pytest

from sklearn.cluster._fcm import FuzzyCMeans, _dist_to_membership
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score


def test_fcm_slide_numerical_example():
    """Verify FCM calculations match the uploaded slide numerical examples.

    Slides Data:
    Points: (1,3), (2,5), (4,8), (7,9)
    Initial membership matrix:
    Cluster 1: [0.8, 0.7, 0.2, 0.1]
    Cluster 2: [0.2, 0.3, 0.8, 0.9]
    """
    X = np.array([[1.0, 3.0], [2.0, 5.0], [4.0, 8.0], [7.0, 9.0]])
    U_init = np.array([[0.8, 0.2], [0.7, 0.3], [0.2, 0.8], [0.1, 0.9]])
    m = 2.0

    # 1. Update centroids using slide formulas
    Um = np.power(U_init, m)
    col_sums = np.sum(Um, axis=0, keepdims=True).T
    centers = (Um.T @ X) / col_sums

    # Centroid 1: (1.568, 4.051)
    # Centroid 2: (5.350, 8.215)
    np.testing.assert_allclose(centers[0], [1.568, 4.051], atol=2e-3)
    np.testing.assert_allclose(centers[1], [5.350, 8.215], atol=2e-3)

    # 2. Update membership matrix U
    from sklearn.metrics.pairwise import euclidean_distances

    D = euclidean_distances(X, centers)
    # Point 1: D_11 = 1.2, D_12 = 6.79 (rounded to 1 decimal place in slide,
    # we use full precision)
    # Let's verify slide distance values
    # D_11 = sqrt((1 - 1.568)^2 + (3 - 4.051)^2) = 1.1946 => 1.2
    # D_12 = sqrt((1 - 5.35)^2 + (3 - 8.215)^2) = 6.7865 => 6.79
    np.testing.assert_allclose(D[0, 0], 1.2, atol=1e-2)
    np.testing.assert_allclose(D[0, 1], 6.79, atol=1e-2)

    # Check slide new memberships for Point 1
    # Y_11 = ( (1.2^2 / 1.2^2)^1 + (1.2^2 / 6.79^2)^1 )^-1 = 0.97
    # Y_12 = 0.03
    # Let's verify that using _dist_to_membership gives exactly this value
    U_new = _dist_to_membership(D, m)
    np.testing.assert_allclose(U_new[0, 0], 0.97, atol=1e-2)
    np.testing.assert_allclose(U_new[0, 1], 0.03, atol=1e-2)


def test_fcm_input_validation():
    """Verify that parameters are validated and appropriate errors raised."""
    X = np.random.randn(10, 2)

    # n_samples < n_clusters
    fcm_est = FuzzyCMeans(n_clusters=12)
    with pytest.raises(ValueError, match="minimum of 12"):
        fcm_est.fit(X)

    # Invalid m
    with pytest.raises(ValueError):
        FuzzyCMeans(m=0.5).fit(X)

    # Invalid n_clusters
    with pytest.raises(ValueError):
        FuzzyCMeans(n_clusters=0).fit(X)


def test_fcm_initializations():
    """Verify that all initialization methods run successfully and behave correctly."""
    X, _ = make_blobs(n_samples=20, n_features=2, centers=3, random_state=42)

    # 1. Random membership (default)
    fcm = FuzzyCMeans(n_clusters=3, init="random_membership", random_state=42).fit(X)
    assert fcm.cluster_centers_.shape == (3, 2)
    assert fcm.u_.shape == (20, 3)

    # 2. Random centers
    fcm_centers = FuzzyCMeans(n_clusters=3, init="random_centers", random_state=42).fit(
        X
    )
    assert fcm_centers.cluster_centers_.shape == (3, 2)

    # 3. K-Means++
    fcm_km = FuzzyCMeans(n_clusters=3, init="k-means++", random_state=42).fit(X)
    assert fcm_km.cluster_centers_.shape == (3, 2)

    # 4. Custom initial centers array
    init_centers = X[:3].copy()
    fcm_custom = FuzzyCMeans(n_clusters=3, init=init_centers, random_state=42).fit(X)
    assert fcm_custom.cluster_centers_.shape == (3, 2)

    # 5. Custom callable
    def custom_init(X, n_clusters, random_state=None):
        return X[:n_clusters]

    fcm_callable = FuzzyCMeans(n_clusters=3, init=custom_init, random_state=42).fit(X)
    assert fcm_callable.cluster_centers_.shape == (3, 2)


def test_fcm_methods_and_attributes():
    """Verify fit, predict, predict_proba, transform, and score methods."""
    X, _ = make_blobs(n_samples=30, n_features=3, centers=2, random_state=0)

    fcm = FuzzyCMeans(n_clusters=2, random_state=0).fit(X)

    # Attributes
    assert hasattr(fcm, "cluster_centers_")
    assert hasattr(fcm, "u_")
    assert hasattr(fcm, "labels_")
    assert hasattr(fcm, "n_iter_")

    # Predict hard labels
    labels = fcm.predict(X)
    assert labels.shape == (30,)
    np.testing.assert_array_equal(labels, fcm.labels_)

    # Predict proba soft labels
    u_pred = fcm.predict_proba(X)
    assert u_pred.shape == (30, 2)
    np.testing.assert_allclose(np.sum(u_pred, axis=1), 1.0)
    np.testing.assert_allclose(u_pred, fcm.u_)

    # Transform (distances)
    distances = fcm.transform(X)
    assert distances.shape == (30, 2)
    # Check that minimum distance corresponds to predicted label
    assert np.all(np.argmin(distances, axis=1) == labels)

    # Score (FCM objective value)
    score = fcm.score(X)
    assert isinstance(score, float)
    assert score <= 0.0  # opposite of objective function


def test_fcm_convergence_to_kmeans():
    """Verify that when m is close to 1.0 (e.g. 1.05), FCM converges to KMeans.

    The FCM clustering results should converge to KMeans clustering results.
    """
    from sklearn.cluster import KMeans

    X, y = make_blobs(n_samples=50, n_features=2, centers=2, random_state=1)

    # Fit FuzzyCMeans with very small m (close to 1)
    fcm = FuzzyCMeans(n_clusters=2, m=1.05, random_state=42).fit(X)

    # Fit standard KMeans
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10).fit(X)

    # Compare clusterings (ARI should be very high, close or equal to 1.0)
    ari = adjusted_rand_score(kmeans.labels_, fcm.labels_)
    assert ari > 0.9  # Extremely high agreement


def test_fcm_overlapping_points_edge_case():
    """Verify that FCM handles duplicate data points and exact centers.

    Check that duplicate data points and exact cluster centers are handled gracefully.
    """
    # Data containing exact duplicate points
    X = np.array([[1.0, 2.0], [1.0, 2.0], [5.0, 5.0], [5.0, 5.0]])

    fcm = FuzzyCMeans(n_clusters=2, m=2.0, random_state=0).fit(X)

    # Check that predictions for duplicate points are well-defined
    u_pred = fcm.predict_proba(X)
    assert not np.any(np.isnan(u_pred))
    assert not np.any(np.isinf(u_pred))
    np.testing.assert_allclose(np.sum(u_pred, axis=1), 1.0)

    # Test single point exactly on a cluster center
    exact_point = np.array([fcm.cluster_centers_[0]])
    u_exact = fcm.predict_proba(exact_point)
    assert u_exact[0, 0] == 1.0
    assert u_exact[0, 1] == 0.0
