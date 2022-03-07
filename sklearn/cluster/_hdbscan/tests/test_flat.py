"""
Simple tests for flat clustering over HDBSCAN hierarchy
"""
import warnings
import numpy as np

from sklearn.cluster import (
    HDBSCAN,
    approximate_predict,
    HDBSCAN_flat,
    approximate_predict_flat,
    membership_vector_flat,
    all_points_membership_vectors_flat,
)

from sklearn.datasets import make_blobs, make_moons
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.utils._testing import assert_array_equal, assert_array_less

# Ignore future warnings thrown by sklearn
warnings.filterwarnings("ignore", category=FutureWarning)

# Create a nice dataset with 6 circular clusters and 2 moons
centers = [(0, 2), (-0.2, 0), (0.2, 0), (1.5, 0), (2.0, 1.0), (2.5, 0.0)]
std = [0.5, 0.08, 0.06, 0.35, 0.35, 0.35]
X0, y0 = make_blobs(
    n_samples=[70, 30, 80, 100, 40, 150],
    centers=centers,
    cluster_std=std,
    random_state=1,
)
X1, y1 = make_moons(n_samples=300, noise=0.07, random_state=42)
X1 += 3.0
y1 += len(centers)
X = np.vstack((X0, X1))
y = np.concatenate((y0, y1))

X, X_test, y, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
scaler = StandardScaler()
X = scaler.fit_transform(X)
X_test = scaler.transform(X_test)


def n_clusters_from_labels(labels_):
    return np.amax(labels_) + 1


def test_flat_base_default():
    """
    Verify that the default clustering of HDBSCAN is preserved.
    """
    # Given, the base HDBSCAN with method 'eom'
    clusterer = HDBSCAN(cluster_selection_method="eom").fit(X)
    n_clusters = n_clusters_from_labels(clusterer.labels_)

    # When we ask for flat clustering with same n_clusters,
    clusterer_flat = HDBSCAN_flat(
        X, n_clusters=n_clusters, cluster_selection_method="eom"
    )

    # Then, the labels and probabilities should match
    assert_array_equal(clusterer_flat.labels_, clusterer.labels_)
    assert_array_equal(clusterer_flat.probabilities_, clusterer.probabilities_)

    # Given, the base HDBSCAN with method 'leaf'
    clusterer = HDBSCAN(cluster_selection_method="leaf").fit(X)
    n_clusters = n_clusters_from_labels(clusterer.labels_)

    # When we ask for flat clustering with same n_clusters,
    clusterer_flat = HDBSCAN_flat(
        X, n_clusters=n_clusters, cluster_selection_method="leaf"
    )

    # Then, the labels and probabilities should match
    assert_array_equal(clusterer_flat.labels_, clusterer.labels_)
    assert_array_equal(clusterer_flat.probabilities_, clusterer.probabilities_)
    return


def test_flat_base_epsilon():
    """
    Verify that a clustering of HDBSCAN specified by
        cluster_selection_epsilon is preserved.
    """
    # Method 'eom'...
    # Given, a flat clustering for required n_clusters,
    n_clusters = 4
    clusterer_flat = HDBSCAN_flat(
        X, n_clusters=n_clusters, cluster_selection_method="eom"
    )

    # When we run the base HDBSCAN using it's epsilon,
    epsilon = clusterer_flat.cluster_selection_epsilon
    clusterer = HDBSCAN(
        cluster_selection_method="eom", cluster_selection_epsilon=epsilon
    ).fit(X)

    # Then, the labels and probabilities should match
    assert_array_equal(clusterer_flat.labels_, clusterer.labels_)
    assert_array_equal(clusterer_flat.probabilities_, clusterer.probabilities_)

    # Method 'leaf'...
    # Given, a flat clustering for required n_clusters,
    n_clusters = 6
    clusterer_flat = HDBSCAN_flat(
        X, n_clusters=n_clusters, cluster_selection_method="leaf"
    )

    # When we run the base HDBSCAN using it's epsilon,
    epsilon = clusterer_flat.cluster_selection_epsilon
    clusterer = HDBSCAN(
        cluster_selection_method="leaf", cluster_selection_epsilon=epsilon
    ).fit(X)

    # Then, the labels and probabilities should match
    assert_array_equal(clusterer_flat.labels_, clusterer.labels_)
    assert_array_equal(clusterer_flat.probabilities_, clusterer.probabilities_)
    return


def test_switch_to_leaf():
    """
    Verify that when we request more clusters than 'eom' can handle,
        method switches to 'leaf' and the results match 'leaf'.
    """
    # Given the max number of clusters that can be produced by 'eom',
    #   (these are produced for epsilon=0) (??? Needs verification)
    clusterer = HDBSCAN(
        cluster_selection_method="eom", cluster_selection_epsilon=0
    ).fit(X)
    max_clusters = n_clusters_from_labels(clusterer.labels_)

    with warnings.catch_warnings(record=True) as w:
        # When we try flat clustering with 'eom' method for more n_clusters,
        clusterer_flat = HDBSCAN_flat(
            X, cluster_selection_method="eom", n_clusters=max_clusters + 2
        )
        # Then, a warning is raised saying 'eom' can't get this clustering,
        assert len(w) > 0
        assert issubclass(w[-1].category, UserWarning)
        assert "Cannot predict" in str(w[-1].message)

    # the resulting clusterer switches to using method 'leaf',
    assert (
        clusterer_flat.cluster_selection_method == "leaf"
    ), "cluster selection method has not switched to 'leaf'"
    # and the resulting probabilities and labels must match
    epsilon = clusterer_flat.cluster_selection_epsilon
    clusterer_leaf = HDBSCAN(
        cluster_selection_method="leaf", cluster_selection_epsilon=epsilon
    ).fit(X)
    assert_array_equal(clusterer_flat.labels_, clusterer_leaf.labels_)
    assert_array_equal(clusterer_flat.probabilities_, clusterer_leaf.probabilities_)
    return


def test_approx_predict_default():
    """
    Verify that approximate_predict_flat produces same results as default
    """
    # Given the base HDBSCAN trained on some data,
    clusterer = HDBSCAN(
        cluster_selection_method="eom",
        cluster_selection_epsilon=0,
        prediction_data=True,
    ).fit(X)

    # When using approximate_predict_flat without specifying n_clusters,
    labels_flat, proba_flat = approximate_predict_flat(
        clusterer, X_test, n_clusters=None
    )

    # Then, the clustering should match that due to approximate_predict,
    labels_base, proba_base = approximate_predict(clusterer, X_test)
    assert_array_equal(labels_flat, labels_base)
    assert_array_equal(proba_flat, proba_base)
    return


def test_approx_predict_same_clusters():
    """
    Verify that approximate_predict_flat produces as many clusters as clusterer
    """
    # Given a flat clustering trained for some n_clusters,
    n_clusters = 5
    clusterer = HDBSCAN_flat(X, cluster_selection_method="eom", n_clusters=n_clusters)

    # When using approximate_predict_flat without specifying n_clusters,
    labels_flat, proba_flat = approximate_predict_flat(
        clusterer, X_test, n_clusters=None
    )

    # Then, the number of clusters produced must match the original n_clusters
    n_clusters_out = n_clusters_from_labels(labels_flat)
    assert n_clusters_out == n_clusters
    # and all probabilities are <= 1.
    assert_array_less(proba_flat, np.ones(len(proba_flat)) + 1.0e-14)
    return


def test_approx_predict_diff_clusters():
    """
    Verify that approximate_predict_flat produces as many clusters as asked
    """
    # Given a flat clustering trained for some n_clusters,
    n_clusters_fit = 5
    clusterer = HDBSCAN_flat(
        X,
        cluster_selection_method="eom",
        n_clusters=n_clusters_fit,
        prediction_data=True,
    )

    # When using approximate_predict_flat with specified n_clusters,
    n_clusters_predict = 3
    labels_flat, proba_flat = approximate_predict_flat(
        clusterer, X_test, n_clusters=n_clusters_predict
    )

    # Then, the requested number of clusters must be produced
    n_clusters_out = n_clusters_from_labels(labels_flat)
    assert n_clusters_out == n_clusters_predict
    # and all probabilities are <= 1.
    assert_array_less(proba_flat, np.ones(len(proba_flat)) + 1.0e-14)

    # When using approximate_predict_flat with more clusters
    #   than 'eom' can handle,
    n_clusters_predict = 12
    with warnings.catch_warnings(record=True) as w:
        labels_flat, proba_flat = approximate_predict_flat(
            clusterer, X_test, n_clusters=n_clusters_predict
        )
        # Then, a warning is raised saying 'eom' can't get this clustering,
        assert len(w) > 0
        assert issubclass(w[-1].category, UserWarning)
        assert "Cannot predict" in str(w[-1].message)
    # But the requested number of clusters must still be produced using 'leaf'
    n_clusters_out = n_clusters_from_labels(labels_flat)
    assert n_clusters_out == n_clusters_predict
    # and all probabilities are <= 1.
    assert_array_less(proba_flat, np.ones(len(proba_flat)) + 1.0e-14)
    return


def test_mem_vec_same_clusters():
    """
    Verify membership vector produces same n_clusters as clusterer
    """
    # Given a flat clustering trained for n_clusters picked by HDBSCAN,
    n_clusters_fit = None
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)

    # When membership_vector_flat is called with new data,
    memberships = membership_vector_flat(clusterer, X_test)

    # Then the number of clusters in memberships matches those of clusterer,
    assert memberships.shape[1] == n_clusters_from_labels(clusterer.labels_)
    # and the number of points should equal those in the test set
    assert len(memberships) == len(X_test)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)

    # ========================================
    # Given a flat clustering for a specified n_clusters,
    n_clusters_fit = n_clusters_from_labels(clusterer.labels_) - 2
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)

    # When membership_vector_flat is called with new data,
    memberships = membership_vector_flat(clusterer, X_test)

    # Then the number of clusters in memberships matches those of clusterer,
    assert memberships.shape[1] == n_clusters_fit
    # and the number of points should equal those in the test set
    assert len(memberships) == len(X_test)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)
    return


def test_mem_vec_diff_clusters():
    """
    Verify membership vector produces as many clusters as requested
    """
    # Ignore user warnings in this function
    warnings.filterwarnings("ignore", category=UserWarning)

    # Given a flat clustering trained for n_clusters picked by HDBSCAN,
    n_clusters_fit = None
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)
    n_clusters_fitted = n_clusters_from_labels(clusterer.labels_)

    # When membership_vector_flat is called with new data for some n_clusters,
    n_clusters_predict = n_clusters_fitted + 3
    memberships = membership_vector_flat(
        clusterer, X_test, n_clusters=n_clusters_predict
    )

    # Then the number of clusters in memberships should be as requested,
    assert memberships.shape[1] == n_clusters_predict
    # and the number of points should equal those in the test set
    assert len(memberships) == len(X_test)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)

    # ========================================
    # Given a flat clustering for a specified n_clusters,
    n_clusters_fit = n_clusters_from_labels(clusterer.labels_) + 2
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)

    # When membership_vector_flat is called with new data for some n_clusters,
    n_clusters_predict = n_clusters_fit + 3
    memberships = membership_vector_flat(
        clusterer, X_test, n_clusters=n_clusters_predict
    )

    # Then the number of clusters in memberships should be as requested,
    assert memberships.shape[1] == n_clusters_predict
    # and the number of points should equal those in the test set
    assert len(memberships) == len(X_test)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)
    return


def test_all_points_mem_vec_same_clusters():
    """
    Verify membership vector for training set produces same n_clusters
        as clusterer
    """
    # Given a flat clustering trained for n_clusters picked by HDBSCAN,
    n_clusters_fit = None
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)

    # When all_points_membership_vectors_flat is called,
    memberships = all_points_membership_vectors_flat(clusterer)

    # Then the number of clusters in memberships matches those of clusterer,
    assert memberships.shape[1] == n_clusters_from_labels(clusterer.labels_)
    # and the number of points should equal those in the training set
    assert len(memberships) == len(X)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)

    # ========================================
    # Given a flat clustering for a specified n_clusters,
    n_clusters_fit = n_clusters_from_labels(clusterer.labels_) - 2
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)

    # When all_points_membership_vectors_flat is called,
    memberships = all_points_membership_vectors_flat(clusterer)

    # Then the number of clusters in memberships matches those of clusterer,
    assert memberships.shape[1] == n_clusters_from_labels(clusterer.labels_)
    # and the number of points should equal those in the training set
    assert len(memberships) == len(X)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)
    return


def test_all_points_mem_vec_diff_clusters():
    """
    Verify membership vector produces as many clusters as requested
    """
    # Ignore user warnings in this function
    warnings.filterwarnings("ignore", category=UserWarning)

    # Given a flat clustering trained for n_clusters picked by HDBSCAN,
    n_clusters_fit = None
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)
    n_clusters_fitted = n_clusters_from_labels(clusterer.labels_)

    # When all_points_membership_vectors_flat is called for some n_clusters,
    n_clusters_predict = n_clusters_fitted + 3
    memberships = all_points_membership_vectors_flat(
        clusterer, n_clusters=n_clusters_predict
    )

    # Then the number of clusters in memberships should be as requested,
    assert memberships.shape[1] == n_clusters_predict
    # and the number of points should equal those in the training set
    assert len(memberships) == len(X)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)

    # ========================================
    # Given a flat clustering for a specified n_clusters,
    n_clusters_fit = n_clusters_from_labels(clusterer.labels_) + 2
    clusterer = HDBSCAN_flat(X, n_clusters=n_clusters_fit)

    # When membership_vector_flat is called for some n_clusters,
    n_clusters_predict = n_clusters_fitted + 3
    memberships = all_points_membership_vectors_flat(
        clusterer, n_clusters=n_clusters_predict
    )

    # Then the number of clusters in memberships should be as requested,
    assert memberships.shape[1] == n_clusters_predict
    # and the number of points should equal those in the training set
    assert len(memberships) == len(X)
    # and all probabilities are <= 1.
    assert_array_less(memberships, np.ones(memberships.shape) + 1.0e-14)
    return
