"""Testing for K-means"""
import sys

import numpy as np
from scipy import sparse as sp

from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils._testing import assert_almost_equal

from sklearn.utils.extmath import row_norms
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster._kmeans import _labels_inertia
from sklearn.cluster._kmeans import _mini_batch_step
from sklearn.datasets import make_blobs
from io import StringIO


# non centered, sparse centers to check the
centers = np.array([
    [0.0, 5.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 4.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 5.0, 1.0],
])
n_samples = 100
n_clusters, n_features = centers.shape
X, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                            cluster_std=1., random_state=42)
X_csr = sp.csr_matrix(X)


def test_minibatch_update_consistency():
    # Check that dense and sparse minibatch update give the same results
    rng = np.random.RandomState(42)
    old_centers = centers + rng.normal(size=centers.shape)

    new_centers = old_centers.copy()
    new_centers_csr = old_centers.copy()

    weight_sums = np.zeros(new_centers.shape[0], dtype=np.double)
    weight_sums_csr = np.zeros(new_centers.shape[0], dtype=np.double)

    x_squared_norms = (X ** 2).sum(axis=1)
    x_squared_norms_csr = row_norms(X_csr, squared=True)

    buffer = np.zeros(centers.shape[1], dtype=np.double)
    buffer_csr = np.zeros(centers.shape[1], dtype=np.double)

    # extract a small minibatch
    X_mb = X[:10]
    X_mb_csr = X_csr[:10]
    x_mb_squared_norms = x_squared_norms[:10]
    x_mb_squared_norms_csr = x_squared_norms_csr[:10]

    sample_weight_mb = np.ones(X_mb.shape[0], dtype=np.double)

    # step 1: compute the dense minibatch update
    old_inertia, incremental_diff = _mini_batch_step(
        X_mb, sample_weight_mb, x_mb_squared_norms, new_centers, weight_sums,
        buffer, 1, np.random.RandomState(0), random_reassign=False)
    assert old_inertia > 0.0

    # compute the new inertia on the same batch to check that it decreased
    labels, new_inertia = _labels_inertia(
        X_mb, sample_weight_mb, x_mb_squared_norms, new_centers)
    assert new_inertia > 0.0
    assert new_inertia < old_inertia

    # check that the incremental difference computation is matching the
    # final observed value
    effective_diff = np.sum((new_centers - old_centers) ** 2)
    assert_almost_equal(incremental_diff, effective_diff)

    # step 2: compute the sparse minibatch update
    old_inertia_csr, incremental_diff_csr = _mini_batch_step(
        X_mb_csr, sample_weight_mb, x_mb_squared_norms_csr, new_centers_csr,
        weight_sums_csr, buffer_csr, 1, np.random.RandomState(0),
        random_reassign=False)
    assert old_inertia_csr > 0.0

    # compute the new inertia on the same batch to check that it decreased
    labels_csr, new_inertia_csr = _labels_inertia(
        X_mb_csr, sample_weight_mb, x_mb_squared_norms_csr, new_centers_csr)
    assert new_inertia_csr > 0.0
    assert new_inertia_csr < old_inertia_csr

    # check that the incremental difference computation is matching the
    # final observed value
    effective_diff = np.sum((new_centers_csr - old_centers) ** 2)
    assert_almost_equal(incremental_diff_csr, effective_diff)

    # step 3: check that sparse and dense updates lead to the same results
    assert_array_equal(labels, labels_csr)
    assert_array_almost_equal(new_centers, new_centers_csr)
    assert_almost_equal(incremental_diff, incremental_diff_csr)
    assert_almost_equal(old_inertia, old_inertia_csr)
    assert_almost_equal(new_inertia, new_inertia_csr)


def test_minibatch_with_many_reassignments():
    # Test for the case that the number of clusters to reassign is bigger
    # than the batch_size
    n_samples = 1000
    rnd = np.random.RandomState(42)
    X = rnd.uniform(size=(n_samples, 10))
    # Check that the fit works if n_clusters is bigger than the batch_size.
    # Run the test with 550 clusters and 550 samples, because it turned out
    # that this values ensure that the number of clusters to reassign
    # is always bigger than the batch_size
    n_clusters = 1000
    MiniBatchKMeans(n_clusters=n_clusters,
                    batch_size=100,
                    init_size=n_samples,
                    random_state=42,
                    verbose=True).fit(X)
    assert False
