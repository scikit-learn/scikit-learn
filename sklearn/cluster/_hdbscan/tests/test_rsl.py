"""
Tests for Robust Single Linkage clustering algorithm
"""
# import pickle
import numpy as np
from scipy.spatial import distance
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils._testing import assert_raises
from sklearn.cluster import RobustSingleLinkage, robust_single_linkage

from sklearn.datasets import make_blobs
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler

import pytest

n_clusters = 3
X, y = make_blobs(n_samples=50, random_state=1)
X, y = shuffle(X, y, random_state=7)
X = StandardScaler().fit_transform(X)
# X = generate_clustered_data(n_clusters=n_clusters, n_samples_per_cluster=50)


def test_rsl_distance_matrix():
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)

    labels = robust_single_linkage(D, 0.4, metric="precomputed")
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)  # ignore noise
    assert n_clusters_1 == 2

    labels = RobustSingleLinkage(metric="precomputed").fit(D).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == 2


def test_rsl_feature_vector():
    labels = robust_single_linkage(X, 0.4)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = RobustSingleLinkage().fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_rsl_callable_metric():
    # metric is the function reference, not the string key.
    metric = distance.euclidean

    labels = robust_single_linkage(X, 0.4, metric=metric)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = RobustSingleLinkage(metric=metric).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_rsl_input_lists():
    X = [[1.0, 2.0], [3.0, 4.0]]
    RobustSingleLinkage().fit(X)  # must not raise exception


def test_rsl_boruvka_balltree():
    labels = robust_single_linkage(X, 0.45, algorithm="boruvka_balltree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = RobustSingleLinkage(cut=0.45, algorithm="boruvka_balltree").fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_rsl_prims_balltree():
    labels = robust_single_linkage(X, 0.4, algorithm="prims_balltree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = RobustSingleLinkage(algorithm="prims_balltree").fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_rsl_prims_kdtree():
    labels = robust_single_linkage(X, 0.4, algorithm="prims_kdtree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = RobustSingleLinkage(algorithm="prims_kdtree").fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_rsl_high_dimensional():
    H, y = make_blobs(n_samples=50, random_state=0, n_features=64)
    # H, y = shuffle(X, y, random_state=7)
    H = StandardScaler().fit_transform(H)
    labels = robust_single_linkage(H, 5.5)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = (
        RobustSingleLinkage(
            cut=5.5,
            algorithm="best",
            metric="seuclidean",
            metric_params={"V": np.ones(H.shape[1])},
        )
        .fit(H)
        .labels_
    )
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_rsl_badargs():
    assert_raises(ValueError, robust_single_linkage, "fail", 0.4)
    assert_raises(ValueError, robust_single_linkage, None, 0.4)
    assert_raises(ValueError, robust_single_linkage, X, 0.4, k="fail")
    assert_raises(ValueError, robust_single_linkage, X, 0.4, k=-1)
    assert_raises(ValueError, robust_single_linkage, X, 0.4, metric="imperial")
    assert_raises(ValueError, robust_single_linkage, X, 0.4, metric=None)
    assert_raises(ValueError, robust_single_linkage, X, 0.4, metric="minkowski", p=-1)
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="minkowski",
        p=-1,
        algorithm="prims_kdtree",
    )
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="minkowski",
        p=-1,
        algorithm="prims_balltree",
    )
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="minkowski",
        p=-1,
        algorithm="boruvka_balltree",
    )
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="precomputed",
        algorithm="boruvka_kdtree",
    )
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="precomputed",
        algorithm="prims_kdtree",
    )
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="precomputed",
        algorithm="prims_balltree",
    )
    assert_raises(
        ValueError,
        robust_single_linkage,
        X,
        0.4,
        metric="precomputed",
        algorithm="boruvka_balltree",
    )
    assert_raises(ValueError, robust_single_linkage, X, 0.4, alpha=-1)
    assert_raises(ValueError, robust_single_linkage, X, 0.4, alpha="fail")
    assert_raises(Exception, robust_single_linkage, X, 0.4, algorithm="something_else")
    assert_raises(TypeError, robust_single_linkage, X, 0.4, metric="minkowski", p=None)
    assert_raises(ValueError, robust_single_linkage, X, 0.4, leaf_size=0)
    assert_raises(ValueError, robust_single_linkage, X, 0.4, gamma=0)


# Disable for now -- need to refactor to meet newer standards
@pytest.mark.skip(reason="need to refactor to meet newer standards")
def test_rsl_is_sklearn_estimator():
    check_estimator(RobustSingleLinkage)
