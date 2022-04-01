"""
Tests for HDBSCAN clustering algorithm
Shamelessly based on (i.e. ripped off from) the DBSCAN test code
"""
import numpy as np
from scipy.spatial import distance
from scipy import sparse
from scipy import stats
from sklearn.utils._testing import (
    assert_array_almost_equal,
    assert_raises,
)
from sklearn.cluster import HDBSCAN, hdbscan
from sklearn.cluster._hdbscan._validity import validity_index

from sklearn.datasets import make_blobs
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from scipy.stats import mode
from sklearn.metrics.pairwise import _VALID_METRICS
from sklearn.neighbors import KDTree, BallTree
from tempfile import mkdtemp
import pytest

from sklearn import datasets

n_clusters = 3
# X = generate_clustered_data(n_clusters=n_clusters, n_samples_per_cluster=50)
X, y = make_blobs(n_samples=200, random_state=10)
X, y = shuffle(X, y, random_state=7)
X = StandardScaler().fit_transform(X)

X_missing_data = X.copy()
X_missing_data[0] = [np.nan, 1]
X_missing_data[5] = [np.nan, np.nan]


def test_missing_data():
    """
    Tests if nan data are treated as infinite distance from all other points
    and assigned to -1 cluster.
    """
    model = HDBSCAN().fit(X_missing_data)
    assert model.labels_[0] == -1
    assert model.labels_[5] == -1
    assert model.probabilities_[0] == 0
    assert model.probabilities_[5] == 0
    assert model.probabilities_[5] == 0
    clean_indices = list(range(1, 5)) + list(range(6, 200))
    clean_model = HDBSCAN().fit(X_missing_data[clean_indices])
    assert np.allclose(clean_model.labels_, model.labels_[clean_indices])


def generate_noisy_data():
    blobs, _ = datasets.make_blobs(
        n_samples=200, centers=[(-0.75, 2.25), (1.0, 2.0)], cluster_std=0.25
    )
    moons, _ = datasets.make_moons(n_samples=200, noise=0.05)
    noise = np.random.uniform(-1.0, 3.0, (50, 2))
    return np.vstack([blobs, moons, noise])


def homogeneity(labels1, labels2):
    num_missed = 0.0
    for label in set(labels1):
        matches = labels2[labels1 == label]
        match_mode = mode(matches)[0][0]
        num_missed += np.sum(matches != match_mode)

    for label in set(labels2):
        matches = labels1[labels2 == label]
        match_mode = mode(matches)[0][0]
        num_missed += np.sum(matches != match_mode)

    return num_missed / 2.0


def test_hdbscan_distance_matrix():
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)

    labels = hdbscan(D, metric="precomputed")[0]
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)  # ignore noise
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(metric="precomputed").fit(D).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    validity = validity_index(D, labels, metric="precomputed", d=2)
    assert validity >= 0.6


def test_hdbscan_sparse_distance_matrix():
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)

    threshold = stats.scoreatpercentile(D.flatten(), 50)

    D[D >= threshold] = 0.0
    D = sparse.csr_matrix(D)
    D.eliminate_zeros()

    labels = hdbscan(D, metric="precomputed")[0]
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)  # ignore noise
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(metric="precomputed").fit(D).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_feature_vector():
    labels = hdbscan(X)[0]
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN().fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    validity = validity_index(X, labels)
    assert validity >= 0.4


@pytest.mark.parametrize(
    "algo",
    [
        "prims_kdtree",
        "prims_balltree",
        "boruvka_kdtree",
        "boruvka_balltree",
        "generic",
        "best",
    ],
)
@pytest.mark.parametrize("metric", _VALID_METRICS)
def test_hdbscan_algorithms(algo, metric):
    labels = hdbscan(X, algorithm=algo)[0]
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(algorithm=algo).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    ALGOS_TREES = {
        "prims_kdtree": KDTree,
        "prims_balltree": BallTree,
        "boruvka_kdtree": KDTree,
        "boruvka_balltree": BallTree,
    }
    METRIC_PARAMS = {
        "mahalanobis": {"V": np.eye(X.shape[1])},
        "seuclidean": {"V": np.ones(X.shape[1])},
        "minkowski": {"p": 2},
        "wminkowski": {"p": 2, "w": np.ones(X.shape[1])},
    }
    if algo not in ("best", "generic"):
        if metric not in ALGOS_TREES[algo].valid_metrics:
            with pytest.raises(ValueError):
                hdbscan(
                    X,
                    algorithm=algo,
                    metric=metric,
                    metric_params=METRIC_PARAMS.get(metric, None),
                )
        elif metric == "wminkowski":
            with pytest.warns(FutureWarning):
                hdbscan(
                    X,
                    algorithm=algo,
                    metric=metric,
                    metric_params=METRIC_PARAMS.get(metric, None),
                )
        else:
            hdbscan(
                X,
                algorithm=algo,
                metric=metric,
                metric_params=METRIC_PARAMS.get(metric, None),
            )


def test_hdbscan_dbscan_clustering():
    clusterer = HDBSCAN().fit(X)
    labels = clusterer.dbscan_clustering(0.3)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters == n_clusters_1


def test_hdbscan_high_dimensional():
    H, y = make_blobs(n_samples=50, random_state=0, n_features=64)
    # H, y = shuffle(X, y, random_state=7)
    H = StandardScaler().fit_transform(H)
    labels = hdbscan(H)[0]
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = (
        HDBSCAN(
            algorithm="best",
            metric="seuclidean",
            metric_params={"V": np.ones(H.shape[1])},
        )
        .fit(H)
        .labels_
    )
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_best_balltree_metric():
    labels = hdbscan(X, metric="seuclidean", metric_params={"V": np.ones(X.shape[1])})[
        0
    ]
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = (
        HDBSCAN(metric="seuclidean", metric_params={"V": np.ones(X.shape[1])})
        .fit(X)
        .labels_
    )
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_no_clusters():
    labels = hdbscan(X, min_cluster_size=len(X) + 1)[0]
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == 0

    labels = HDBSCAN(min_cluster_size=len(X) + 1).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == 0


def test_hdbscan_min_cluster_size():
    for min_cluster_size in range(2, len(X) + 1, 1):
        labels = hdbscan(X, min_cluster_size=min_cluster_size)[0]
        true_labels = [label for label in labels if label != -1]
        if len(true_labels) != 0:
            assert np.min(np.bincount(true_labels)) >= min_cluster_size

        labels = HDBSCAN(min_cluster_size=min_cluster_size).fit(X).labels_
        true_labels = [label for label in labels if label != -1]
        if len(true_labels) != 0:
            assert np.min(np.bincount(true_labels)) >= min_cluster_size


def test_hdbscan_callable_metric():
    metric = distance.euclidean

    labels = hdbscan(X, metric=metric)[0]
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(metric=metric).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_input_lists():
    X = [[1.0, 2.0], [3.0, 4.0]]
    HDBSCAN().fit(X)


@pytest.mark.parametrize("tree", ["kdtree", "balltree"])
def test_hdbscan_boruvka_matches(tree):

    data = generate_noisy_data()

    labels_prims = hdbscan(data, algorithm="generic")[0]
    labels_boruvka = hdbscan(data, algorithm=f"boruvka_{tree}")[0]

    num_mismatches = homogeneity(labels_prims, labels_boruvka)

    assert (num_mismatches / float(data.shape[0])) < 0.15

    labels_prims = HDBSCAN(algorithm="generic").fit_predict(data)
    labels_boruvka = HDBSCAN(algorithm=f"boruvka_{tree}").fit_predict(data)

    num_mismatches = homogeneity(labels_prims, labels_boruvka)

    assert (num_mismatches / float(data.shape[0])) < 0.15


def test_hdbscan_badargs():
    assert_raises(ValueError, hdbscan, X="fail")
    assert_raises(ValueError, hdbscan, X=None)
    assert_raises(ValueError, hdbscan, X, min_cluster_size="fail")
    assert_raises(ValueError, hdbscan, X, min_samples="fail")
    assert_raises(ValueError, hdbscan, X, min_samples=-1)
    assert_raises(ValueError, hdbscan, X, metric="imperial")
    assert_raises(ValueError, hdbscan, X, metric=None)
    assert_raises(
        ValueError, hdbscan, X, metric="precomputed", algorithm="boruvka_kdtree"
    )
    assert_raises(
        ValueError, hdbscan, X, metric="precomputed", algorithm="prims_kdtree"
    )
    assert_raises(
        ValueError, hdbscan, X, metric="precomputed", algorithm="prims_balltree"
    )
    assert_raises(
        ValueError, hdbscan, X, metric="precomputed", algorithm="boruvka_balltree"
    )
    assert_raises(ValueError, hdbscan, X, alpha=-1)
    assert_raises(ValueError, hdbscan, X, alpha="fail")
    assert_raises(Exception, hdbscan, X, algorithm="something_else")
    assert_raises(TypeError, hdbscan, X, metric="minkowski", p=None)
    assert_raises(ValueError, hdbscan, X, leaf_size=0)


def test_hdbscan_sparse():

    sparse_X = sparse.csr_matrix(X)

    labels = HDBSCAN().fit(sparse_X).labels_
    n_clusters = len(set(labels)) - int(-1 in labels)
    assert n_clusters == 3


def test_hdbscan_caching():

    cachedir = mkdtemp()
    labels1 = HDBSCAN(memory=cachedir, min_samples=5).fit(X).labels_
    labels2 = HDBSCAN(memory=cachedir, min_samples=5, min_cluster_size=6).fit(X).labels_
    n_clusters1 = len(set(labels1)) - int(-1 in labels1)
    n_clusters2 = len(set(labels2)) - int(-1 in labels2)
    assert n_clusters1 == n_clusters2


def test_hdbscan_centroids_medoids():
    centers = [(0.0, 0.0), (3.0, 3.0)]
    H, _ = make_blobs(n_samples=1000, random_state=0, centers=centers, cluster_std=0.5)
    clusterer = HDBSCAN().fit(H)

    for idx, center in enumerate(centers):
        centroid = clusterer.weighted_cluster_centroid(idx)
        assert_array_almost_equal(centroid, center, decimal=1)

        medoid = clusterer.weighted_cluster_medoid(idx)
        assert_array_almost_equal(medoid, center, decimal=1)


def test_hdbscan_no_centroid_medoid_for_noise():
    clusterer = HDBSCAN().fit(X)
    assert_raises(ValueError, clusterer.weighted_cluster_centroid, -1)
    assert_raises(ValueError, clusterer.weighted_cluster_medoid, -1)


def test_hdbscan_allow_single_cluster_with_epsilon():
    np.random.seed(0)
    no_structure = np.random.rand(150, 2)
    # without epsilon we should see many noise points as children of root.
    labels = HDBSCAN(
        min_cluster_size=5,
        cluster_selection_epsilon=0.0,
        cluster_selection_method="eom",
        allow_single_cluster=True,
    ).fit_predict(no_structure)
    unique_labels, counts = np.unique(labels, return_counts=True)
    assert len(unique_labels) == 2
    assert counts[unique_labels == -1] == 46

    # for this random seed an epsilon of 0.2 will produce exactly 2 noise
    # points at that cut in single linkage.
    labels = HDBSCAN(
        min_cluster_size=5,
        cluster_selection_epsilon=0.2,
        cluster_selection_method="eom",
        allow_single_cluster=True,
    ).fit_predict(no_structure)
    unique_labels, counts = np.unique(labels, return_counts=True)
    assert len(unique_labels) == 2
    assert counts[unique_labels == -1] == 2
