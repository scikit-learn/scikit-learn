"""
Tests for HDBSCAN clustering algorithm
Based on the DBSCAN test code
"""
import numpy as np
import pytest
from scipy import sparse, stats
from scipy.spatial import distance
from scipy.stats import mode

from sklearn import datasets
from sklearn.cluster import HDBSCAN
from sklearn.datasets import make_blobs
from sklearn.metrics import fowlkes_mallows_score
from sklearn.metrics.pairwise import _VALID_METRICS
from sklearn.neighbors import BallTree, KDTree
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle
from sklearn.utils._testing import assert_array_almost_equal

n_clusters_true = 3
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
    rng = np.random.RandomState(0)
    blobs, _ = datasets.make_blobs(
        n_samples=200, centers=[(-0.75, 2.25), (1.0, 2.0)], cluster_std=0.25
    )
    moons, _ = datasets.make_moons(n_samples=200, noise=0.05)
    noise = rng.uniform(-1.0, 3.0, (50, 2))
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

    labels = HDBSCAN(metric="precomputed").fit(D).labels_
    n_clusters_ = len(set(labels)) - int(-1 in labels)
    assert n_clusters_ == n_clusters_true

    # Check that clustering is arbitrarily good
    # This is a heuristic to guard against regression
    score = fowlkes_mallows_score(y, labels)
    assert score >= 0.98


def test_hdbscan_sparse_distance_matrix():
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)

    threshold = stats.scoreatpercentile(D.flatten(), 50)

    D[D >= threshold] = 0.0
    D = sparse.csr_matrix(D)
    D.eliminate_zeros()

    labels = HDBSCAN(metric="precomputed").fit(D).labels_
    n_clusters = len(set(labels)) - int(-1 in labels)
    assert n_clusters == n_clusters_true


def test_hdbscan_feature_vector():
    labels = HDBSCAN().fit_predict(X)
    n_clusters = len(set(labels)) - int(-1 in labels)
    assert n_clusters == n_clusters_true

    # Check that clustering is arbitrarily good
    # This is a heuristic to guard against regression
    score = fowlkes_mallows_score(y, labels)
    assert score >= 0.98


@pytest.mark.parametrize(
    "algo",
    [
        "kdtree",
        "balltree",
        "brute",
        "auto",
    ],
)
@pytest.mark.parametrize("metric", _VALID_METRICS)
def test_hdbscan_algorithms(algo, metric):
    labels = HDBSCAN(algorithm=algo).fit_predict(X)
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters_true

    ALGOS_TREES = {
        "kdtree": KDTree,
        "balltree": BallTree,
    }
    metric_params = {
        "mahalanobis": {"V": np.eye(X.shape[1])},
        "seuclidean": {"V": np.ones(X.shape[1])},
        "minkowski": {"p": 2},
        "wminkowski": {"p": 2, "w": np.ones(X.shape[1])},
    }.get(metric, None)
    if algo not in ("auto", "brute"):
        if metric not in ALGOS_TREES[algo].valid_metrics:
            with pytest.raises(ValueError):
                HDBSCAN(
                    algorithm=algo,
                    metric=metric,
                    metric_params=metric_params,
                ).fit(X)
        elif metric == "wminkowski":
            with pytest.warns(FutureWarning):
                HDBSCAN(
                    algorithm=algo,
                    metric=metric,
                    metric_params=metric_params,
                ).fit(X)
        else:
            HDBSCAN(
                algorithm=algo,
                metric=metric,
                metric_params=metric_params,
            ).fit(X)


def test_hdbscan_dbscan_clustering():
    clusterer = HDBSCAN().fit(X)
    labels = clusterer.dbscan_clustering(0.3)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_true == n_clusters_1


def test_hdbscan_high_dimensional():
    H, y = make_blobs(n_samples=50, random_state=0, n_features=64)
    H = StandardScaler().fit_transform(H)
    labels = HDBSCAN(
        algorithm="auto",
        metric="seuclidean",
        metric_params={"V": np.ones(H.shape[1])},
    ).fit_predict(H)
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters_true


def test_hdbscan_best_balltree_metric():
    labels = HDBSCAN(
        metric="seuclidean", metric_params={"V": np.ones(X.shape[1])}
    ).fit_predict(X)
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters_true


def test_hdbscan_no_clusters():
    labels = HDBSCAN(min_cluster_size=len(X) - 1).fit_predict(X)
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == 0


def test_hdbscan_min_cluster_size():
    """
    Test that the smallest non-noise cluster has at least `min_cluster_size`
    many points
    """
    for min_cluster_size in range(2, len(X), 1):
        labels = HDBSCAN(min_cluster_size=min_cluster_size).fit_predict(X)
        true_labels = [label for label in labels if label != -1]
        if len(true_labels) != 0:
            assert np.min(np.bincount(true_labels)) >= min_cluster_size


def test_hdbscan_callable_metric():
    metric = distance.euclidean
    labels = HDBSCAN(metric=metric).fit_predict(X)
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters_true


def test_hdbscan_input_lists():
    X = [[1.0, 2.0], [3.0, 4.0]]
    HDBSCAN(min_samples=1).fit(X)


@pytest.mark.parametrize("tree", ["kd", "ball"])
def test_hdbscan_precomputed_non_brute(tree):
    hdb = HDBSCAN(metric="precomputed", algorithm=f"prims_{tree}tree")
    with pytest.raises(ValueError):
        hdb.fit(X)


def test_hdbscan_sparse():

    sparse_X = sparse.csr_matrix(X)

    labels = HDBSCAN().fit(sparse_X).labels_
    n_clusters = len(set(labels)) - int(-1 in labels)
    assert n_clusters == 3

    sparse_X_nan = sparse_X.copy()
    sparse_X_nan[0, 0] = np.nan
    labels = HDBSCAN().fit(sparse_X_nan).labels_
    n_clusters = len(set(labels)) - int(-1 in labels)
    assert n_clusters == 3

    msg = "Sparse data matrices only support algorithm `brute`."
    with pytest.raises(ValueError, match=msg):
        HDBSCAN(metric="euclidean", algorithm="balltree").fit(sparse_X)


def test_hdbscan_caching(tmp_path):

    labels1 = HDBSCAN(memory=tmp_path, min_samples=5).fit_predict(X)
    labels2 = HDBSCAN(memory=tmp_path, min_samples=5, min_cluster_size=6).fit_predict(X)
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
    with pytest.raises(ValueError):
        clusterer.weighted_cluster_centroid(-1)
    with pytest.raises(ValueError):
        clusterer.weighted_cluster_medoid(-1)


def test_hdbscan_allow_single_cluster_with_epsilon():
    rng = np.random.RandomState(0)
    no_structure = rng.rand(150, 2)
    # without epsilon we should see many noise points as children of root.
    labels = HDBSCAN(
        min_cluster_size=5,
        cluster_selection_epsilon=0.0,
        cluster_selection_method="eom",
        allow_single_cluster=True,
    ).fit_predict(no_structure)
    unique_labels, counts = np.unique(labels, return_counts=True)
    assert len(unique_labels) == 2

    # Arbitrary heuristic. Would prefer something more precise.
    assert counts[unique_labels == -1] > 30

    # for this random seed an epsilon of 0.18 will produce exactly 2 noise
    # points at that cut in single linkage.
    labels = HDBSCAN(
        min_cluster_size=5,
        cluster_selection_epsilon=0.18,
        cluster_selection_method="eom",
        allow_single_cluster=True,
        algorithm="kdtree",
    ).fit_predict(no_structure)
    unique_labels, counts = np.unique(labels, return_counts=True)
    assert len(unique_labels) == 2
    assert counts[unique_labels == -1] == 2


def test_hdbscan_better_than_dbscan():
    """
    Validate that HDBSCAN can properly cluster this difficult synthetic
    dataset. Note that DBSCAN fails on this (see HDBSCAN plotting
    example)
    """
    centers = [[-0.85, -0.85], [-0.85, 0.85], [3, 3], [3, -3]]
    X, _ = make_blobs(
        n_samples=750,
        centers=centers,
        cluster_std=[0.2, 0.35, 1.35, 1.35],
        random_state=0,
    )
    hdb = HDBSCAN().fit(X)
    n_clusters = len(set(hdb.labels_)) - int(-1 in hdb.labels_)
    assert n_clusters == 4


def test_hdbscan_unfit_centers_errors():
    hdb = HDBSCAN()
    msg = "Model has not been fit to data"
    with pytest.raises(AttributeError, match=msg):
        hdb.weighted_cluster_centroid(0)
    with pytest.raises(AttributeError, match=msg):
        hdb.weighted_cluster_medoid(0)


def test_hdbscan_precomputed_array_like():
    X = np.array([[1, np.inf], [np.inf, 1]])
    HDBSCAN(metric="precomputed").fit(X)


def test_hdbscan_sparse_distances_too_few_nonzero():
    X = sparse.csr_matrix(np.zeros((10, 10)))

    msg = "There exists points with less than"
    with pytest.raises(ValueError, match=msg):
        HDBSCAN(metric="precomputed").fit(X)
