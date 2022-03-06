"""
Tests for HDBSCAN clustering algorithm
Shamelessly based on (i.e. ripped off from) the DBSCAN test code
"""
import numpy as np
from scipy.spatial import distance
from scipy import sparse
from scipy import stats
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils._testing import (
    assert_array_equal,
    assert_array_almost_equal,
    assert_raises,
)
from sklearn.cluster import (
    HDBSCAN,
    hdbscan,
    validity_index,
    approximate_predict,
    approximate_predict_scores,
    all_points_membership_vectors,
)

# from sklearn.cluster.tests.common import generate_clustered_data
from sklearn.datasets import make_blobs
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from scipy.stats import mode

from tempfile import mkdtemp
import pytest

from sklearn import datasets

import warnings

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

    labels, p, persist, ctree, ltree, mtree = hdbscan(D, metric="precomputed")
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

    labels, p, persist, ctree, ltree, mtree = hdbscan(D, metric="precomputed")
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)  # ignore noise
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(metric="precomputed", gen_min_span_tree=True).fit(D).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_feature_vector():
    labels, p, persist, ctree, ltree, mtree = hdbscan(X)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN().fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    validity = validity_index(X, labels)
    assert validity >= 0.4


def test_hdbscan_prims_kdtree():
    labels, p, persist, ctree, ltree, mtree = hdbscan(X, algorithm="prims_kdtree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(algorithm="prims_kdtree", gen_min_span_tree=True).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    assert_raises(ValueError, hdbscan, X, algorithm="prims_kdtree", metric="russelrao")


def test_hdbscan_prims_balltree():
    labels, p, persist, ctree, ltree, mtree = hdbscan(X, algorithm="prims_balltree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(algorithm="prims_balltree", gen_min_span_tree=True).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    assert_raises(ValueError, hdbscan, X, algorithm="prims_balltree", metric="cosine")


def test_hdbscan_boruvka_kdtree():
    labels, p, persist, ctree, ltree, mtree = hdbscan(X, algorithm="boruvka_kdtree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(algorithm="boruvka_kdtree", gen_min_span_tree=True).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    assert_raises(
        ValueError, hdbscan, X, algorithm="boruvka_kdtree", metric="russelrao"
    )


def test_hdbscan_boruvka_balltree():
    labels, p, persist, ctree, ltree, mtree = hdbscan(X, algorithm="boruvka_balltree")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = (
        HDBSCAN(algorithm="boruvka_balltree", gen_min_span_tree=True).fit(X).labels_
    )
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    assert_raises(ValueError, hdbscan, X, algorithm="boruvka_balltree", metric="cosine")


def test_hdbscan_generic():
    labels, p, persist, ctree, ltree, mtree = hdbscan(X, algorithm="generic")
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(algorithm="generic", gen_min_span_tree=True).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_dbscan_clustering():
    clusterer = HDBSCAN().fit(X)
    labels = clusterer.dbscan_clustering(0.3)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters == n_clusters_1


def test_hdbscan_high_dimensional():
    H, y = make_blobs(n_samples=50, random_state=0, n_features=64)
    # H, y = shuffle(X, y, random_state=7)
    H = StandardScaler().fit_transform(H)
    labels, p, persist, ctree, ltree, mtree = hdbscan(H)
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
    labels, p, persist, ctree, ltree, mtree = hdbscan(
        X, metric="seuclidean", metric_params={"V": np.ones(X.shape[1])}
    )
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
    labels, p, persist, ctree, ltree, mtree = hdbscan(X, min_cluster_size=len(X) + 1)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == 0

    labels = HDBSCAN(min_cluster_size=len(X) + 1).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == 0


def test_hdbscan_min_cluster_size():
    for min_cluster_size in range(2, len(X) + 1, 1):
        labels, p, persist, ctree, ltree, mtree = hdbscan(
            X, min_cluster_size=min_cluster_size
        )
        true_labels = [label for label in labels if label != -1]
        if len(true_labels) != 0:
            assert np.min(np.bincount(true_labels)) >= min_cluster_size

        labels = HDBSCAN(min_cluster_size=min_cluster_size).fit(X).labels_
        true_labels = [label for label in labels if label != -1]
        if len(true_labels) != 0:
            assert np.min(np.bincount(true_labels)) >= min_cluster_size


def test_hdbscan_callable_metric():
    # metric is the function reference, not the string key.
    metric = distance.euclidean

    labels, p, persist, ctree, ltree, mtree = hdbscan(X, metric=metric)
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    labels = HDBSCAN(metric=metric).fit(X).labels_
    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_hdbscan_input_lists():
    X = [[1.0, 2.0], [3.0, 4.0]]
    HDBSCAN().fit(X)  # must not raise exception


def test_hdbscan_boruvka_kdtree_matches():

    data = generate_noisy_data()

    labels_prims, p, persist, ctree, ltree, mtree = hdbscan(data, algorithm="generic")
    labels_boruvka, p, persist, ctree, ltree, mtree = hdbscan(
        data, algorithm="boruvka_kdtree"
    )

    num_mismatches = homogeneity(labels_prims, labels_boruvka)

    assert (num_mismatches / float(data.shape[0])) < 0.15

    labels_prims = HDBSCAN(algorithm="generic").fit_predict(data)
    labels_boruvka = HDBSCAN(algorithm="boruvka_kdtree").fit_predict(data)

    num_mismatches = homogeneity(labels_prims, labels_boruvka)

    assert (num_mismatches / float(data.shape[0])) < 0.15


def test_hdbscan_boruvka_balltree_matches():

    data = generate_noisy_data()

    labels_prims, p, persist, ctree, ltree, mtree = hdbscan(data, algorithm="generic")
    labels_boruvka, p, persist, ctree, ltree, mtree = hdbscan(
        data, algorithm="boruvka_balltree"
    )

    num_mismatches = homogeneity(labels_prims, labels_boruvka)

    assert (num_mismatches / float(data.shape[0])) < 0.15

    labels_prims = HDBSCAN(algorithm="generic").fit_predict(data)
    labels_boruvka = HDBSCAN(algorithm="boruvka_balltree").fit_predict(data)

    num_mismatches = homogeneity(labels_prims, labels_boruvka)

    assert (num_mismatches / float(data.shape[0])) < 0.15


def test_tree_numpy_output_formats():

    clusterer = HDBSCAN(gen_min_span_tree=True).fit(X)

    clusterer.single_linkage_tree_.to_numpy()
    clusterer.condensed_tree_.to_numpy()
    clusterer.minimum_spanning_tree_.to_numpy()


def test_hdbscan_outliers():
    clusterer = HDBSCAN(gen_min_span_tree=True).fit(X)
    scores = clusterer.outlier_scores_
    assert scores is not None


# def test_hdbscan_unavailable_attributes():
#     clusterer = HDBSCAN(gen_min_span_tree=False)
#     with warnings.catch_warnings(record=True) as w:
#         tree = clusterer.condensed_tree_
#         assert len(w) > 0
#         assert tree is None
#     with warnings.catch_warnings(record=True) as w:
#         tree = clusterer.single_linkage_tree_
#         assert len(w) > 0
#         assert tree is None
#     with warnings.catch_warnings(record=True) as w:
#         scores = clusterer.outlier_scores_
#         assert len(w) > 0
#         assert scores is None
#     with warnings.catch_warnings(record=True) as w:
#         tree = clusterer.minimum_spanning_tree_
#         assert len(w) > 0
#         assert tree is None


# def test_hdbscan_min_span_tree_availability():
#     clusterer = HDBSCAN().fit(X)
#     tree = clusterer.minimum_spanning_tree_
#     assert tree is None
#     D = distance.squareform(distance.pdist(X))
#     D /= np.max(D)
#     HDBSCAN(metric='precomputed').fit(D)
#     tree = clusterer.minimum_spanning_tree_
#     assert tree is None


def test_hdbscan_approximate_predict():
    clusterer = HDBSCAN(prediction_data=True).fit(X)
    cluster, prob = approximate_predict(clusterer, np.array([[-1.5, -1.0]]))
    assert cluster == 2
    cluster, prob = approximate_predict(clusterer, np.array([[1.5, -1.0]]))
    assert cluster == 1
    cluster, prob = approximate_predict(clusterer, np.array([[0.0, 0.0]]))
    assert cluster == -1


def test_hdbscan_approximate_predict_score():
    clusterer = HDBSCAN(min_cluster_size=200).fit(X)
    # no prediction data error
    assert_raises(ValueError, approximate_predict_scores, clusterer, X)
    clusterer.generate_prediction_data()
    # wrong dimensions error
    assert_raises(
        ValueError, approximate_predict_scores, clusterer, np.array([[1, 2, 3]])
    )
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        approximate_predict_scores(clusterer, np.array([[1.5, -1.0]]))
        # no clusters warning
        assert "Clusterer does not have any defined clusters" in str(w[-1].message)
    clusterer = HDBSCAN(prediction_data=True).fit(X)
    scores = approximate_predict_scores(clusterer, X)
    assert_array_almost_equal(scores, clusterer.outlier_scores_)
    assert scores.min() >= 0
    assert scores.max() <= 1


# def test_hdbscan_membership_vector():
#     clusterer = HDBSCAN(prediction_data=True).fit(X)
#     vector = membership_vector(clusterer, np.array([[-1.5, -1.0]]))
#     assert_array_almost_equal(
#         vector,
#         np.array([[ 0.05705305,  0.05974177,  0.12228153]]))
#     vector = membership_vector(clusterer, np.array([[1.5, -1.0]]))
#     assert_array_almost_equal(
#         vector,
#         np.array([[ 0.09462176,  0.32061556,  0.10112905]]))
#     vector = membership_vector(clusterer, np.array([[0.0, 0.0]]))
#     assert_array_almost_equal(
#         vector,
#         np.array([[ 0.03545607,  0.03363318,  0.04643177]]))
#
# def test_hdbscan_all_points_membership_vectors():
#     clusterer = HDBSCAN(prediction_data=True).fit(X)
#     vects = all_points_membership_vectors(clusterer)
#     assert_array_almost_equal(vects[0], np.array([7.86400992e-002,
#                                                    2.52734246e-001,
#                                                    8.38299608e-002]))
#     assert_array_almost_equal(vects[-1], np.array([8.09055344e-001,
#                                                    8.35882503e-002,
#                                                    1.07356406e-001]))


def test_hdbscan_all_points_membership_vectors():
    clusterer = HDBSCAN(prediction_data=True, min_cluster_size=200).fit(X)
    vects = all_points_membership_vectors(clusterer)
    assert_array_equal(vects, np.zeros(clusterer.prediction_data_.raw_data.shape[0]))


def test_hdbscan_badargs():
    assert_raises(ValueError, hdbscan, X="fail")
    assert_raises(ValueError, hdbscan, X=None)
    assert_raises(ValueError, hdbscan, X, min_cluster_size="fail")
    assert_raises(ValueError, hdbscan, X, min_samples="fail")
    assert_raises(ValueError, hdbscan, X, min_samples=-1)
    assert_raises(ValueError, hdbscan, X, metric="imperial")
    assert_raises(ValueError, hdbscan, X, metric=None)
    assert_raises(ValueError, hdbscan, X, metric="minkowski", p=-1)
    assert_raises(
        ValueError, hdbscan, X, metric="minkowski", p=-1, algorithm="prims_kdtree"
    )
    assert_raises(
        ValueError, hdbscan, X, metric="minkowski", p=-1, algorithm="prims_balltree"
    )
    assert_raises(
        ValueError, hdbscan, X, metric="minkowski", p=-1, algorithm="boruvka_balltree"
    )
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
    H, y = make_blobs(n_samples=1000, random_state=0, centers=centers, cluster_std=0.5)
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


# Disable for now -- need to refactor to meet newer standards
@pytest.mark.skip(reason="need to refactor to meet newer standards")
def test_hdbscan_is_sklearn_estimator():
    check_estimator(HDBSCAN)


# Probably not applicable now #
# def test_dbscan_sparse():
# def test_dbscan_balltree():
# def test_pickle():
# def test_dbscan_core_samples_toy():
# def test_boundaries():
