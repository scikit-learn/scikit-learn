import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.sparse import issparse

from sklearn import datasets
from sklearn.metrics import pairwise_distances
from sklearn.metrics.cluster import (
    calinski_harabasz_score,
    davies_bouldin_score,
    dbcv_score,
    silhouette_samples,
    silhouette_score,
)
from sklearn.metrics.cluster._dbcv_helper import (
    _check_duplicated_samples,
    compute_cluster_core_distance,
    compute_mutual_reach_dists,
    compute_pair_to_pair_dists,
    fn_density_separation,
    fn_density_sparseness,
    get_internal_objects,
    get_subarray,
)
from sklearn.metrics.cluster._unsupervised import _silhouette_reduce
from sklearn.utils._testing import assert_array_equal
from sklearn.utils.fixes import (
    CSC_CONTAINERS,
    CSR_CONTAINERS,
    DOK_CONTAINERS,
    LIL_CONTAINERS,
)


@pytest.mark.parametrize(
    "sparse_container",
    [None] + CSR_CONTAINERS + CSC_CONTAINERS + DOK_CONTAINERS + LIL_CONTAINERS,
)
@pytest.mark.parametrize("sample_size", [None, "half"])
def test_silhouette(sparse_container, sample_size):
    # Tests the Silhouette Coefficient.
    dataset = datasets.load_iris()
    X, y = dataset.data, dataset.target
    if sparse_container is not None:
        X = sparse_container(X)
    sample_size = int(X.shape[0] / 2) if sample_size == "half" else sample_size

    D = pairwise_distances(X, metric="euclidean")
    # Given that the actual labels are used, we can assume that S would be positive.
    score_precomputed = silhouette_score(
        D, y, metric="precomputed", sample_size=sample_size, random_state=0
    )
    score_euclidean = silhouette_score(
        X, y, metric="euclidean", sample_size=sample_size, random_state=0
    )
    assert score_precomputed > 0
    assert score_euclidean > 0
    assert score_precomputed == pytest.approx(score_euclidean)


def test_cluster_size_1():
    # Assert Silhouette Coefficient == 0 when there is 1 sample in a cluster
    # (cluster 0). We also test the case where there are identical samples
    # as the only members of a cluster (cluster 2). To our knowledge, this case
    # is not discussed in reference material, and we choose for it a sample
    # score of 1.
    X = [[0.0], [1.0], [1.0], [2.0], [3.0], [3.0]]
    labels = np.array([0, 1, 1, 1, 2, 2])

    # Cluster 0: 1 sample -> score of 0 by Rousseeuw's convention
    # Cluster 1: intra-cluster = [.5, .5, 1]
    #            inter-cluster = [1, 1, 1]
    #            silhouette    = [.5, .5, 0]
    # Cluster 2: intra-cluster = [0, 0]
    #            inter-cluster = [arbitrary, arbitrary]
    #            silhouette    = [1., 1.]

    silhouette = silhouette_score(X, labels)
    assert not np.isnan(silhouette)
    ss = silhouette_samples(X, labels)
    assert_array_equal(ss, [0, 0.5, 0.5, 0, 1, 1])


def test_silhouette_paper_example():
    # Explicitly check per-sample results against Rousseeuw (1987)
    # Data from Table 1
    lower = [
        5.58,
        7.00,
        6.50,
        7.08,
        7.00,
        3.83,
        4.83,
        5.08,
        8.17,
        5.83,
        2.17,
        5.75,
        6.67,
        6.92,
        4.92,
        6.42,
        5.00,
        5.58,
        6.00,
        4.67,
        6.42,
        3.42,
        5.50,
        6.42,
        6.42,
        5.00,
        3.92,
        6.17,
        2.50,
        4.92,
        6.25,
        7.33,
        4.50,
        2.25,
        6.33,
        2.75,
        6.08,
        6.67,
        4.25,
        2.67,
        6.00,
        6.17,
        6.17,
        6.92,
        6.17,
        5.25,
        6.83,
        4.50,
        3.75,
        5.75,
        5.42,
        6.08,
        5.83,
        6.67,
        3.67,
        4.75,
        3.00,
        6.08,
        6.67,
        5.00,
        5.58,
        4.83,
        6.17,
        5.67,
        6.50,
        6.92,
    ]
    D = np.zeros((12, 12))
    D[np.tril_indices(12, -1)] = lower
    D += D.T

    names = [
        "BEL",
        "BRA",
        "CHI",
        "CUB",
        "EGY",
        "FRA",
        "IND",
        "ISR",
        "USA",
        "USS",
        "YUG",
        "ZAI",
    ]

    # Data from Figure 2
    labels1 = [1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1]
    expected1 = {
        "USA": 0.43,
        "BEL": 0.39,
        "FRA": 0.35,
        "ISR": 0.30,
        "BRA": 0.22,
        "EGY": 0.20,
        "ZAI": 0.19,
        "CUB": 0.40,
        "USS": 0.34,
        "CHI": 0.33,
        "YUG": 0.26,
        "IND": -0.04,
    }
    score1 = 0.28

    # Data from Figure 3
    labels2 = [1, 2, 3, 3, 1, 1, 2, 1, 1, 3, 3, 2]
    expected2 = {
        "USA": 0.47,
        "FRA": 0.44,
        "BEL": 0.42,
        "ISR": 0.37,
        "EGY": 0.02,
        "ZAI": 0.28,
        "BRA": 0.25,
        "IND": 0.17,
        "CUB": 0.48,
        "USS": 0.44,
        "YUG": 0.31,
        "CHI": 0.31,
    }
    score2 = 0.33

    for labels, expected, score in [
        (labels1, expected1, score1),
        (labels2, expected2, score2),
    ]:
        expected = [expected[name] for name in names]
        # we check to 2dp because that's what's in the paper
        pytest.approx(
            expected,
            silhouette_samples(D, np.array(labels), metric="precomputed"),
            abs=1e-2,
        )
        pytest.approx(
            score, silhouette_score(D, np.array(labels), metric="precomputed"), abs=1e-2
        )


def test_correct_labelsize():
    # Assert 1 < n_labels < n_samples
    dataset = datasets.load_iris()
    X = dataset.data

    # n_labels = n_samples
    y = np.arange(X.shape[0])
    err_msg = (
        r"Number of labels is %d\. Valid values are 2 "
        r"to n_samples - 1 \(inclusive\)" % len(np.unique(y))
    )
    with pytest.raises(ValueError, match=err_msg):
        silhouette_score(X, y)

    # n_labels = 1
    y = np.zeros(X.shape[0])
    err_msg = (
        r"Number of labels is %d\. Valid values are 2 "
        r"to n_samples - 1 \(inclusive\)" % len(np.unique(y))
    )
    with pytest.raises(ValueError, match=err_msg):
        silhouette_score(X, y)


def test_non_encoded_labels():
    dataset = datasets.load_iris()
    X = dataset.data
    labels = dataset.target
    assert silhouette_score(X, labels * 2 + 10) == silhouette_score(X, labels)
    assert_array_equal(
        silhouette_samples(X, labels * 2 + 10), silhouette_samples(X, labels)
    )


def test_non_numpy_labels():
    dataset = datasets.load_iris()
    X = dataset.data
    y = dataset.target
    assert silhouette_score(list(X), list(y)) == silhouette_score(X, y)


@pytest.mark.parametrize("dtype", (np.float32, np.float64))
def test_silhouette_nonzero_diag(dtype):
    # Make sure silhouette_samples requires diagonal to be zero.
    # Non-regression test for #12178

    # Construct a zero-diagonal matrix
    dists = pairwise_distances(
        np.array([[0.2, 0.1, 0.12, 1.34, 1.11, 1.6]], dtype=dtype).T
    )
    labels = [0, 0, 0, 1, 1, 1]

    # small values on the diagonal are OK
    dists[2][2] = np.finfo(dists.dtype).eps * 10
    silhouette_samples(dists, labels, metric="precomputed")

    # values bigger than eps * 100 are not
    dists[2][2] = np.finfo(dists.dtype).eps * 1000
    with pytest.raises(ValueError, match="contains non-zero"):
        silhouette_samples(dists, labels, metric="precomputed")


@pytest.mark.parametrize(
    "sparse_container",
    CSC_CONTAINERS + CSR_CONTAINERS + DOK_CONTAINERS + LIL_CONTAINERS,
)
def test_silhouette_samples_precomputed_sparse(sparse_container):
    """Check that silhouette_samples works for sparse matrices correctly."""
    X = np.array([[0.2, 0.1, 0.1, 0.2, 0.1, 1.6, 0.2, 0.1]], dtype=np.float32).T
    y = [0, 0, 0, 0, 1, 1, 1, 1]
    pdist_dense = pairwise_distances(X)
    pdist_sparse = sparse_container(pdist_dense)
    assert issparse(pdist_sparse)
    output_with_sparse_input = silhouette_samples(pdist_sparse, y, metric="precomputed")
    output_with_dense_input = silhouette_samples(pdist_dense, y, metric="precomputed")
    assert_allclose(output_with_sparse_input, output_with_dense_input)


@pytest.mark.parametrize(
    "sparse_container",
    CSC_CONTAINERS + CSR_CONTAINERS + DOK_CONTAINERS + LIL_CONTAINERS,
)
def test_silhouette_samples_euclidean_sparse(sparse_container):
    """Check that silhouette_samples works for sparse matrices correctly."""
    X = np.array([[0.2, 0.1, 0.1, 0.2, 0.1, 1.6, 0.2, 0.1]], dtype=np.float32).T
    y = [0, 0, 0, 0, 1, 1, 1, 1]
    pdist_dense = pairwise_distances(X)
    pdist_sparse = sparse_container(pdist_dense)
    assert issparse(pdist_sparse)
    output_with_sparse_input = silhouette_samples(pdist_sparse, y)
    output_with_dense_input = silhouette_samples(pdist_dense, y)
    assert_allclose(output_with_sparse_input, output_with_dense_input)


@pytest.mark.parametrize(
    "sparse_container", CSC_CONTAINERS + DOK_CONTAINERS + LIL_CONTAINERS
)
def test_silhouette_reduce(sparse_container):
    """Check for non-CSR input to private method `_silhouette_reduce`."""
    X = np.array([[0.2, 0.1, 0.1, 0.2, 0.1, 1.6, 0.2, 0.1]], dtype=np.float32).T
    pdist_dense = pairwise_distances(X)
    pdist_sparse = sparse_container(pdist_dense)
    y = [0, 0, 0, 0, 1, 1, 1, 1]
    label_freqs = np.bincount(y)
    with pytest.raises(
        TypeError,
        match="Expected CSR matrix. Please pass sparse matrix in CSR format.",
    ):
        _silhouette_reduce(pdist_sparse, start=0, labels=y, label_freqs=label_freqs)


def assert_raises_on_only_one_label(func):
    """Assert message when there is only one label"""
    rng = np.random.RandomState(seed=0)
    with pytest.raises(ValueError, match="Number of labels is"):
        func(rng.rand(10, 2), np.zeros(10))


def assert_raises_on_all_points_same_cluster(func):
    """Assert message when all point are in different clusters"""
    rng = np.random.RandomState(seed=0)
    with pytest.raises(ValueError, match="Number of labels is"):
        func(rng.rand(10, 2), np.arange(10))


def test_calinski_harabasz_score():
    assert_raises_on_only_one_label(calinski_harabasz_score)

    assert_raises_on_all_points_same_cluster(calinski_harabasz_score)

    # Assert the value is 1. when all samples are equals
    assert 1.0 == calinski_harabasz_score(np.ones((10, 2)), [0] * 5 + [1] * 5)

    # Assert the value is 0. when all the mean cluster are equal
    assert 0.0 == calinski_harabasz_score([[-1, -1], [1, 1]] * 10, [0] * 10 + [1] * 10)

    # General case (with non numpy arrays)
    X = (
        [[0, 0], [1, 1]] * 5
        + [[3, 3], [4, 4]] * 5
        + [[0, 4], [1, 3]] * 5
        + [[3, 1], [4, 0]] * 5
    )
    labels = [0] * 10 + [1] * 10 + [2] * 10 + [3] * 10
    pytest.approx(calinski_harabasz_score(X, labels), 45 * (40 - 4) / (5 * (4 - 1)))


def test_davies_bouldin_score():
    assert_raises_on_only_one_label(davies_bouldin_score)
    assert_raises_on_all_points_same_cluster(davies_bouldin_score)

    # Assert the value is 0. when all samples are equals
    assert davies_bouldin_score(np.ones((10, 2)), [0] * 5 + [1] * 5) == pytest.approx(
        0.0
    )

    # Assert the value is 0. when all the mean cluster are equal
    assert davies_bouldin_score(
        [[-1, -1], [1, 1]] * 10, [0] * 10 + [1] * 10
    ) == pytest.approx(0.0)

    # General case (with non numpy arrays)
    X = (
        [[0, 0], [1, 1]] * 5
        + [[3, 3], [4, 4]] * 5
        + [[0, 4], [1, 3]] * 5
        + [[3, 1], [4, 0]] * 5
    )
    labels = [0] * 10 + [1] * 10 + [2] * 10 + [3] * 10
    pytest.approx(davies_bouldin_score(X, labels), 2 * np.sqrt(0.5) / 3)

    # Ensure divide by zero warning is not raised in general case
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        davies_bouldin_score(X, labels)

    # General case - cluster have one sample
    X = [[0, 0], [2, 2], [3, 3], [5, 5]]
    labels = [0, 0, 1, 2]
    pytest.approx(davies_bouldin_score(X, labels), (5.0 / 4) / 3)


def test_silhouette_score_integer_precomputed():
    """Check that silhouette_score works for precomputed metrics that are integers.

    Non-regression test for #22107.
    """
    result = silhouette_score(
        [[0, 1, 2], [1, 0, 1], [2, 1, 0]], [0, 0, 1], metric="precomputed"
    )
    assert result == pytest.approx(1 / 6)

    # non-zero on diagonal for ints raises an error
    with pytest.raises(ValueError, match="contains non-zero"):
        silhouette_score(
            [[1, 1, 2], [1, 0, 1], [2, 1, 0]], [0, 0, 1], metric="precomputed"
        )


@pytest.fixture
def sample_data():
    return np.array([[0, 1], [2, 3], [4, 5]])


@pytest.fixture
def sample_data_large():
    return np.arange(1000000).reshape(1000, 1000)


def test_compute_pair_to_pair_dists(sample_data):
    dists = compute_pair_to_pair_dists(sample_data, metric="euclidean")
    assert dists.shape == (3, 3)


def test_get_subarray(sample_data):
    subarray = get_subarray(sample_data, inds_a=[0, 2], inds_b=[1])
    expected_subarray = np.array([[1, 5]])

    assert subarray.shape == (1, 2)
    assert np.array_equal(subarray, expected_subarray)


def test_get_internal_objects(sample_data):
    mutual_reach_dists = compute_pair_to_pair_dists(sample_data, metric="euclidean")
    internal_node_inds, internal_edge_weights = get_internal_objects(mutual_reach_dists)

    expected_internal_node_inds = np.array([1])
    expected_internal_edge_weights = np.array([[0.0]])

    assert np.array_equal(internal_node_inds, expected_internal_node_inds)
    assert np.array_equal(internal_edge_weights, expected_internal_edge_weights)


def test_compute_cluster_core_distance(sample_data):
    dists = compute_pair_to_pair_dists(sample_data, metric="euclidean")
    core_dists = compute_cluster_core_distance(dists, d=2)
    assert core_dists.shape == (3, 1)


def test_compute_cluster_core_distance_large(sample_data_large):
    dists = compute_pair_to_pair_dists(sample_data_large, metric="euclidean")
    core_dists = compute_cluster_core_distance(dists, d=2)
    assert core_dists.shape == (1000, 1)


def test_compute_mutual_reach_dists(sample_data):
    dists = compute_pair_to_pair_dists(sample_data, metric="euclidean")
    mutual_reach_dists = compute_mutual_reach_dists(dists, d=2, is_symmetric=True)
    assert mutual_reach_dists.shape == (3, 3)


def test_fn_density_sparseness():
    cls_inds = np.array([1, 2, 3])
    dists = np.array([[np.inf, 1.0, 2.0], [1.0, np.inf, 1.0], [2.0, 1.0, np.inf]])
    dsc, internal_node_inds = fn_density_sparseness(cls_inds, dists, d=2)
    assert np.isclose(dsc, 0, rtol=1e-8)
    np.testing.assert_array_equal(internal_node_inds, np.array([]))

    assert cls_inds.size == 3
    assert isinstance(dsc, float)
    assert isinstance(internal_node_inds, np.ndarray)


def test_fn_density_sparseness_larger():
    cls_inds = np.array([1, 2, 3, 4, 5])
    dists = np.array(
        [
            [np.inf, 1.0, 2.0, 3.0, 4.0],
            [1.0, np.inf, 1.0, 2.0, 3.0],
            [2.0, 1.0, np.inf, 1.0, 2.0],
            [3.0, 2.0, 1.0, np.inf, 1.0],
            [4.0, 3.0, 2.0, 1.0, np.inf],
        ]
    )

    dsc, internal_node_inds = fn_density_sparseness(cls_inds, dists, d=2)
    assert np.isclose(dsc, 1.2649110640675099, rtol=1e-8)
    np.testing.assert_array_equal(internal_node_inds, np.array([2, 3, 4]))

    assert cls_inds.size > 3
    assert isinstance(dsc, float)
    assert isinstance(internal_node_inds, np.ndarray)
    assert internal_node_inds.size > 0


def test_fn_density_separation(sample_data):
    dists = compute_pair_to_pair_dists(sample_data, metric="euclidean")
    cls_i, cls_j, dspc_ij = fn_density_separation(0, 1, dists, d=2)
    assert isinstance(cls_i, int)
    assert isinstance(cls_j, int)
    assert isinstance(dspc_ij, float)


def test_check_duplicated_samples_unique_value():
    X = np.array([1])
    result = _check_duplicated_samples(X)
    assert result is None


def test_check_duplicated_samples():
    X = np.array([[1, 2], [3, 4], [1, 2]])
    with pytest.raises(ValueError, match="Duplicated samples have been found in X."):
        _check_duplicated_samples(X)


def test_dbcv_one_dimension():
    X, _ = datasets.make_moons(n_samples=30, noise=0.05, random_state=1782)
    X = X.flatten()
    y = np.zeros((60,))

    actual_dbcv_score = dbcv_score(X, y)
    expected_dbcv_score = 0.9999999999998341

    assert actual_dbcv_score == pytest.approx(expected_dbcv_score, rel=1e-6)


def test_dbcv_value_error_on_dimension_mismatch():
    X = np.random.rand(10, 3)
    y = np.random.randint(0, 2, size=(15,))

    with pytest.raises(ValueError, match=r"Mismatch in .* and .* dimensions."):
        dbcv_score(X, y)


def test_dbcv_noise_id_equals_all_y_values():
    X = np.random.rand(10, 3)
    y = np.zeros(10, dtype=int)

    result = dbcv_score(X, y, noise_id=0)

    assert result == 0.0


def test_dbcv_kmeans_dbscan():
    X, _ = datasets.make_moons(n_samples=100, noise=0.05, random_state=1782)

    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=2, algorithm="lloyd", n_init=10)
    kmeans_labels = kmeans.fit_predict(X)

    from sklearn.cluster import DBSCAN

    dbscanner = DBSCAN(algorithm="ball_tree")
    dbscan_labels = dbscanner.fit_predict(X)

    actual_kmeans_score = dbcv_score(X, kmeans_labels)
    actual_dbscan_score = dbcv_score(X, dbscan_labels)

    expected_dbscan_score = 0.9999999999997616

    # Kmeans_score is randomly between -0.2 and -0.5 so we can't give a unique value
    assert -0.5 <= actual_kmeans_score <= 0.2
    assert actual_dbscan_score == pytest.approx(expected_dbscan_score, rel=1e-6)
