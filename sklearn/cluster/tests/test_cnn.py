"""Tests for CNN clustering
"""

import pickle

import numpy as np

from scipy.spatial import distance
from scipy import sparse

import pytest

from sklearn.utils._testing import assert_array_equal
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import CNN
from sklearn.cluster import cnn
from sklearn.cluster.tests.common import generate_clustered_data
from sklearn.metrics.pairwise import pairwise_distances


n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)


def test_cnn_similarity():
    # Tests the CNN algorithm with a similarity array.
    # Parameters chosen specifically for this task.
    eps = 0.15
    min_samples = 5
    # Compute similarities
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)
    # Compute CNN
    labels = cnn(D, metric="precomputed", eps=eps,
                 min_samples=min_samples)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - (1 if -1 in labels else 0)

    assert n_clusters_1 == n_clusters

    cobj = CNN(metric="precomputed", eps=eps, min_samples=min_samples)
    labels = cobj.fit(D).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_cnn_feature():
    # Tests the CNN algorithm with a feature vector array.
    # Parameters chosen specifically for this task.
    # Different eps to other test, because distance is not normalised.
    eps = 0.8
    min_samples = 5
    metric = 'euclidean'
    # Compute CNN
    # parameters chosen for task
    labels = cnn(X, metric=metric, eps=eps,
                 min_samples=min_samples)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    cobj = CNN(metric=metric, eps=eps, min_samples=min_samples)
    labels = cobj.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_cnn_sparse():
    labels_sparse = cnn(sparse.lil_matrix(X), eps=.8, min_samples=5)
    labels_dense = cnn(X, eps=.8, min_samples=5)
    assert_array_equal(labels_dense, labels_sparse)


@pytest.mark.parametrize('include_self', [False, True])
def test_cnn_sparse_precomputed(include_self):
    D = pairwise_distances(X)
    nn = NearestNeighbors(radius=.9).fit(X)
    X_ = X if include_self else None
    D_sparse = nn.radius_neighbors_graph(X=X_, mode='distance')
    # Ensure it is sparse not merely on diagonals:
    assert D_sparse.nnz < D.shape[0] * (D.shape[0] - 1)
    labels_sparse = cnn(D_sparse,
                        eps=.8,
                        min_samples=5,
                        metric='precomputed')
    labels_dense = cnn(D, eps=.8, min_samples=5,
                       metric='precomputed')
    assert_array_equal(labels_dense, labels_sparse)


def test_cnn_sparse_precomputed_different_eps():
    # test that precomputed neighbors graph is filtered if computed with
    # a radius larger than CNN's eps.
    lower_eps = 0.2
    nn = NearestNeighbors(radius=lower_eps).fit(X)
    D_sparse = nn.radius_neighbors_graph(X, mode='distance')
    cnn_lower = cnn(D_sparse, eps=lower_eps, metric='precomputed')

    higher_eps = lower_eps + 0.7
    nn = NearestNeighbors(radius=higher_eps).fit(X)
    D_sparse = nn.radius_neighbors_graph(X, mode='distance')
    cnn_higher = cnn(D_sparse, eps=lower_eps, metric='precomputed')

    assert_array_equal(cnn_lower, cnn_higher)


@pytest.mark.parametrize('use_sparse', [True, False])
@pytest.mark.parametrize('metric', ['precomputed', 'minkowski'])
def test_cnn_input_not_modified(use_sparse, metric):
    # test that the input is not modified by cnn
    X = np.random.RandomState(0).rand(10, 10)
    X = sparse.csr_matrix(X) if use_sparse else X
    X_copy = X.copy()
    cnn(X, metric=metric)

    if use_sparse:
        assert_array_equal(X.toarray(), X_copy.toarray())
    else:
        assert_array_equal(X, X_copy)


def test_cnn_callable():
    # Tests the CNN algorithm with a callable metric.
    # Parameters chosen specifically for this task.
    # Different eps to other test, because distance is not normalised.
    eps = 0.8
    min_samples = 5
    # metric is the function reference, not the string key.
    metric = distance.euclidean
    # Compute CNN
    # parameters chosen for task
    labels = cnn(X, metric=metric, eps=eps,
                 min_samples=min_samples,
                 algorithm='ball_tree')

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    cobj = CNN(metric=metric, eps=eps, min_samples=min_samples,
               algorithm='ball_tree')
    labels = cobj.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters


def test_cnn_metric_params():
    # Tests that CNN works with the metrics_params argument.
    eps = 0.8
    min_samples = 5
    p = 1

    # Compute DBSCAN with metric_params arg
    cobj = CNN(metric='minkowski', metric_params={'p': p}, eps=eps,
               min_samples=min_samples, algorithm='ball_tree').fit(X)
    labels_1 = cobj.labels_

    # Test that sample labels are the same as passing Minkowski 'p' directly
    cobj = CNN(metric='minkowski', eps=eps, min_samples=min_samples,
               algorithm='ball_tree', p=p).fit(X)
    labels_2 = cobj.labels_

    assert_array_equal(labels_1, labels_2)

    # Minkowski with p=1 should be equivalent to Manhattan distance
    cobj = CNN(metric='manhattan', eps=eps, min_samples=min_samples,
               algorithm='ball_tree').fit(X)
    labels_3 = cobj.labels_

    assert_array_equal(labels_1, labels_3)


def test_cnn_balltree():
    # Tests the CNN algorithm with balltree for neighbor calculation.
    eps = 0.8
    min_samples = 5

    D = pairwise_distances(X)
    labels = cnn(D, metric="precomputed", eps=eps,
                 min_samples=min_samples)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_1 == n_clusters

    cobj = CNN(p=2.0, eps=eps, min_samples=min_samples, algorithm='ball_tree')
    labels = cobj.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_2 == n_clusters

    cobj = CNN(p=2.0, eps=eps, min_samples=min_samples, algorithm='kd_tree')
    labels = cobj.fit(X).labels_

    n_clusters_3 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_3 == n_clusters

    cobj = CNN(p=1.0, eps=eps, min_samples=min_samples, algorithm='ball_tree')
    labels = cobj.fit(X).labels_

    n_clusters_4 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_4 == n_clusters

    cobj = CNN(leaf_size=20, eps=eps, min_samples=min_samples,
               algorithm='ball_tree')
    labels = cobj.fit(X).labels_

    n_clusters_5 = len(set(labels)) - int(-1 in labels)
    assert n_clusters_5 == n_clusters


def test_input_validation():
    # CNN.fit should accept a list of lists.
    X = [[1., 2.], [3., 4.]]
    CNN().fit(X)             # must not raise exception


@pytest.mark.parametrize(
    "args",
    [{'eps': -1.0}, {'algorithm': 'blah'}, {'metric': 'blah'},
     {'leaf_size': -1}, {'p': -1}]
)
def test_cnn_badargs(args):
    # Test bad argument values: these should all raise ValueErrors
    with pytest.raises(ValueError):
        cnn(X, **args)


def test_pickle():
    obj = CNN()
    s = pickle.dumps(obj)
    assert type(pickle.loads(s)) == obj.__class__


def test_boundaries():
    # ensure min_samples is inclusive of core point
    core = np.where(cnn([[0], [1]], eps=2, min_samples=0) >= 0)[0]
    assert 0 in core
    # ensure eps is inclusive of circumference
    core = np.where(cnn([[0], [1], [1]], eps=1, min_samples=0) >= 0)[0]
    assert 0 in core
    core = np.where(cnn([[0], [1], [1]], eps=.99, min_samples=0) >= 0)[0]
    assert 0 not in core


@pytest.mark.skip(reason="Sample weights have no effect on the clustering yet")
def test_weighted_cnn():
    # ensure sample_weight is validated
    with pytest.raises(ValueError):
        cnn([[0], [1]], sample_weight=[2])
    with pytest.raises(ValueError):
        cnn([[0], [1]], sample_weight=[2, 3, 4])

    # ensure sample_weight has an effect
    assert_array_equal([], cnn([[0], [1]], sample_weight=None,
                               min_samples=6)[0])
    assert_array_equal([], cnn([[0], [1]], sample_weight=[5, 5],
                               min_samples=6)[0])
    assert_array_equal([0], cnn([[0], [1]], sample_weight=[6, 5],
                                min_samples=6)[0])
    assert_array_equal([0, 1], cnn([[0], [1]], sample_weight=[6, 6],
                                   min_samples=6)[0])

    # points within eps of each other:
    assert_array_equal([0, 1], cnn([[0], [1]], eps=1.5,
                                   sample_weight=[5, 1], min_samples=6)[0])
    # and effect of non-positive and non-integer sample_weight:
    assert_array_equal([], cnn([[0], [1]], sample_weight=[5, 0],
                               eps=1.5, min_samples=6)[0])
    assert_array_equal([0, 1], cnn([[0], [1]], sample_weight=[5.9, 0.1],
                                   eps=1.5, min_samples=6)[0])
    assert_array_equal([0, 1], cnn([[0], [1]], sample_weight=[6, 0],
                                   eps=1.5, min_samples=6)[0])
    assert_array_equal([], cnn([[0], [1]], sample_weight=[6, -1],
                               eps=1.5, min_samples=6)[0])

    # for non-negative sample_weight, cores should be identical to repetition
    rng = np.random.RandomState(42)
    sample_weight = rng.randint(0, 5, X.shape[0])
    label1 = cnn(X, sample_weight=sample_weight)
    assert len(label1) == len(X)

    X_repeated = np.repeat(X, sample_weight, axis=0)
    label_repeated = cnn(X_repeated)
    core_repeated_mask = np.zeros(X_repeated.shape[0], dtype=bool)
    core_repeated_mask[np.where(label_repeated >= 0)[0]] = True
    core_mask = np.zeros(X.shape[0], dtype=bool)
    core_mask[np.where(label1 >= 0)[0]] = True
    assert_array_equal(np.repeat(core_mask, sample_weight), core_repeated_mask)

    # sample_weight should work with precomputed distance matrix
    D = pairwise_distances(X)
    core3, label3 = cnn(D, sample_weight=sample_weight,
                        metric='precomputed')
    assert_array_equal(label1, label3)

    # sample_weight should work with estimator
    est = CNN().fit(X, sample_weight=sample_weight)
    label4 = est.labels_
    assert_array_equal(label1, label4)

    est = CNN()
    label5 = est.fit_predict(X, sample_weight=sample_weight)

    assert_array_equal(label1, label5)
    assert_array_equal(label1, est.labels_)


@pytest.mark.parametrize('algorithm', ['brute', 'kd_tree', 'ball_tree'])
def test_cnn_core_samples_toy_1(algorithm):
    X = [[0], [2], [3], [4], [6], [8], [10]]

    # Within eps = 1, only points at 2, 3, and 4
    # are neighbours. Valid clusters need to have more than one
    # members, so all other points are isolated and considered
    # noise.
    labels = cnn(X, algorithm=algorithm, eps=1, min_samples=0)
    assert_array_equal(labels, [-1, 0, 0, 0, -1, -1, -1])

    # With eps=1 and min_samples=1 the 3 samples from the
    # denser area are no core samples anymore (2 and 4 share 3 as
    # common neighbour but are not neighbours of each other)
    labels = cnn(X, algorithm=algorithm, eps=1, min_samples=1)
    assert_array_equal(labels, [-1, -1, -1, -1, -1, -1, -1])


@pytest.mark.parametrize('algorithm', ['brute', 'kd_tree', 'ball_tree'])
def test_cnn_core_samples_toy_2(algorithm):
    X = [[0, 0],       # 0
         [1, 1],       # 1
         [1, 0],       # 2
         [0, -1],      # 3
         [0.5, -0.5],  # 4
         [2,  1.5],    # 5
         [2.5, -0.5],  # 6
         [4, 2],       # 7
         [4.5, 2.5],   # 8
         [5, -1],      # 9
         [5.5, -0.5],  # 10
         [5.5, -1.5]]  # 11

    labels = cnn(X, algorithm=algorithm, eps=1.5, min_samples=0)
    assert_array_equal(labels, [0, 0, 0, 0, 0, 0, -1, 1, 1, 2, 2, 2])

    labels = cnn(X, algorithm=algorithm, eps=1.5, min_samples=1)
    assert_array_equal(labels, [0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, 1])

    labels = cnn(X, algorithm=algorithm, eps=1.5, min_samples=2)
    assert_array_equal(labels, [0, -1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1])

    labels = cnn(X, algorithm=algorithm, eps=1.5, min_samples=3)
    assert_array_equal(labels, [0, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1])

    labels = cnn(X, algorithm=algorithm, eps=1.5, min_samples=4)
    assert_array_equal(labels,
                       [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])


def test_cnn_precomputed_metric_with_degenerate_input_arrays():
    # see https://github.com/scikit-learn/scikit-learn/issues/4641 for
    # more details
    X = np.eye(10)
    labels = CNN(eps=0.5, metric='precomputed').fit(X).labels_
    assert len(set(labels)) == 1

    X = np.zeros((10, 10))
    labels = CNN(eps=0.5, metric='precomputed').fit(X).labels_
    assert len(set(labels)) == 1


def test_cnn_precomputed_metric_with_initial_rows_zero():
    # sample matrix with initial two row all zero
    ar = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0],
        [0.0, 0.0, 0.1, 0.1, 0.0, 0.0, 0.3],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1],
        [0.0, 0.0, 0.0, 0.0, 0.3, 0.1, 0.0]
    ])
    matrix = sparse.csr_matrix(ar)
    labels = CNN(eps=0.2, metric='precomputed',
                 min_samples=0).fit(matrix).labels_
    assert_array_equal(labels, [-1, -1,  0,  0,  0,  1,  1])
