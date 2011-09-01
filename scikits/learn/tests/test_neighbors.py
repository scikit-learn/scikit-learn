import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from scipy.sparse import (bsr_matrix, coo_matrix, csc_matrix, csr_matrix,
                          dok_matrix, lil_matrix)

from scikits.learn import neighbors, datasets, ball_tree

# load and shuffle iris dataset
iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

SPARSE_TYPES = (bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dok_matrix,
                lil_matrix)
SPARSE_OR_DENSE = SPARSE_TYPES + (np.asarray,)


def test_neighbors_1D():
    """
    Nearest Neighbors in a line.

    Samples are are set of n two-category equally spaced points.
    """
    # some constants
    n = 6
    X = [[x] for x in range(0, n)]
    Y = [0] * (n / 2) + [1] * (n / 2)

    for s in ('auto', 'ball_tree', 'brute', 'inplace'):
        # n_neighbors = 1
        knn = neighbors.NeighborsClassifier(n_neighbors=1, algorithm=s)
        knn.fit(X, Y)
        test = [[i + 0.01] for i in range(0, n / 2)] + \
               [[i - 0.01] for i in range(n / 2, n)]
        assert_array_equal(knn.predict(test), [0] * 3 + [1] * 3)

        # n_neighbors = 2
        knn = neighbors.NeighborsClassifier(n_neighbors=2, algorithm=s)
        knn.fit(X, Y)
        assert_array_equal(knn.predict(test), [0] * 4 + [1] * 2)

        # n_neighbors = 3
        knn = neighbors.NeighborsClassifier(n_neighbors=3, algorithm=s)
        knn.fit(X, Y)
        assert_array_equal(knn.predict([[i + 0.01] for i in range(0, n / 2)]),
                            [0 for i in range(n / 2)])
        assert_array_equal(knn.predict([[i - 0.01] for i in range(n / 2, n)]),
                            [1 for i in range(n / 2)])


def test_neighbors_high_dimension():
    """ Nearest Neighbors on high-dimensional data.
    """
    # some constants
    n = 20
    p = 40
    X = 2 * np.random.random(size=(n, p)) - 1
    Y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    for s in ('auto', 'ball_tree', 'brute', 'inplace'):
        knn = neighbors.NeighborsClassifier(n_neighbors=1, algorithm=s)
        knn.fit(X, Y)
        for i, (x, y) in enumerate(zip(X[:10], Y[:10])):
            epsilon = 1e-5 * (2 * np.random.random(size=p) - 1)
            assert_array_equal(knn.predict(x + epsilon), y)
            dist, idxs = knn.kneighbors(x + epsilon, n_neighbors=1)


def test_neighbors_sparse_classification():
    """Test k-NN classifier on sparse matrices"""

    # Like the above, but with various types of sparse matrices
    n = 10
    p = 30
    X = 2 * np.random.random(size=(n, p)) - 1
    Y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    SPARSE_TYPES = (bsr_matrix, coo_matrix, csc_matrix, csr_matrix,
                    dok_matrix, lil_matrix)
    for sparsemat in SPARSE_TYPES:
        # 'ball_tree' option should be overridden automatically
        knn = neighbors.NeighborsClassifier(n_neighbors=1,
                                            algorithm='ball_tree')
        knn.fit(sparsemat(X), Y)

        for i, (x, y) in enumerate(zip(X[:5], Y[:5])):
            epsilon = 1e-5 * (2 * np.random.random(size=p) - 1)
            for sparsev in SPARSE_TYPES + (np.asarray,):
                x_eps = sparsev(np.atleast_2d(x) + epsilon)
                assert_array_equal(knn.predict(x_eps), y)
                dist, idxs = knn.kneighbors(x_eps, n_neighbors=1)


def test_neighbors_iris():
    """
    Sanity checks on the iris dataset

    Puts three points of each label in the plane and performs a
    nearest neighbor query on points near the decision boundary.
    """

    for s in ('auto', 'ball_tree', 'brute', 'inplace'):
        clf = neighbors.NeighborsClassifier(n_neighbors=1, algorithm=s)
        clf.fit(iris.data, iris.target)
        assert_array_equal(clf.predict(iris.data), iris.target)

        clf.set_params(n_neighbors=9, algorithm=s)
        clf.fit(iris.data, iris.target)
        assert np.mean(clf.predict(iris.data) == iris.target) > 0.95

        for m in ('barycenter', 'mean'):
            rgs = neighbors.NeighborsRegressor(mode=m)
            rgs.fit(iris.data, iris.target)
            assert np.mean(
                rgs.predict(iris.data).round() == iris.target) > 0.95


# Test disabled because sparse regression is not yet supported.
# Remove the leading _ to enable it.
def _test_neighbors_sparse_regression():
    """Test k-NN regression on sparse matrices

    Repeats part of the iris test.
    """

    for sparse1 in SPARSE_TYPES:
        for m in ('barycenter', 'mean'):
            rgs = neighbors.NeighborsRegressor(mode=m)
            rgs.fit(sparse1(iris.data), iris.target)
            for sparse2 in SPARSE_OR_DENSE:
                data2 = sparse2(iris.data)
                assert (np.mean(rgs.predict(data2).round() == iris.target)
                        > 0.95)


def test_kneighbors_graph():
    """
    Test kneighbors_graph to build the k-Nearest Neighbor graph.
    """
    X = np.array([[0, 1], [1.01, 1.], [2, 0]])

    # n_neighbors = 1
    A = neighbors.kneighbors_graph(X, 1, mode='connectivity')
    assert_array_equal(A.todense(), np.eye(A.shape[0]))

    A = neighbors.kneighbors_graph(X, 1, mode='distance')
    assert_array_almost_equal(
        A.todense(),
        [[ 0.        ,  1.01      ,  0.        ],
         [ 1.01      ,  0.        ,  0.        ],
         [ 0.        ,  1.40716026,  0.        ]])

    A = neighbors.kneighbors_graph(X, 1, mode='barycenter')
    assert_array_almost_equal(
        A.todense(),
        [[ 0.,  1.,  0.],
         [ 1.,  0.,  0.],
         [ 0.,  1.,  0.]])

    # n_neighbors = 2
    A = neighbors.kneighbors_graph(X, 2, mode='connectivity')
    assert_array_equal(
        A.todense(),
        [[ 1.,  1.,  0.],
         [ 1.,  1.,  0.],
         [ 0.,  1.,  1.]])

    A = neighbors.kneighbors_graph(X, 2, mode='distance')
    assert_array_almost_equal(
        A.todense(),
        [[ 0.        ,  1.01      ,  2.23606798],
         [ 1.01      ,  0.        ,  1.40716026],
         [ 2.23606798,  1.40716026,  0.        ]])

    A = neighbors.kneighbors_graph(X, 2, mode='barycenter')
    # check that columns sum to one
    assert_array_almost_equal(np.sum(A.todense(), 1), np.ones((3, 1)))
    pred = np.dot(A.todense(), X)
    assert np.linalg.norm(pred - X) / X.shape[0] < 1

    # n_neighbors = 3
    A = neighbors.kneighbors_graph(X, 3, mode='connectivity')
    assert_array_almost_equal(
        A.todense(),
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]])


def test_radius_neighbors_graph():
    """
    Test radius_neighbors_graph to build the Nearest Neighbor graph.
    """
    X = np.array([[0, 1], [1.01, 1.], [2, 0]])

    A = neighbors.radius_neighbors_graph(X, 1.5, mode='connectivity')
    assert_array_equal(
        A.todense(),
        [[ 1.,  1.,  0.],
         [ 1.,  1.,  1.],
         [ 0.,  1.,  1.]])

    A = neighbors.radius_neighbors_graph(X, 1.5, mode='distance')
    assert_array_almost_equal(
        A.todense(),
        [[ 0.        ,  1.01      ,  0.        ],
         [ 1.01      ,  0.        ,  1.40716026],
         [ 0.        ,  1.40716026,  0.        ]])


def test_kneighbors_iris():

    # make sure reconstruction error is kept small using a real datasets
    # note that we choose neighbors < n_dim and n_neighbors > n_dim
    for i in range(1, 8):
        for data in (iris.data, neighbors.BallTree(iris.data)):
            # check for both input as numpy array and as BallTree
            A = neighbors.kneighbors_graph(data, i, mode='barycenter')
            if hasattr(data, 'query'):
                data = data.data
            pred_data = np.dot(A.todense(), data)
            assert np.linalg.norm(pred_data - data) / data.shape[0] < 0.1


def test_ball_tree_query_radius(n_samples=100, n_features=10):
    X = 2 * np.random.random(size=(n_samples, n_features)) - 1
    query_pt = np.zeros(n_features, dtype=float)

    eps = 1E-15  # roundoff error can cause test to fail
    bt = ball_tree.BallTree(X)
    rad = np.sqrt(((X - query_pt) ** 2).sum(1))

    for r in np.linspace(rad[0], rad[-1], 100):
        ind = bt.query_radius(query_pt, r + eps)[0]
        i = np.where(rad <= r + eps)[0]

        ind.sort()
        i.sort()

        assert np.all(i == ind)


def test_ball_tree_query_radius_distance(n_samples=100, n_features=10):
    X = 2 * np.random.random(size=(n_samples, n_features)) - 1
    query_pt = np.zeros(n_features, dtype=float)

    eps = 1E-15  # roundoff error can cause test to fail
    bt = ball_tree.BallTree(X)
    rad = np.sqrt( ((X - query_pt) ** 2).sum(1) )

    for r in np.linspace(rad[0], rad[-1], 100):
        ind, dist = bt.query_radius(query_pt, r + eps, return_distance=True)

        ind = ind[0]
        dist = dist[0]

        d = np.sqrt(((query_pt - X[ind]) ** 2).sum(1))

        assert_array_almost_equal(d, dist)


def test_ball_tree_pickle():
    import pickle
    X = np.random.random(size=(10, 3))
    bt1 = ball_tree.BallTree(X, 1)
    ind1, dist1 = bt1.query(X)
    for protocol in (0, 1, 2):
        s = pickle.dumps(bt1, protocol=protocol)
        bt2 = pickle.loads(s)
        ind2, dist2 = bt2.query(X)
        assert np.all(ind1 == ind2)
        assert_array_almost_equal(dist1, dist2)


if __name__ == '__main__':
    import nose
    nose.runmodule()
