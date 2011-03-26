import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

from scikits.learn import neighbors, datasets

# load and shuffle iris dataset
iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_neighbors_1D():
    """
    Nearest Neighbors in a line.

    Samples are are set of n two-category equally spaced points.
    """
    # some constants
    n = 6
    X = [[x] for x in range(0, n)]
    Y = [0]*(n/2) + [1]*(n/2)

    for s in ('auto', 'ball_tree', 'brute', 'inplace'):
        # n_neighbors = 1
        knn = neighbors.NeighborsClassifier(n_neighbors=1, algorithm=s)
        knn.fit(X, Y)
        test = [[i + 0.01] for i in range(0, n/2)] + \
               [[i - 0.01] for i in range(n/2, n)]
        assert_array_equal(knn.predict(test), [0]*3 + [1]*3)

        # n_neighbors = 2
        knn = neighbors.NeighborsClassifier(n_neighbors=2, algorithm=s)
        knn.fit(X, Y)
        assert_array_equal(knn.predict(test), [0]*4 + [1]*2)

        # n_neighbors = 3
        knn = neighbors.NeighborsClassifier(n_neighbors=3, algorithm=s)
        knn.fit(X, Y)
        assert_array_equal(knn.predict([[i +0.01] for i in range(0, n/2)]),
                            [0 for i in range(n/2)])
        assert_array_equal(knn.predict([[i-0.01] for i in range(n/2, n)]),
                            [1 for i in range(n/2)])


def test_neighbors_high_dimension():
    """ Nearest Neighbors on high-dimensional data.
    """
    # some constants
    n = 20
    p = 40
    X = 2*np.random.random(size=(n, p)) - 1
    Y = ((X**2).sum(axis=1) < .25).astype(np.int)

    for s in ('auto', 'ball_tree', 'brute', 'inplace'):
        knn = neighbors.NeighborsClassifier(n_neighbors=1, algorithm=s)
        knn.fit(X, Y)
        for i, (x, y) in enumerate(zip(X[:10], Y[:10])):
            epsilon = 1e-5*(2*np.random.random(size=p)-1)
            assert_array_equal(knn.predict(x+epsilon), y)
            dist, idxs = knn.kneighbors(x+epsilon, n_neighbors=1)


def test_neighbors_iris():
    """
    Sanity checks on the iris dataset

    Puts three points of each label in the plane and performs a
    nearest neighbor query on points near the decision boundary.
    """

    for s in ('auto', 'ball_tree', 'brute', 'inplace'):
        clf = neighbors.NeighborsClassifier()
        clf.fit(iris.data, iris.target, n_neighbors=1, algorithm=s)
        assert_array_equal(clf.predict(iris.data), iris.target)

        clf.fit(iris.data, iris.target, n_neighbors=9, algorithm=s)
        assert np.mean(clf.predict(iris.data)== iris.target) > 0.95

        for m in ('barycenter', 'mean'):
            rgs = neighbors.NeighborsRegressor()
            rgs.fit(iris.data, iris.target, mode=m, algorithm=s)
            assert np.mean(
                rgs.predict(iris.data).round() == iris.target) > 0.95


def test_kneighbors_graph():
    """
    Test kneighbors_graph to build the k-Nearest Neighbor graph.
    """
    X = [[0, 1], [1.01, 1.], [2, 0]]

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
    assert_array_almost_equal(
        A.todense(),
        [[ 0.        ,  1.5049745 , -0.5049745 ],
        [ 0.596     ,  0.        ,  0.404     ],
        [-0.98019802,  1.98019802,  0.        ]])

    # n_neighbors = 3
    A = neighbors.kneighbors_graph(X, 3, mode='connectivity')
    assert_array_almost_equal(
        A.todense(),
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]])


if __name__ == '__main__':
    import nose
    nose.runmodule()
