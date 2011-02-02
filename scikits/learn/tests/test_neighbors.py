import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
                          assert_equal

from scikits.learn import neighbors


def test_neighbors_1D():
    """
    Nearest Neighbors in a line.

    Samples are are set of n two-category equally spaced points.
    """
    # some constants
    n = 6
    X = [[x] for x in range(0, n)]
    Y = [0]*(n/2) + [1]*(n/2)

    # n_neighbors = 1
    knn = neighbors.Neighbors(n_neighbors=1)
    knn.fit(X, Y)
    test = [[i + 0.01] for i in range(0, n/2)] + \
           [[i - 0.01] for i in range(n/2, n)]
    assert_array_equal(knn.predict(test), [0]*3 + [1]*3)

    # n_neighbors = 2
    knn = neighbors.Neighbors(n_neighbors=2)
    knn.fit(X, Y)
    assert_array_equal(knn.predict(test), [0]*4 + [1]*2)


    # n_neighbors = 3
    knn = neighbors.Neighbors(n_neighbors=3)
    knn.fit(X, Y)
    assert_array_equal(knn.predict([[i +0.01] for i in range(0, n/2)]),
                        [0 for i in range(n/2)])
    assert_array_equal(knn.predict([[i-0.01] for i in range(n/2, n)]),
                        [1 for i in range(n/2)])


def test_neighbors_2D():
    """
    Nearest Neighbor in the plane.

    Puts three points of each label in the plane and performs a
    nearest neighbor query on points near the decision boundary.
    """
    X = (
        (0, 1), (1, 1), (1, 0), # label 0
        (-1, 0), (-1, -1), (0, -1)) # label 1
    n_2 = len(X)/2
    Y = [0]*n_2 + [1]*n_2
    knn = neighbors.Neighbors()
    knn.fit(X, Y)

    prediction = knn.predict([[0, .1], [0, -.1], [.1, 0], [-.1, 0]])
    assert_array_equal(prediction, [0, 1, 0, 1])


def test_neighbors_barycenter():
    """
    NeighborsBarycenter for regression using k-NN
    """
    X = [[0], [1], [2], [3]]
    y = [0, 0, 1, 1]
    neigh = neighbors.NeighborsBarycenter(n_neighbors=2)
    neigh.fit(X, y)
    assert_equal(neigh.predict([[1.5]]), 0.5)


def test_kneighbors_graph():
    """
    Test kneighbors_graph to build the k-Nearest Neighbor graph.
    """
    X = [[0, 0], [1.01, 0], [2, 0]]

    # n_neighbors = 1
    A = neighbors.kneighbors_graph(X, 1, mode='connectivity')
    assert_array_equal(A.todense(), np.eye(A.shape[0]))

    A = neighbors.kneighbors_graph(X, 1, mode='distance')
    assert_array_equal(
        A.todense(),
        [[ 0.  ,  1.01,  0.  ],
         [ 0.  ,  0.  ,  0.99],
         [ 0.  ,  0.99,  0.  ]])

    A = neighbors.kneighbors_graph(X, 1, mode='barycenter')
    assert_array_equal(
        A.todense(),
        [[ 0.,  1.,  0.],
        [ 0.,  0.,  1.],
        [ 0.,  1.,  0.]])

    # n_neigbors = 2
    A = neighbors.kneighbors_graph(X, 2, mode='connectivity')
    assert_array_equal(
        A.todense(),
        [[ 1.,  1.,  0.],
         [ 0.,  1.,  1.],
         [ 0.,  1.,  1.]])

    A = neighbors.kneighbors_graph(X, 2, mode='distance')
    assert_array_almost_equal(
        A.todense(),
        [[ 0.  ,  1.01,  2.  ],
        [ 1.01,  0.  ,  0.99],
        [ 2.  ,  0.99,  0.  ]])

    A = neighbors.kneighbors_graph(X, 2, mode='barycenter')
    # check that columns sum to one
    assert_array_almost_equal(np.sum(A.todense(), 1), np.ones((3, 1)))
    assert_array_almost_equal(
        A.todense(),
        [[ 0.        ,  2.02018645, -1.02018645],
        [ 0.49500001,  0.        ,  0.50499999],
        [-0.98018357,  1.98018357,  0.        ]])
    # check that we can reconstruct X from A
    assert_array_almost_equal(
        X, A.dot(X), decimal=3)

    # n_neighbors = 3
    A = neighbors.kneighbors_graph(X, 3, mode='connectivity')
    assert_array_almost_equal(
        A.todense(),
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]])


if __name__ == '__main__':
    import nose
    nose.runmodule()
