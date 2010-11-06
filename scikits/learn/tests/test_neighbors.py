
from numpy.testing import assert_array_equal

from .. import neighbors

def test_neighbors_1D():
    """
    Nearest Neighbors in a line.

    Samples are are set of n two-category equally spaced points.
    """
    # some constants
    n = 6
    n_2 = n/2
    X = [[x] for x in range(0, n)]
    Y  = [0]*n_2 + [1]*n_2

    # k = 1
    knn = neighbors.Neighbors(k=1)
    knn.fit(X, Y)
    test = [[i +0.01] for i in range(0,n_2)] + [[i-0.01] for i in range(n_2, n)]
    assert_array_equal( knn.predict(test), [0, 0, 0, 1, 1, 1])
    # same as before, but using predict() instead of Neighbors object

    # k = 3
    knn = neighbors.Neighbors(k=3)
    knn.fit(X, Y)
    assert_array_equal( knn.predict([ [i +0.01] for i in range(0,n_2)]),
                        [0 for i in range(n_2)])
    assert_array_equal( knn.predict([ [i-0.01] for i in range(n_2, n)]),
                        [1 for i in range(n_2)])


def test_neighbors_2D():
    """
    Nearest Neighbor in the plane.

    Puts three points of each label in the plane and performs a
    nearest neighbor query on points near the decision boundary.
    """
    X = (
        (0, 1), (1, 1), (1, 0), # label 0
        (-1,0), (-1,-1),(0,-1)) # label 1
    n_2 = len(X)/2
    Y = [0]*n_2 + [1]*n_2
    knn = neighbors.Neighbors()
    knn.fit(X, Y)

    prediction =  knn.predict([[0, .1], [0, -.1], [.1, 0], [-.1, 0]])
    assert_array_equal(prediction, [0, 1, 0, 1])

