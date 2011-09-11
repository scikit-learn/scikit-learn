import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from scipy.sparse import (bsr_matrix, coo_matrix, csc_matrix, csr_matrix,
                          dok_matrix, lil_matrix)

from sklearn import neighbors, datasets

# load and shuffle iris dataset
iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

SPARSE_TYPES = (bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dok_matrix,
                lil_matrix)
SPARSE_OR_DENSE = SPARSE_TYPES + (np.asarray,)

ALGORITHMS = ('ball_tree', 'brute', 'kd_tree', 'auto')


def test_unsupervised_kneighbors(n_samples=20, n_features=5,
                                 n_query_pts=2, n_neighbors=5,
                                 random_state=0):
    """
    Samples are are set of n two-category equally spaced points.
    """
    rng = np.random.RandomState(random_state)

    X = rng.rand(n_samples, n_features)

    test = rng.rand(n_query_pts, n_features)

    results = []

    for algorithm in ALGORITHMS:
        neigh = neighbors.NearestNeighbors(n_neighbors=n_neighbors,
                                           algorithm=algorithm)
        neigh.fit(X)

        results.append(neigh.kneighbors(test, return_distance=True))

    for i in range(len(results) - 1):
        assert_array_almost_equal(results[i][0], results[i + 1][0])
        assert_array_almost_equal(results[i][1], results[i + 1][1])


def test_unsupervised_radius_neighbors(n_samples=20, n_features=5,
                                       n_query_pts=2, radius=0.5,
                                       random_state=0):
    """
    Samples are are set of n two-category equally spaced points.
    """
    rng = np.random.RandomState(random_state)

    X = rng.rand(n_samples, n_features)

    test = rng.rand(n_query_pts, n_features)

    results = []

    for algorithm in ALGORITHMS:
        neigh = neighbors.NearestNeighbors(radius=radius,
                                           algorithm=algorithm)
        neigh.fit(X)

        # sort the results: this is not done automatically for
        # radius searches
        dist, ind = neigh.radius_neighbors(test, return_distance=True)
        for (d, i) in zip(dist, ind):
            j = d.argsort()
            d[:] = d[j]
            i[:] = i[j]
        results.append((dist, ind))

    for i in range(len(results) - 1):
        assert_array_almost_equal(np.concatenate(list(results[i][0])),
                                  np.concatenate(list(results[i + 1][0]))),
        assert_array_almost_equal(np.concatenate(list(results[i][1])),
                                  np.concatenate(list(results[i + 1][1])))


def test_kneighbors_classifier(n_samples=40,
                               n_features=5,
                               n_test_pts=10,
                               n_neighbors=5,
                               random_state=0):
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    weight_func = lambda d: d ** -2

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors,
                                                 weights=weights,
                                                 algorithm=algorithm)
            knn.fit(X, y)
            epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = knn.predict(X[:n_test_pts] + epsilon)
            assert_array_equal(y_pred , y[:n_test_pts])


def test_radius_neighbors_classifier(n_samples=40,
                                     n_features=5,
                                     n_test_pts=10,
                                     radius=0.5,
                                     random_state=0):
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    weight_func = lambda d: d ** -2

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            neigh = neighbors.RadiusNeighborsClassifier(radius=radius,
                                                        weights=weights,
                                                        algorithm=algorithm)
            neigh.fit(X, y)
            epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = neigh.predict(X[:n_test_pts] + epsilon)
            assert_array_equal(y_pred, y[:n_test_pts])


def test_kneighbors_classifier_sparse(n_samples=40,
                                      n_features=5,
                                      n_test_pts=10,
                                      n_neighbors=5,
                                      random_state=0):
    """Test k-NN classifier on sparse matrices"""
    # Like the above, but with various types of sparse matrices
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    SPARSE_TYPES = (bsr_matrix, coo_matrix, csc_matrix, csr_matrix,
                    dok_matrix, lil_matrix)
    for sparsemat in SPARSE_TYPES:
        knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors,
                                             algorithm='auto')
        knn.fit(sparsemat(X), y)
        epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
        for sparsev in SPARSE_TYPES + (np.asarray,):
            X_eps = sparsev(X[:n_test_pts] + epsilon)
            y_pred = knn.predict(X_eps)
            assert_array_equal(y_pred, y[:n_test_pts])


def test_kneighbors_regressor(n_samples=40,
                              n_features=5,
                              n_test_pts=10,
                              n_neighbors=3,
                              random_state=0):
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = np.sqrt((X ** 2).sum(1))
    y /= y.max()

    y_target = y[:n_test_pts]

    weight_func = lambda d: d ** -2

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            knn = neighbors.KNeighborsRegressor(n_neighbors=n_neighbors,
                                                weights=weights,
                                                algorithm=algorithm)
            knn.fit(X, y)
            epsilon = 1E-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = knn.predict(X[:n_test_pts] + epsilon)
            assert np.all(abs(y_pred - y_target) < 0.3)


def test_radius_neighbors_regressor(n_samples=40,
                                    n_features=3,
                                    n_test_pts=10,
                                    radius=0.5,
                                    random_state=0):
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = np.sqrt((X ** 2).sum(1))
    y /= y.max()

    y_target = y[:n_test_pts]

    weight_func = lambda d: d ** -2

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            neigh = neighbors.RadiusNeighborsRegressor(radius=radius,
                                                       algorithm=algorithm)
            neigh.fit(X, y)
            epsilon = 1E-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = neigh.predict(X[:n_test_pts] + epsilon)
            assert np.all(abs(y_pred - y_target) < radius / 2)


def test_kneighbors_regressor_sparse(n_samples=40,
                                     n_features=5,
                                     n_test_pts=10,
                                     n_neighbors=5,
                                     random_state=0):
    """Test radius-neighbors classifier on sparse matrices"""
    # Like the above, but with various types of sparse matrices
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    SPARSE_TYPES = (bsr_matrix, coo_matrix, csc_matrix, csr_matrix,
                    dok_matrix, lil_matrix)
    for sparsemat in SPARSE_TYPES:
        knn = neighbors.KNeighborsRegressor(n_neighbors=n_neighbors,
                                            algorithm='auto')
        knn.fit(sparsemat(X), y)
        epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
        for sparsev in SPARSE_OR_DENSE:
            X2 = sparsev(X)
            assert (np.mean(knn.predict(X2).round() == y)
                    > 0.95)


def test_neighbors_iris():
    """
    Sanity checks on the iris dataset

    Puts three points of each label in the plane and performs a
    nearest neighbor query on points near the decision boundary.
    """

    for algorithm in ALGORITHMS:
        clf = neighbors.KNeighborsClassifier(n_neighbors=1,
                                             algorithm=algorithm)
        clf.fit(iris.data, iris.target)
        assert_array_equal(clf.predict(iris.data), iris.target)

        clf.set_params(n_neighbors=9, algorithm=algorithm)
        clf.fit(iris.data, iris.target)
        assert np.mean(clf.predict(iris.data) == iris.target) > 0.95

        rgs = neighbors.KNeighborsRegressor(n_neighbors=5,
                                            algorithm=algorithm)
        rgs.fit(iris.data, iris.target)
        assert np.mean(
            rgs.predict(iris.data).round() == iris.target) > 0.95


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

if __name__ == '__main__':
    import nose
    nose.runmodule()
