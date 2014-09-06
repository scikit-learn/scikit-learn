from itertools import product

import numpy as np
from scipy.sparse import (bsr_matrix, coo_matrix, csc_matrix, csr_matrix,
                          dok_matrix, lil_matrix)

from sklearn.cross_validation import train_test_split
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.validation import check_random_state
from sklearn import neighbors, datasets

rng = np.random.RandomState(0)
# load and shuffle iris dataset
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# load and shuffle digits
digits = datasets.load_digits()
perm = rng.permutation(digits.target.size)
digits.data = digits.data[perm]
digits.target = digits.target[perm]

SPARSE_TYPES = (bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dok_matrix,
                lil_matrix)
SPARSE_OR_DENSE = SPARSE_TYPES + (np.asarray,)

ALGORITHMS = ('ball_tree', 'brute', 'kd_tree', 'auto')
P = (1, 2, 3, 4, np.inf)


def _weight_func(dist):
    """ Weight function to replace lambda d: d ** -2.
    The lambda function is not valid because:
    if d==0 then 0^-2 is not valid. """

    # Dist could be multidimensional, flatten it so all values
    # can be looped
    with np.errstate(divide='ignore'):
        retval = 1. / dist
    return retval ** 2


def test_unsupervised_kneighbors(n_samples=20, n_features=5,
                                 n_query_pts=2, n_neighbors=5):
    """Test unsupervised neighbors methods"""
    X = rng.rand(n_samples, n_features)

    test = rng.rand(n_query_pts, n_features)

    for p in P:
        results_nodist = []
        results = []

        for algorithm in ALGORITHMS:
            neigh = neighbors.NearestNeighbors(n_neighbors=n_neighbors,
                                               algorithm=algorithm,
                                               p=p)
            neigh.fit(X)

            results_nodist.append(neigh.kneighbors(test,
                                                   return_distance=False))
            results.append(neigh.kneighbors(test, return_distance=True))

        for i in range(len(results) - 1):
            assert_array_almost_equal(results_nodist[i], results[i][1])
            assert_array_almost_equal(results[i][0], results[i + 1][0])
            assert_array_almost_equal(results[i][1], results[i + 1][1])


def test_unsupervised_inputs():
    """test the types of valid input into NearestNeighbors"""
    X = rng.random_sample((10, 3))

    nbrs_fid = neighbors.NearestNeighbors(n_neighbors=1)
    nbrs_fid.fit(X)

    dist1, ind1 = nbrs_fid.kneighbors(X)

    nbrs = neighbors.NearestNeighbors(n_neighbors=1)

    for input in (nbrs_fid, neighbors.BallTree(X), neighbors.KDTree(X)):
        nbrs.fit(input)
        dist2, ind2 = nbrs.kneighbors(X)

        assert_array_almost_equal(dist1, dist2)
        assert_array_almost_equal(ind1, ind2)


def test_unsupervised_radius_neighbors(n_samples=20, n_features=5,
                                       n_query_pts=2, radius=0.5,
                                       random_state=0):
    """Test unsupervised radius-based query"""
    rng = np.random.RandomState(random_state)

    X = rng.rand(n_samples, n_features)

    test = rng.rand(n_query_pts, n_features)

    for p in P:
        results = []

        for algorithm in ALGORITHMS:
            neigh = neighbors.NearestNeighbors(radius=radius,
                                               algorithm=algorithm,
                                               p=p)
            neigh.fit(X)

            ind1 = neigh.radius_neighbors(test, return_distance=False)

            # sort the results: this is not done automatically for
            # radius searches
            dist, ind = neigh.radius_neighbors(test, return_distance=True)
            for (d, i, i1) in zip(dist, ind, ind1):
                j = d.argsort()
                d[:] = d[j]
                i[:] = i[j]
                i1[:] = i1[j]
            results.append((dist, ind))

            assert_array_almost_equal(np.concatenate(list(ind)),
                                      np.concatenate(list(ind1)))

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
    """Test k-neighbors classification"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .5).astype(np.int)
    y_str = y.astype(str)

    weight_func = _weight_func

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors,
                                                 weights=weights,
                                                 algorithm=algorithm)
            knn.fit(X, y)
            epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = knn.predict(X[:n_test_pts] + epsilon)
            assert_array_equal(y_pred, y[:n_test_pts])
            # Test prediction with y_str
            knn.fit(X, y_str)
            y_pred = knn.predict(X[:n_test_pts] + epsilon)
            assert_array_equal(y_pred, y_str[:n_test_pts])


def test_kneighbors_classifier_float_labels(n_samples=40, n_features=5,
                                            n_test_pts=10, n_neighbors=5,
                                            random_state=0):
    """Test k-neighbors classification"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .5).astype(np.int)

    knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors)
    knn.fit(X, y.astype(np.float))
    epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
    y_pred = knn.predict(X[:n_test_pts] + epsilon)
    assert_array_equal(y_pred, y[:n_test_pts])


def test_kneighbors_classifier_predict_proba():
    """Test KNeighborsClassifier.predict_proba() method"""
    X = np.array([[0, 2, 0],
                  [0, 2, 1],
                  [2, 0, 0],
                  [2, 2, 0],
                  [0, 0, 2],
                  [0, 0, 1]])
    y = np.array([4, 4, 5, 5, 1, 1])
    cls = neighbors.KNeighborsClassifier(n_neighbors=3, p=1)  # cityblock dist
    cls.fit(X, y)
    y_prob = cls.predict_proba(X)
    real_prob = np.array([[0, 2. / 3, 1. / 3],
                          [1. / 3, 2. / 3, 0],
                          [1. / 3, 0, 2. / 3],
                          [0, 1. / 3, 2. / 3],
                          [2. / 3, 1. / 3, 0],
                          [2. / 3, 1. / 3, 0]])
    assert_array_equal(real_prob, y_prob)
    # Check that it also works with non integer labels
    cls.fit(X, y.astype(str))
    y_prob = cls.predict_proba(X)
    assert_array_equal(real_prob, y_prob)
    # Check that it works with weights='distance'
    cls = neighbors.KNeighborsClassifier(
        n_neighbors=2, p=1, weights='distance')
    cls.fit(X, y)
    y_prob = cls.predict_proba(np.array([[0, 2, 0], [2, 2, 2]]))
    real_prob = np.array([[0, 1, 0], [0, 0.4, 0.6]])
    assert_array_almost_equal(real_prob, y_prob)


def test_radius_neighbors_classifier(n_samples=40,
                                     n_features=5,
                                     n_test_pts=10,
                                     radius=0.5,
                                     random_state=0):
    """Test radius-based classification"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .5).astype(np.int)
    y_str = y.astype(str)

    weight_func = _weight_func

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            neigh = neighbors.RadiusNeighborsClassifier(radius=radius,
                                                        weights=weights,
                                                        algorithm=algorithm)
            neigh.fit(X, y)
            epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = neigh.predict(X[:n_test_pts] + epsilon)
            assert_array_equal(y_pred, y[:n_test_pts])
            neigh.fit(X, y_str)
            y_pred = neigh.predict(X[:n_test_pts] + epsilon)
            assert_array_equal(y_pred, y_str[:n_test_pts])


def test_radius_neighbors_classifier_when_no_neighbors():
    """ Test radius-based classifier when no neighbors found.
    In this case it should rise an informative exception """

    X = np.array([[1.0, 1.0], [2.0, 2.0]])
    y = np.array([1, 2])
    radius = 0.1

    z1 = np.array([[1.01, 1.01], [2.01, 2.01]])  # no outliers
    z2 = np.array([[1.01, 1.01], [1.4, 1.4]])    # one outlier

    weight_func = _weight_func

    for outlier_label in [0, -1, None]:
        for algorithm in ALGORITHMS:
            for weights in ['uniform', 'distance', weight_func]:
                rnc = neighbors.RadiusNeighborsClassifier
                clf = rnc(radius=radius, weights=weights, algorithm=algorithm,
                          outlier_label=outlier_label)
                clf.fit(X, y)
                assert_array_equal(np.array([1, 2]),
                                   clf.predict(z1))
                if outlier_label is None:
                    assert_raises(ValueError, clf.predict, z2)
                elif False:
                    assert_array_equal(np.array([1, outlier_label]),
                                       clf.predict(z2))


def test_radius_neighbors_classifier_outlier_labeling():
    """ Test radius-based classifier when no neighbors found and outliers
    are labeled. """

    X = np.array([[1.0, 1.0], [2.0, 2.0]])
    y = np.array([1, 2])
    radius = 0.1

    z1 = np.array([[1.01, 1.01], [2.01, 2.01]])  # no outliers
    z2 = np.array([[1.01, 1.01], [1.4, 1.4]])    # one outlier
    correct_labels1 = np.array([1, 2])
    correct_labels2 = np.array([1, -1])

    weight_func = _weight_func

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            clf = neighbors.RadiusNeighborsClassifier(radius=radius,
                                                      weights=weights,
                                                      algorithm=algorithm,
                                                      outlier_label=-1)
            clf.fit(X, y)
            assert_array_equal(correct_labels1, clf.predict(z1))
            assert_array_equal(correct_labels2, clf.predict(z2))


def test_radius_neighbors_classifier_zero_distance():
    """ Test radius-based classifier, when distance to a sample is zero. """

    X = np.array([[1.0, 1.0], [2.0, 2.0]])
    y = np.array([1, 2])
    radius = 0.1

    z1 = np.array([[1.01, 1.01], [2.0, 2.0]])
    correct_labels1 = np.array([1, 2])

    weight_func = _weight_func

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            clf = neighbors.RadiusNeighborsClassifier(radius=radius,
                                                      weights=weights,
                                                      algorithm=algorithm)
            clf.fit(X, y)
            assert_array_equal(correct_labels1, clf.predict(z1))


def test_RadiusNeighborsClassifier_multioutput():
    """Test k-NN classifier on multioutput data"""
    rng = check_random_state(0)
    n_features = 2
    n_samples = 40
    n_output = 3

    X = rng.rand(n_samples, n_features)
    y = rng.randint(0, 3, (n_samples, n_output))

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    weights = [None, 'uniform', 'distance', _weight_func]

    for algorithm, weights in product(ALGORITHMS, weights):
        # Stack single output prediction
        y_pred_so = []
        for o in range(n_output):
            rnn = neighbors.RadiusNeighborsClassifier(weights=weights,
                                                      algorithm=algorithm)
            rnn.fit(X_train, y_train[:, o])
            y_pred_so.append(rnn.predict(X_test))

        y_pred_so = np.vstack(y_pred_so).T
        assert_equal(y_pred_so.shape, y_test.shape)

        # Multioutput prediction
        rnn_mo = neighbors.RadiusNeighborsClassifier(weights=weights,
                                                     algorithm=algorithm)
        rnn_mo.fit(X_train, y_train)
        y_pred_mo = rnn_mo.predict(X_test)

        assert_equal(y_pred_mo.shape, y_test.shape)
        assert_array_almost_equal(y_pred_mo, y_pred_so)


def test_kneighbors_classifier_sparse(n_samples=40,
                                      n_features=5,
                                      n_test_pts=10,
                                      n_neighbors=5,
                                      random_state=0):
    """Test k-NN classifier on sparse matrices"""
    # Like the above, but with various types of sparse matrices
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    X *= X > .2
    y = ((X ** 2).sum(axis=1) < .5).astype(np.int)

    for sparsemat in SPARSE_TYPES:
        knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors,
                                             algorithm='auto')
        knn.fit(sparsemat(X), y)
        epsilon = 1e-5 * (2 * rng.rand(1, n_features) - 1)
        for sparsev in SPARSE_TYPES + (np.asarray,):
            X_eps = sparsev(X[:n_test_pts] + epsilon)
            y_pred = knn.predict(X_eps)
            assert_array_equal(y_pred, y[:n_test_pts])


def test_KNeighborsClassifier_multioutput():
    """Test k-NN classifier on multioutput data"""
    rng = check_random_state(0)
    n_features = 5
    n_samples = 50
    n_output = 3

    X = rng.rand(n_samples, n_features)
    y = rng.randint(0, 3, (n_samples, n_output))

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    weights = [None, 'uniform', 'distance', _weight_func]

    for algorithm, weights in product(ALGORITHMS, weights):
        # Stack single output prediction
        y_pred_so = []
        y_pred_proba_so = []
        for o in range(n_output):
            knn = neighbors.KNeighborsClassifier(weights=weights,
                                                 algorithm=algorithm)
            knn.fit(X_train, y_train[:, o])
            y_pred_so.append(knn.predict(X_test))
            y_pred_proba_so.append(knn.predict_proba(X_test))

        y_pred_so = np.vstack(y_pred_so).T
        assert_equal(y_pred_so.shape, y_test.shape)
        assert_equal(len(y_pred_proba_so), n_output)

        # Multioutput prediction
        knn_mo = neighbors.KNeighborsClassifier(weights=weights,
                                                algorithm=algorithm)
        knn_mo.fit(X_train, y_train)
        y_pred_mo = knn_mo.predict(X_test)

        assert_equal(y_pred_mo.shape, y_test.shape)
        assert_array_almost_equal(y_pred_mo, y_pred_so)

        # Check proba
        y_pred_proba_mo = knn_mo.predict_proba(X_test)
        assert_equal(len(y_pred_proba_mo), n_output)

        for proba_mo, proba_so in zip(y_pred_proba_mo, y_pred_proba_so):
            assert_array_almost_equal(proba_mo, proba_so)


def test_kneighbors_regressor(n_samples=40,
                              n_features=5,
                              n_test_pts=10,
                              n_neighbors=3,
                              random_state=0):
    """Test k-neighbors regression"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = np.sqrt((X ** 2).sum(1))
    y /= y.max()

    y_target = y[:n_test_pts]

    weight_func = _weight_func

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            knn = neighbors.KNeighborsRegressor(n_neighbors=n_neighbors,
                                                weights=weights,
                                                algorithm=algorithm)
            knn.fit(X, y)
            epsilon = 1E-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = knn.predict(X[:n_test_pts] + epsilon)
            assert_true(np.all(abs(y_pred - y_target) < 0.3))


def test_KNeighborsRegressor_multioutput_uniform_weight():
    """Test k-neighbors in multi-output regression with uniform weight"""
    rng = check_random_state(0)
    n_features = 5
    n_samples = 40
    n_output = 4

    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples, n_output)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    for algorithm, weights in product(ALGORITHMS, [None, 'uniform']):
        knn = neighbors.KNeighborsRegressor(weights=weights,
                                            algorithm=algorithm)
        knn.fit(X_train, y_train)

        neigh_idx = knn.kneighbors(X_test, return_distance=False)
        y_pred_idx = np.array([np.mean(y_train[idx], axis=0)
                               for idx in neigh_idx])

        y_pred = knn.predict(X_test)

        assert_equal(y_pred.shape, y_test.shape)
        assert_equal(y_pred_idx.shape, y_test.shape)
        assert_array_almost_equal(y_pred, y_pred_idx)


def test_kneighbors_regressor_multioutput(n_samples=40,
                                          n_features=5,
                                          n_test_pts=10,
                                          n_neighbors=3,
                                          random_state=0):
    """Test k-neighbors in multi-output regression"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = np.sqrt((X ** 2).sum(1))
    y /= y.max()
    y = np.vstack([y, y]).T

    y_target = y[:n_test_pts]

    weights = ['uniform', 'distance', _weight_func]
    for algorithm, weights in product(ALGORITHMS, weights):
        knn = neighbors.KNeighborsRegressor(n_neighbors=n_neighbors,
                                            weights=weights,
                                            algorithm=algorithm)
        knn.fit(X, y)
        epsilon = 1E-5 * (2 * rng.rand(1, n_features) - 1)
        y_pred = knn.predict(X[:n_test_pts] + epsilon)
        assert_equal(y_pred.shape, y_target.shape)

        assert_true(np.all(np.abs(y_pred - y_target) < 0.3))


def test_radius_neighbors_regressor(n_samples=40,
                                    n_features=3,
                                    n_test_pts=10,
                                    radius=0.5,
                                    random_state=0):
    """Test radius-based neighbors regression"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = np.sqrt((X ** 2).sum(1))
    y /= y.max()

    y_target = y[:n_test_pts]

    weight_func = _weight_func

    for algorithm in ALGORITHMS:
        for weights in ['uniform', 'distance', weight_func]:
            neigh = neighbors.RadiusNeighborsRegressor(radius=radius,
                                                       weights=weights,
                                                       algorithm=algorithm)
            neigh.fit(X, y)
            epsilon = 1E-5 * (2 * rng.rand(1, n_features) - 1)
            y_pred = neigh.predict(X[:n_test_pts] + epsilon)
            assert_true(np.all(abs(y_pred - y_target) < radius / 2))


def test_RadiusNeighborsRegressor_multioutput_with_uniform_weight():
    """Test radius neighbors in multi-output regression (uniform weight)"""

    rng = check_random_state(0)
    n_features = 5
    n_samples = 40
    n_output = 4

    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples, n_output)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    for algorithm, weights in product(ALGORITHMS, [None, 'uniform']):

        rnn = neighbors. RadiusNeighborsRegressor(weights=weights,
                                                  algorithm=algorithm)
        rnn.fit(X_train, y_train)

        neigh_idx = rnn.radius_neighbors(X_test, return_distance=False)
        y_pred_idx = np.array([np.mean(y_train[idx], axis=0)
                               for idx in neigh_idx])

        y_pred_idx = np.array(y_pred_idx)
        y_pred = rnn.predict(X_test)

        assert_equal(y_pred_idx.shape, y_test.shape)
        assert_equal(y_pred.shape, y_test.shape)
        assert_array_almost_equal(y_pred, y_pred_idx)


def test_RadiusNeighborsRegressor_multioutput(n_samples=40,
                                              n_features=5,
                                              n_test_pts=10,
                                              n_neighbors=3,
                                              random_state=0):
    """Test k-neighbors in multi-output regression with various weight"""
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = np.sqrt((X ** 2).sum(1))
    y /= y.max()
    y = np.vstack([y, y]).T

    y_target = y[:n_test_pts]
    weights = ['uniform', 'distance', _weight_func]

    for algorithm, weights in product(ALGORITHMS, weights):
        rnn = neighbors.RadiusNeighborsRegressor(n_neighbors=n_neighbors,
                                                 weights=weights,
                                                 algorithm=algorithm)
        rnn.fit(X, y)
        epsilon = 1E-5 * (2 * rng.rand(1, n_features) - 1)
        y_pred = rnn.predict(X[:n_test_pts] + epsilon)

        assert_equal(y_pred.shape, y_target.shape)
        assert_true(np.all(np.abs(y_pred - y_target) < 0.3))


def test_kneighbors_regressor_sparse(n_samples=40,
                                     n_features=5,
                                     n_test_pts=10,
                                     n_neighbors=5,
                                     random_state=0):
    """Test radius-based regression on sparse matrices"""
    # Like the above, but with various types of sparse matrices
    rng = np.random.RandomState(random_state)
    X = 2 * rng.rand(n_samples, n_features) - 1
    y = ((X ** 2).sum(axis=1) < .25).astype(np.int)

    for sparsemat in SPARSE_TYPES:
        knn = neighbors.KNeighborsRegressor(n_neighbors=n_neighbors,
                                            algorithm='auto')
        knn.fit(sparsemat(X), y)
        for sparsev in SPARSE_OR_DENSE:
            X2 = sparsev(X)
            assert_true(np.mean(knn.predict(X2).round() == y) > 0.95)


def test_neighbors_iris():
    """Sanity checks on the iris dataset

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
        assert_true(np.mean(clf.predict(iris.data) == iris.target) > 0.95)

        rgs = neighbors.KNeighborsRegressor(n_neighbors=5, algorithm=algorithm)
        rgs.fit(iris.data, iris.target)
        assert_true(np.mean(rgs.predict(iris.data).round() == iris.target)
                    > 0.95)


def test_neighbors_digits():
    """Sanity check on the digits dataset

    the 'brute' algorithm has been observed to fail if the input
    dtype is uint8 due to overflow in distance calculations.
    """

    X = digits.data.astype('uint8')
    Y = digits.target
    (n_samples, n_features) = X.shape
    train_test_boundary = int(n_samples * 0.8)
    train = np.arange(0, train_test_boundary)
    test = np.arange(train_test_boundary, n_samples)
    (X_train, Y_train, X_test, Y_test) = X[train], Y[train], X[test], Y[test]

    clf = neighbors.KNeighborsClassifier(n_neighbors=1, algorithm='brute')
    score_uint8 = clf.fit(X_train, Y_train).score(X_test, Y_test)
    score_float = clf.fit(X_train.astype(float), Y_train).score(
        X_test.astype(float), Y_test)
    assert_equal(score_uint8, score_float)


def test_kneighbors_graph():
    """Test kneighbors_graph to build the k-Nearest Neighbor graph."""
    X = np.array([[0, 1], [1.01, 1.], [2, 0]])

    # n_neighbors = 1
    A = neighbors.kneighbors_graph(X, 1, mode='connectivity')
    assert_array_equal(A.toarray(), np.eye(A.shape[0]))

    A = neighbors.kneighbors_graph(X, 1, mode='distance')
    assert_array_almost_equal(
        A.toarray(),
        [[0.00, 1.01, 0.],
         [1.01, 0., 0.],
         [0.00, 1.40716026, 0.]])

    # n_neighbors = 2
    A = neighbors.kneighbors_graph(X, 2, mode='connectivity')
    assert_array_equal(
        A.toarray(),
        [[1., 1., 0.],
         [1., 1., 0.],
         [0., 1., 1.]])

    A = neighbors.kneighbors_graph(X, 2, mode='distance')
    assert_array_almost_equal(
        A.toarray(),
        [[0., 1.01, 2.23606798],
         [1.01, 0., 1.40716026],
         [2.23606798, 1.40716026, 0.]])

    # n_neighbors = 3
    A = neighbors.kneighbors_graph(X, 3, mode='connectivity')
    assert_array_almost_equal(
        A.toarray(),
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]])


def test_kneighbors_graph_sparse(seed=36):
    """Test kneighbors_graph to build the k-Nearest Neighbor graph
    for sparse input."""
    rng = np.random.RandomState(seed)
    X = rng.randn(10, 10)
    Xcsr = csr_matrix(X)

    for n_neighbors in [1, 2, 3]:
        for mode in ["connectivity", "distance"]:
            assert_array_almost_equal(
                neighbors.kneighbors_graph(X,
                                           n_neighbors,
                                           mode=mode).toarray(),
                neighbors.kneighbors_graph(Xcsr,
                                           n_neighbors,
                                           mode=mode).toarray())


def test_radius_neighbors_graph():
    """Test radius_neighbors_graph to build the Nearest Neighbor graph."""
    X = np.array([[0, 1], [1.01, 1.], [2, 0]])

    A = neighbors.radius_neighbors_graph(X, 1.5, mode='connectivity')
    assert_array_equal(
        A.toarray(),
        [[1., 1., 0.],
         [1., 1., 1.],
         [0., 1., 1.]])

    A = neighbors.radius_neighbors_graph(X, 1.5, mode='distance')
    assert_array_almost_equal(
        A.toarray(),
        [[0., 1.01, 0.],
         [1.01, 0., 1.40716026],
         [0., 1.40716026, 0.]])


def test_radius_neighbors_graph_sparse(seed=36):
    """Test radius_neighbors_graph to build the Nearest Neighbor graph
    for sparse input."""
    rng = np.random.RandomState(seed)
    X = rng.randn(10, 10)
    Xcsr = csr_matrix(X)

    for n_neighbors in [1, 2, 3]:
        for mode in ["connectivity", "distance"]:
            assert_array_almost_equal(
                neighbors.radius_neighbors_graph(X,
                                                 n_neighbors,
                                                 mode=mode).toarray(),
                neighbors.radius_neighbors_graph(Xcsr,
                                                 n_neighbors,
                                                 mode=mode).toarray())


def test_neighbors_badargs():
    """Test bad argument values: these should all raise ValueErrors"""
    assert_raises(ValueError,
                  neighbors.NearestNeighbors,
                  algorithm='blah')

    X = rng.random_sample((10, 2))
    Xsparse = csr_matrix(X)
    y = np.ones(10)

    for cls in (neighbors.KNeighborsClassifier,
                neighbors.RadiusNeighborsClassifier,
                neighbors.KNeighborsRegressor,
                neighbors.RadiusNeighborsRegressor):
        assert_raises(ValueError,
                      cls,
                      weights='blah')
        assert_raises(ValueError,
                      cls, p=-1)
        assert_raises(ValueError,
                      cls, algorithm='blah')
        nbrs = cls(algorithm='ball_tree', metric='haversine')
        assert_raises(ValueError,
                      nbrs.predict,
                      X)
        assert_raises(ValueError,
                      ignore_warnings(nbrs.fit),
                      Xsparse, y)
        nbrs = cls()
        assert_raises(ValueError,
                      nbrs.fit,
                      np.ones((0, 2)), np.ones(0))
        assert_raises(ValueError,
                      nbrs.fit,
                      X[:, :, None], y)
        nbrs.fit(X, y)
        assert_raises(ValueError,
                      nbrs.predict,
                      [])

    nbrs = neighbors.NearestNeighbors().fit(X)

    assert_raises(ValueError,
                  nbrs.kneighbors_graph,
                  X, mode='blah')
    assert_raises(ValueError,
                  nbrs.radius_neighbors_graph,
                  X, mode='blah')


def test_neighbors_metrics(n_samples=20, n_features=3,
                           n_query_pts=2, n_neighbors=5):
    """Test computing the neighbors for various metrics"""
    # create a symmetric matrix
    V = rng.rand(n_features, n_features)
    VI = np.dot(V, V.T)

    metrics = [('euclidean', {}),
               ('manhattan', {}),
               ('minkowski', dict(p=1)),
               ('minkowski', dict(p=2)),
               ('minkowski', dict(p=3)),
               ('minkowski', dict(p=np.inf)),
               ('chebyshev', {}),
               ('seuclidean', dict(V=rng.rand(n_features))),
               ('wminkowski', dict(p=3, w=rng.rand(n_features))),
               ('mahalanobis', dict(VI=VI))]
    algorithms = ['brute', 'ball_tree', 'kd_tree']
    X = rng.rand(n_samples, n_features)

    test = rng.rand(n_query_pts, n_features)

    for metric, metric_params in metrics:
        results = []
        p = metric_params.pop('p', 2)
        for algorithm in algorithms:
            # KD tree doesn't support all metrics
            if (algorithm == 'kd_tree' and
                    metric not in neighbors.KDTree.valid_metrics):
                assert_raises(ValueError,
                              neighbors.NearestNeighbors,
                              algorithm=algorithm,
                              metric=metric, metric_params=metric_params)
                continue

            neigh = neighbors.NearestNeighbors(n_neighbors=n_neighbors,
                                               algorithm=algorithm,
                                               metric=metric, p=p,
                                               metric_params=metric_params)
            neigh.fit(X)
            results.append(neigh.kneighbors(test, return_distance=True))

        assert_array_almost_equal(results[0][0], results[1][0])
        assert_array_almost_equal(results[0][1], results[1][1])


def test_callable_metric():
    metric = lambda x1, x2: np.sqrt(np.sum(x1 ** 2 + x2 ** 2))

    X = np.random.RandomState(42).rand(20, 2)
    nbrs1 = neighbors.NearestNeighbors(3, algorithm='auto', metric=metric)
    nbrs2 = neighbors.NearestNeighbors(3, algorithm='brute', metric=metric)

    nbrs1.fit(X)
    nbrs2.fit(X)

    dist1, ind1 = nbrs1.kneighbors(X)
    dist2, ind2 = nbrs2.kneighbors(X)

    assert_array_almost_equal(dist1, dist2)


def test_metric_params_interface():
    assert_warns(DeprecationWarning, neighbors.KNeighborsClassifier,
                 metric='wminkowski', w=np.ones(10))
    assert_warns(SyntaxWarning, neighbors.KNeighborsClassifier,
                 metric_params={'p': 3})


if __name__ == '__main__':
    import nose
    nose.runmodule()
