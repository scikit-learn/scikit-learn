from math import sqrt

import numpy as np
from scipy.sparse import csr_matrix

from sklearn import neighbors
import pytest
from numpy.testing import assert_array_equal


from sklearn.utils import check_random_state
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils.estimator_checks import check_outlier_corruption
from sklearn.utils.estimator_checks import parametrize_with_checks

from sklearn.datasets import load_iris


# load the iris dataset
# and randomly permute it
rng = check_random_state(0)
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_knnod():
    # Toy sample
    # Normal value are between [-2,-2] and [2,2]
    # The last two samples are outliers
    X = [
        [0, 0],
        [0, 2],
        [-2, -1],
        [-1, -1],
        [-1, -2],
        [1, 1],
        [1, 2],
        [2, 1],
        [4, 4],
        [9, -1],
    ]
    Y = [1, 1, 1, 1, 1, 1, 1, 1, -1, -1]

    # Test
    clf = neighbors.KNeighborsOutlierDetection(n_neighbors=3, contamination=0.2)
    clf.fit(X)
    y = clf.predict(X)
    assert_array_equal(Y, y)

    # Faster code test
    clf = neighbors.KNeighborsOutlierDetection(n_neighbors=3, contamination=0.2)
    y = clf.fit_predict(X)
    assert_array_equal(Y, y)


def test_knnod_large_k():
    X = [[0]]

    # Should not crash when k==len(X)
    clf = neighbors.KNeighborsOutlierDetection(n_neighbors=1)
    y = clf.fit_predict(X)
    assert y is not None

    # Should not crash when k>len(X)
    clf = neighbors.KNeighborsOutlierDetection(n_neighbors=10)
    y = clf.fit_predict(X)
    assert y is not None


def test_knn_distances():
    # Check the knn distances
    X = [[0.1], [0.2], [0.3], [999.9], [0.4]]
    clf = neighbors.KNeighborsOutlierDetection(n_neighbors=2)
    clf.fit(X)
    d = clf.cached_knn_distances

    # Check the position=3 is the largest
    assert max(d) == d[3]


def test_hasattr_prediction():
    # check availability of prediction methods
    clf = neighbors.KNeighborsOutlierDetection()
    assert hasattr(clf, "fit")
    assert hasattr(clf, "predict")
    assert hasattr(clf, "fit_predict")
    assert hasattr(clf, "decision_function")
    assert hasattr(clf, "knn_distance")
    assert hasattr(clf, "cached_knn_distances")


def test_csr():
    # KNeighborsOutlierDetection must support CSR inputs
    # compare results on dense and sparse data (CSR format)
    dense_X = iris.data
    csr_X = csr_matrix(dense_X)

    knnod = neighbors.KNeighborsOutlierDetection(n_neighbors=3, contamination=0.1)
    knnod.fit_predict(csr_X)
    knnod.fit(csr_X)
    csr_y = knnod.predict(csr_X)
    knnod.knn_distance(csr_X)
    knnod.decision_function(csr_X)

    knnod = neighbors.KNeighborsOutlierDetection(n_neighbors=3, contamination=0.1)
    dense_y = knnod.fit_predict(dense_X)

    assert_array_equal(dense_y, csr_y)


def test_list():
    # KNeighborsOutlierDetection must support list inputs
    # compare results on array and list
    dense_X = iris.data
    list_X = list(dense_X)

    knnod = neighbors.KNeighborsOutlierDetection(n_neighbors=3, contamination=0.1)
    knnod.fit_predict(list_X)
    knnod.fit(list_X)
    list_y = knnod.predict(list_X)
    knnod.knn_distance(list_X)
    knnod.decision_function(list_X)

    knnod = neighbors.KNeighborsOutlierDetection(n_neighbors=3, contamination=0.1)
    dense_y = knnod.fit_predict(dense_X)

    assert_array_equal(dense_y, list_y)


def test_predicted_outliers_number():
    # the number of predicted outliers should be equal to the number of
    # expected outliers. However, several ties may exist, so we accept a margin of error (+/- delta)
    X = iris.data
    n_samples = X.shape[0]
    contamination = 0.1

    clf = neighbors.KNeighborsOutlierDetection(contamination=contamination)
    y_pred = clf.fit_predict(X)

    delta = 0.01
    num_outliers = np.average(y_pred != 1)
    assert contamination - delta < num_outliers < contamination + delta
