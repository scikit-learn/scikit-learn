import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises

import collections
import sklearn.ecoc_utils as ecoc_utils
from sklearn import datasets

iris = datasets.load_iris()
rng = np.random.RandomState(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_calculate_parzen_estimate():
    x = np.array([1, 1, 1, 1])
    y = np.array([1, 2, 3, 4])
    sigma = 4

    parzen_estimate = ecoc_utils.calculate_parzen_estimate(x, y, sigma)
    assert_almost_equal(parzen_estimate, 0.280, 3)


def test_calculate_parzen_estimate_different_sizes():
    x = np.array([1, 1, 1, 1])
    y = np.array([1, 2, 3])
    sigma = 4

    assert_raises(ValueError,
                  ecoc_utils.calculate_parzen_estimate,
                  x, y, sigma)


def test_calculate_sigma():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [4, 4, 4, 4]])

    sigma = ecoc_utils.calculate_sigma(X)
    assert_equal(sigma, 0.5*36)


def test_calculate_vbtw():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 2, 3])

    y_left = np.array([0, 2])
    y_right = np.array([1, 3])

    sigma = 6

    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    vbtw = ecoc_utils.calculate_vbtw(y_left, y_right,
                                     mtInfMat, classFrequencies)
    assert_almost_equal(vbtw, 0.0737, 4)


def test_calculate_vall():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 2, 3])

    y_left = np.array([0, 2])
    y_right = np.array([1, 3])

    sigma = 6

    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    vall = ecoc_utils.calculate_vall(y_left, y_right,
                                     mtInfMat, classFrequencies)
    assert_almost_equal(vall, 0.0737, 4)


def test_calculate_vin():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 2, 3])

    y_left = np.array([0, 2])
    y_right = np.array([1, 3])

    sigma = 6

    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    vin = ecoc_utils.calculate_vin(y_left, y_right,
                                   mtInfMat, classFrequencies)
    assert_almost_equal(vin, 0.0807, 4)


def test_calculate_qm_information():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 2, 3])

    y_left = np.array([0, 2])
    y_right = np.array([1, 3])

    sigma = 6

    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    qmi = ecoc_utils.calculate_qm_information(y_left, y_right,
                                              mtInfMat, classFrequencies)

    qmi_expected = 0.0807 + 0.0737 - 2*0.0737

    assert_almost_equal(qmi, qmi_expected, 4)


def test_random_split():
    y = np.array([1, 2, 3])

    y_left, y_right = ecoc_utils.random_split(y, rng)

    assert_array_equal(y_left, np.array([1]))
    assert_array_equal(y_right, np.array([3, 2]))


def test_add_class_to_binary_partition():
    X = np.array([[1, 1, 1, 1],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [2, 2, 2, 2],
                  [4, 4, 4, 4],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 2, 3])
    y_left = np.array([0, 2])
    y_right = np.array([1, 3])

    sigma = ecoc_utils.calculate_sigma(X)
    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    current_qmi = ecoc_utils.calculate_qm_information(y_left, y_right,
                                                      mtInfMat,
                                                      classFrequencies)

    assert_almost_equal(sigma, 18.0, 1)

    # Binary-partition parameters
    bp_params = collections.namedtuple('BinaryPartitionParams',
                                       ['y_left', 'y_right',
                                        'qmi', 'sigma'], verbose=False)

    bp_params.y_left = y_left
    bp_params.y_right = y_right
    bp_params.sigma = sigma
    bp_params.qmi = current_qmi

    ecoc_utils.add_class_to_binary_partition(bp_params, rng,
                                             mtInfMat, classFrequencies)

    assert_array_equal(bp_params.y_left, np.array([0, 2, 3]))
    assert_array_equal(bp_params.y_right, np.array([1]))


def test_remove_class_to_binary_partition():

    X = np.array([[1, 1, 1, 1],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [2, 2, 2, 2],
                  [4, 4, 4, 4],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 2, 3])
    y_left = np.array([0, 2])
    y_right = np.array([1, 3])

    sigma = ecoc_utils.calculate_sigma(X)
    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    current_qmi = ecoc_utils.calculate_qm_information(y_left, y_right,
                                                      mtInfMat,
                                                      classFrequencies)

    assert_almost_equal(sigma, 18.0, 1)

    # Binary-partition parameters
    bp_params = collections.namedtuple('BinaryPartitionParams',
                                       ['y_left', 'y_right',
                                        'qmi', 'sigma'], verbose=False)

    bp_params.y_left = y_left
    bp_params.y_right = y_right
    bp_params.sigma = sigma
    bp_params.qmi = current_qmi

    ecoc_utils.remove_class_to_binary_partition(bp_params, rng,
                                                mtInfMat, classFrequencies)

    assert_array_equal(bp_params.y_left, np.array([2]))
    assert_array_equal(bp_params.y_right, np.array([1, 3, 0]))


def test_sffs():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4],
                  [4, 4, 4, 4],
                  [4, 4, 4, 4],
                  [4, 4, 4, 4]])

    y = np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3])

    sigma = ecoc_utils.calculate_sigma(X)
    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)
    y_left, y_right = ecoc_utils.sffs(np.unique(y),
                                      rng, mtInfMat, classFrequencies)

    assert_array_equal(y_left, np.array([1, 0, 2]))
    assert_array_equal(y_right, np.array([3]))


def test_parse_decoc():
    X = np.array([[0, 0, 0, 0],
                  [1, 1, 1, 1],
                  [1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [2, 2, 2, 2],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3]])

    y = np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3])
    n_classes = 4

    code_book = np.ones(shape=(n_classes, n_classes-1), dtype=np.int)

    class_from = 0
    class_to = n_classes-1

    sigma = ecoc_utils.calculate_sigma(X)
    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)

    ecoc_utils.parse_decoc(np.unique(y), code_book, class_from, class_to, rng,
                           mtInfMat, classFrequencies)

    assert_array_equal(code_book, np.array([[1, -1, -1],
                                            [1, 1, 0],
                                            [1, -1, 1],
                                            [-1, 0, 0]]))


def test_create_decoc_codebook():
    X = np.array([[0, 0, 0, 0],
                  [1, 1, 1, 1],
                  [1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [2, 2, 2, 2],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3]])

    y = np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3])
    n_classes = 4

    code_book = ecoc_utils.create_decoc_codebook(n_classes, X, y, rng)
    assert_array_equal(code_book, np.array([[1, 1, -1],
                                            [1, -1, 0],
                                            [1, 1, 1],
                                            [-1, 0, 0]]))


def test_create_random_codebook():
    n_classes = 4
    code_book = ecoc_utils.create_random_codebook(n_classes, n_classes, rng)
    assert_array_almost_equal(code_book,
                              np.array([[1., 0., 1., 1.],
                                        [1., 1., 1., 1.],
                                        [1., 1., 0., 1.],
                                        [0., 0., 1., 0.]]), 2)


def test_calculate_class_frequencies():
    y = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4])
    classFrequencies = ecoc_utils.calculate_class_frequencies(y)
    expected_counts = np.array([1, 2, 3, 4])

    assert_array_equal(classFrequencies, expected_counts)


def test_calculate_mutual_information_matrix():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3]])

    y = np.array([1, 2, 3])

    sigma = 1
    mtInfMat = ecoc_utils.calculate_mutual_information_matrix(X, y, sigma)

    expected_mtInfMat = np.array([[2.420e-01, 7.500e-02, 1.916e-02],
                                  [7.498e-02, 4.432e-03, 1.994e-04],
                                  [1.916e-02, 1.994e-04, 1.487e-06]])

    assert_array_almost_equal(mtInfMat, expected_mtInfMat, 3)

def test_calculate_mutual_information_atom():
    X_left = np.array([[1, 1, 1, 1],
                       [3, 3, 3, 3],
                       [3, 3, 3, 3]])

    X_right = np.array([[2, 2, 2, 2],
                        [4, 4, 4, 4],
                        [4, 4, 4, 4]])

    sigma = 0.05

    mutInfAt = ecoc_utils.calculate_mutual_information_atom(X_left,
                                                            X_right, sigma)

    assert_almost_equal(mutInfAt, 0.0087, 4)
