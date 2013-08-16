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

    y = np.array([-1, 1, -1, 1])
    sigma = 6

    vbtw = ecoc_utils.calculate_vbtw(X, y, sigma)
    assert_almost_equal(vbtw, 0.0736, 4)


def test_calculate_vall():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([-1, 1, -1, 1])
    sigma = 6

    vall = ecoc_utils.calculate_vall(X, y, sigma)
    assert_almost_equal(vall, 0.0736, 4)


def test_calculate_vin():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([-1, 1, -1, 1])
    sigma = 6

    vin = ecoc_utils.calculate_vin(X, y, sigma)
    assert_almost_equal(vin, 0.3228, 4)


def test_calculate_quadratic_mutal_information():
    X = np.array([[1, 1, 1, 1],
                  [2, 2, 2, 2],
                  [3, 3, 3, 3],
                  [4, 4, 4, 4]])

    y = np.array([-1, 1, -1, 1])
    sigma = 6

    qmi = ecoc_utils.calculate_quadratic_mutal_information(X, y, sigma)
    qmi_expected = 0.3228 + 0.0736 - 2*0.0736

    assert_almost_equal(qmi, qmi_expected, 4)


def test_find_quadratic_mutal_information():
    X_left = np.array([[2, 2, 2, 2],
                       [4, 4, 4, 4]])

    X_right = np.array([[1, 1, 1, 1],
                        [3, 3, 3, 3]])

    sigma = 6

    qmi = ecoc_utils.find_quadratic_mutual_information(X_left, X_right, sigma)
    qmi_expected = 9.7564

    assert_almost_equal(qmi, qmi_expected, 4)


def test_random_split():
    X = np.array([[1, 1, 1, 1],
                  [3, 3, 3, 3],
                  [3, 3, 3, 3]])
    y = np.array([1, 3, 3])

    X_left, y_left, X_right, y_right = ecoc_utils.random_split(X, y, rng)

    print X_left
    print y_left
    print X_right

    assert_array_equal(X_left, np.array([[1, 1, 1, 1]]))
    assert_array_equal(y_left, np.array([1]))

    assert_array_equal(X_right, np.array([[3, 3, 3, 3],
                                          [3, 3, 3, 3]]))
    assert_array_equal(y_right, np.array([3, 3]))


def test_add_class_to_binary_partition():
    X_left = np.array([[1, 1, 1, 1],
                       [3, 3, 3, 3],
                       [3, 3, 3, 3]])
    y_left = np.array([1, 3, 3])

    X_right = np.array([[2, 2, 2, 2],
                        [4, 4, 4, 4],
                        [4, 4, 4, 4]])
    y_right = np.array([2, 4, 4])

    X = np.concatenate((X_left, X_right), axis=0)

    sigma = ecoc_utils.calculate_sigma(X)
    current_qmi = ecoc_utils.find_quadratic_mutual_information(
        X_left, X_right, sigma)

    print(current_qmi)

    assert_almost_equal(sigma, 18.0, 1)

    # Binary-partition parameters
    bp_params = collections.namedtuple('BinaryPartitionParams',
                                       ['X_left', 'y_left',
                                        'X_right', 'y_right',
                                        'qmi', 'sigma'], verbose=False)
    bp_params.X_left = X_left
    bp_params.y_left = y_left
    bp_params.X_right = X_right
    bp_params.y_right = y_right
    bp_params.sigma = sigma
    bp_params.qmi = current_qmi

    ecoc_utils.add_class_to_binary_partition(bp_params, rng)

    assert_array_equal(bp_params.X_left, np.array([[1, 1, 1, 1],
                                                   [3, 3, 3, 3],
                                                   [3, 3, 3, 3],
                                                   [4, 4, 4, 4],
                                                   [4, 4, 4, 4]]))
    assert_array_equal(bp_params.y_left, np.array([1, 3, 3, 4, 4]))
    assert_array_equal(bp_params.X_right, np.array([[2, 2, 2, 2]]))
    assert_array_equal(bp_params.y_right, np.array([2]))


def test_remove_class_to_binary_partition():
    X_left = np.array([[1, 1, 1, 1],
                       [3, 3, 3, 3],
                       [3, 3, 3, 3]])
    y_left = np.array([1, 3, 3])

    X_right = np.array([[2, 2, 2, 2],
                        [4, 4, 4, 4],
                        [4, 4, 4, 4]])
    y_right = np.array([2, 4, 4])

    X = np.concatenate((X_left, X_right), axis=0)

    sigma = ecoc_utils.calculate_sigma(X)
    current_qmi = ecoc_utils.find_quadratic_mutual_information(
        X_left, X_right, sigma)

    print(current_qmi)

    assert_almost_equal(sigma, 18.0, 1)

    # Binary-partition parameters
    bp_params = collections.namedtuple('BinaryPartitionParams',
                                       ['X_left', 'y_left',
                                        'X_right', 'y_right',
                                        'qmi', 'sigma'], verbose=False)
    bp_params.X_left = X_left
    bp_params.y_left = y_left
    bp_params.X_right = X_right
    bp_params.y_right = y_right
    bp_params.sigma = sigma
    bp_params.qmi = current_qmi

    ecoc_utils.remove_class_to_binary_partition(bp_params, rng)

    assert_array_equal(bp_params.X_left, np.array([[3, 3, 3, 3],
                                                   [3, 3, 3, 3]]))
    assert_array_equal(bp_params.y_left, np.array([3, 3]))
    assert_array_equal(bp_params.X_right, np.array([[2, 2, 2, 2],
                                                    [4, 4, 4, 4],
                                                    [4, 4, 4, 4],
                                                    [1, 1, 1, 1]]))
    assert_array_equal(bp_params.y_right, np.array([2, 4, 4, 1]))


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

    y = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4])

    bp_params = ecoc_utils.sffs(X, y, rng)

    assert_array_equal(bp_params.X_left, np.array([[1, 1, 1, 1]]))
    assert_array_equal(bp_params.y_left, np.array([1]))
    assert_array_equal(bp_params.X_right, np.array([[2, 2, 2, 2],
                                                    [2, 2, 2, 2],
                                                    [3, 3, 3, 3],
                                                    [3, 3, 3, 3],
                                                    [3, 3, 3, 3],
                                                    [4, 4, 4, 4],
                                                    [4, 4, 4, 4],
                                                    [4, 4, 4, 4],
                                                    [4, 4, 4, 4]]))
    assert_array_equal(bp_params.y_right,
                       np.array([2, 2, 3, 3, 3, 4, 4, 4, 4]))


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

    ecoc_utils.parse_decoc(X, y, code_book, class_from, class_to, rng)

    assert_array_equal(code_book, np.array([[1, 0, 0],
                                            [-1, 1, 0],
                                            [-1, -1, 1],
                                            [-1, -1, -1]]))


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
    assert_array_equal(code_book, np.array([[1, 0, 0],
                                            [-1, 1, 0],
                                            [-1, -1, 1],
                                            [-1, -1, -1]]))


def test_create_random_codebook():
    n_classes = 4
    code_book = ecoc_utils.create_random_codebook(n_classes, n_classes, rng)
    assert_array_almost_equal(code_book,
                              np.array([[0., 1., 0., 1.],
                                        [0., 1., 0., 0.],
                                        [1., 1., 1., 1.],
                                        [0., 0., 1., 0.]]), 2)
