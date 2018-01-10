import sys
import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false

from sklearn import datasets
from sklearn.neighbors import LargeMarginNearestNeighbor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors.lmnn import make_lmnn_pipeline
from sklearn.neighbors.lmnn import _paired_distances_blockwise
from sklearn.neighbors.lmnn import _euclidean_distances_without_checks
from sklearn.metrics.pairwise import paired_euclidean_distances
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import train_test_split
from sklearn.utils.extmath import row_norms
from sklearn.externals.six.moves import cStringIO as StringIO


rng = np.random.RandomState(0)
# load and shuffle iris dataset
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris_data = iris.data[perm]
iris_target = iris.target[perm]

# load and shuffle digits
digits = datasets.load_digits()
perm = rng.permutation(digits.target.size)
digits_data = digits.data[perm]
digits_target = digits.target[perm]


def test_neighbors_iris():
    # Sanity checks on the iris dataset
    # Puts three points of each label in the plane and performs a
    # nearest neighbor query on points near the decision boundary.

    lmnn = LargeMarginNearestNeighbor(n_neighbors=1)
    lmnn.fit(iris_data, iris_target)
    knn = KNeighborsClassifier(n_neighbors=lmnn.n_neighbors_)
    LX = lmnn.transform(iris_data)
    knn.fit(LX, iris_target)
    y_pred = knn.predict(LX)

    assert_array_equal(y_pred, iris_target)

    lmnn.set_params(n_neighbors=9)
    lmnn.fit(iris_data, iris_target)
    knn = KNeighborsClassifier(n_neighbors=lmnn.n_neighbors_)
    knn.fit(LX, iris_target)

    assert_true(knn.score(LX, iris_target) > 0.95)


def test_neighbors_digits():
    # Sanity check on the digits dataset
    # the 'brute' algorithm has been observed to fail if the input
    # dtype is uint8 due to overflow in distance calculations.

    X = digits_data.astype('uint8')
    y = digits_target
    n_samples, n_features = X.shape
    train_test_boundary = int(n_samples * 0.8)
    train = np.arange(0, train_test_boundary)
    test = np.arange(train_test_boundary, n_samples)
    X_train, y_train, X_test, y_test = X[train], y[train], X[test], y[test]

    k = 1
    lmnn = LargeMarginNearestNeighbor(n_neighbors=k, max_iter=30)
    lmnn.fit(X_train, y_train)
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(lmnn.transform(X_train), y_train)
    score_uint8 = knn.score(lmnn.transform(X_test), y_test)

    knn.fit(lmnn.transform(X_train.astype(float)), y_train)
    score_float = knn.score(lmnn.transform(X_test.astype(float)), y_test)

    assert_equal(score_uint8, score_float)


def test_params_validation():
    # Test that invalid parameters raise value error
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]
    LMNN = LargeMarginNearestNeighbor

    # TypeError
    assert_raises(TypeError, LMNN(n_neighbors=1.3).fit, X, y)
    assert_raises(TypeError, LMNN(max_iter='21').fit, X, y)
    assert_raises(TypeError, LMNN(verbose='true').fit, X, y)
    assert_raises(TypeError, LMNN(max_impostors=23.1).fit, X, y)
    assert_raises(TypeError, LMNN(tol=1).fit, X, y)
    assert_raises(TypeError, LMNN(n_components='invalid').fit, X, y)
    assert_raises(TypeError, LMNN(n_jobs='yes').fit, X, y)
    assert_raises(TypeError, LMNN(warm_start=1).fit, X, y)
    assert_raises(TypeError, LMNN(impostor_store=0.5).fit, X, y)
    assert_raises(TypeError, LMNN(neighbors_params=65).fit, X, y)
    assert_raises(TypeError, LMNN(weight_push_loss='0.3').fit, X, y)

    # ValueError
    assert_raise_message(ValueError,
                         "`init` must be 'pca', 'identity', or a numpy "
                         "array of shape (n_components, n_features).",
                         LMNN(init=1).fit, X, y)

    assert_raise_message(ValueError,
                         '`n_neighbors`= -1, must be >= 1.',
                         LMNN(n_neighbors=-1).fit, X, y)

    assert_raise_message(ValueError,
                         '`n_neighbors`= {}, must be <= {}.'
                         .format(X.shape[0], X.shape[0] - 1),
                         LMNN(n_neighbors=X.shape[0]).fit, X, y)

    assert_raise_message(ValueError,
                         '`max_iter`= -1, must be >= 1.',
                         LMNN(max_iter=-1).fit, X, y)
    assert_raise_message(ValueError,
                         '`max_impostors`= -1, must be >= 1.',
                         LMNN(max_impostors=-1).fit, X, y)
    assert_raise_message(ValueError,
                         "`impostor_store` must be 'auto', 'sparse' "
                         "or 'list'.",
                         LMNN(impostor_store='dense').fit, X, y)

    assert_raise_message(ValueError,
                         '`weight_push_loss`= 2.0, must be <= 1.0.',
                         LMNN(weight_push_loss=2.).fit, X, y)

    assert_raise_message(ValueError,
                         '`weight_push_loss` cannot be zero.',
                         LMNN(weight_push_loss=0.).fit, X, y)

    init = np.random.rand(5, 3)
    assert_raise_message(ValueError,
                         'The output dimensionality ({}) of the given linear '
                         'transformation `init` cannot be greater than its '
                         'input dimensionality ({}).'
                         .format(init.shape[0], init.shape[1]),
                         LMNN(init=init).fit, X, y)

    n_components = 10
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) cannot be greater '
                         'than the given data dimensionality ({})!'
                         .format(n_components, X.shape[1]),
                         LMNN(n_components=n_components).fit, X, y)

    n_jobs = 0
    assert_raise_message(ValueError,
                         'Parameter n_jobs == 0 has no meaning.',
                         LMNN(n_jobs=n_jobs).fit, X, y)

    # test min_class_size < 2
    y = [1, 1, 1, 2]
    assert_raise_message(ValueError,
                         'LargeMarginNearestNeighbor needs at least 2 '
                         'non-singleton classes, got 1.',
                         LMNN(n_neighbors=1).fit, X, y)


def test_same_lmnn_parallel():
    X, y = datasets.make_classification(n_samples=30, n_features=5,
                                        n_redundant=0, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y)

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3)
    lmnn.fit(X_train, y_train)
    components = lmnn.components_

    lmnn.set_params(n_jobs=3)
    lmnn.fit(X_train, y_train)
    components_parallel = lmnn.components_

    assert_array_almost_equal(components, components_parallel)


def test_init_transformation_validation():
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]

    # Fail if transformation input dimension does not match inputs dimensions
    init = np.array([[1, 2], [3, 4]])
    assert_raise_message(ValueError,
                         'The input dimensionality ({}) of the given '
                         'linear transformation `init` must match the '
                         'dimensionality of the given inputs `X` ({}).'
                         .format(init.shape[1], X.shape[1]),
                         LargeMarginNearestNeighbor(init=init,
                                                    n_neighbors=1).fit,
                         X, y)

    # Fail if transformation output dimension is larger than
    # transformation input dimension
    init = np.random.rand(4, 3)
    # len(init) > len(init[0])
    assert_raise_message(ValueError,
                         'The output dimensionality ({}) of the given '
                         'linear transformation `init` cannot be '
                         'greater than its input dimensionality ({}).'
                         .format(init.shape[0], init.shape[1]),
                         LargeMarginNearestNeighbor(init=init,
                                                    n_neighbors=1).fit,
                         X, y)

    # Pass otherwise
    init = np.arange(9).reshape(3, 3)
    LargeMarginNearestNeighbor(init=init, n_neighbors=1).fit(X, y)


def test_n_neighbors():
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]

    lmnn = LargeMarginNearestNeighbor(n_neighbors=2)
    assert_warns_message(UserWarning,
                         '`n_neighbors` (=2) is not less than the number of '
                         'samples in the smallest non-singleton class (=2). '
                         '`n_neighbors_` will be set to 1 for estimation.',
                         lmnn.fit, X, y)


def test_n_components():
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]

    init = np.random.rand(X.shape[1] - 1, 3)

    # n_components = X.shape[1] != transformation.shape[0]
    n_components = X.shape[1]
    lmnn = LargeMarginNearestNeighbor(init=init, n_components=n_components)
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) does not match '
                         'the output dimensionality of the given '
                         'linear transformation `init` ({})!'
                         .format(n_components, init.shape[0]),
                         lmnn.fit, X, y)

    # n_components > X.shape[1]
    n_components = X.shape[1] + 2
    lmnn = LargeMarginNearestNeighbor(init=init, n_components=n_components)
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) cannot be greater '
                         'than the given data dimensionality ({})!'
                         .format(n_components, X.shape[1]),
                         lmnn.fit, X, y)

    # n_components < X.shape[1]
    lmnn = LargeMarginNearestNeighbor(n_components=2, init='identity')
    lmnn.fit(X, y)


def test_init_transformation():
    X, y = datasets.make_classification(n_samples=30, n_features=5,
                                        n_redundant=0, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y)

    # Start learning from scratch
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, init='identity')
    lmnn.fit(X_train, y_train)

    # Initialize with PCA
    lmnn_pca = LargeMarginNearestNeighbor(n_neighbors=3, init='pca')
    lmnn_pca.fit(X_train, y_train)

    # Not always True
    # assert_true(lmnn_pca.n_iter_ <= lmnn.n_iter_)

    init = np.random.rand(X.shape[1], X.shape[1])
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, init=init)
    lmnn.fit(X_train, y_train)

    # init.shape[1] must match X.shape[1]
    init = np.random.rand(X.shape[1], X.shape[1] + 1)
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, init=init)
    assert_raise_message(ValueError,
                         'The input dimensionality ({}) of the given '
                         'linear transformation `init` must match the '
                         'dimensionality of the given inputs `X` ({}).'
                         .format(init.shape[1], X.shape[1]),
                         lmnn.fit, X_train, y_train)

    # init.shape[0] must be <= init.shape[1]
    init = np.random.rand(X.shape[1] + 1, X.shape[1])
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, init=init)
    assert_raise_message(ValueError,
                         'The output dimensionality ({}) of the given '
                         'linear transformation `init` cannot be '
                         'greater than its input dimensionality ({}).'
                         .format(init.shape[0], init.shape[1]),
                         lmnn.fit, X_train, y_train)

    # init.shape[0] must match n_components
    init = np.random.rand(X.shape[1], X.shape[1])
    n_components = X.shape[1] - 2
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, init=init,
                                      n_components=n_components)
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) does not match '
                         'the output dimensionality of the given '
                         'linear transformation `init` ({})!'
                         .format(n_components, init.shape[0]),
                         lmnn.fit, X_train, y_train)


def test_max_impostors():
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, max_impostors=1,
                                      impostor_store='list')
    lmnn.fit(iris_data, iris_target)

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, max_impostors=1,
                                      impostor_store='sparse')
    lmnn.fit(iris_data, iris_target)


def test_neighbors_params():
    from scipy.spatial.distance import hamming

    params = {'algorithm': 'brute', 'metric': hamming}
    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, neighbors_params=params)
    lmnn.fit(iris_data, iris_target)
    components_hamming = lmnn.components_

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3)
    lmnn.fit(iris_data, iris_target)
    components_euclidean = lmnn.components_

    assert_false(np.allclose(components_hamming, components_euclidean))


def test_impostor_store():
    X = iris_data
    y = iris_target
    n_samples, n_features = X.shape
    train_test_boundary = int(n_samples * 0.8)
    train = np.arange(0, train_test_boundary)
    test = np.arange(train_test_boundary, n_samples)
    X_train, y_train, X_test, y_test = X[train], y[train], X[test], y[test]

    k = 3
    lmnn = LargeMarginNearestNeighbor(n_neighbors=k, impostor_store='list')
    lmnn.fit(X_train, y_train)
    knn = KNeighborsClassifier(n_neighbors=k)
    knn.fit(lmnn.fit_transform(X_train, y_train), y_train)
    acc_sparse = knn.score(lmnn.transform(X_test), y_test)

    lmnn = LargeMarginNearestNeighbor(n_neighbors=k, impostor_store='sparse')
    lmnn.fit(X_train, y_train)
    knn.fit(lmnn.fit_transform(X_train, y_train), y_train)
    acc_dense = knn.score(lmnn.transform(X_test), y_test)

    err_msg = 'Toggling `impostor_store` results in different accuracy.'
    assert_equal(acc_dense, acc_sparse, msg=err_msg)


def test_warm_start():
    # A 1-iteration second fit on same data should give almost same result
    # with warm starting, and quite different result without warm starting.

    X, y = datasets.make_classification(n_samples=30, n_features=5,
                                        n_redundant=0, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    n_iter = 10

    lmnn_warm = LargeMarginNearestNeighbor(n_neighbors=3, warm_start=True,
                                           max_iter=n_iter, random_state=0)
    lmnn_warm.fit(X_train, y_train)
    transformation_warm = lmnn_warm.components_
    lmnn_warm.max_iter = 1
    lmnn_warm.fit(X_train, y_train)
    transformation_warm_plus_one = lmnn_warm.components_

    lmnn_cold = LargeMarginNearestNeighbor(n_neighbors=3, warm_start=False,
                                           max_iter=n_iter, random_state=0)
    lmnn_cold.fit(X_train, y_train)
    transformation_cold = lmnn_cold.components_
    lmnn_cold.max_iter = 1
    lmnn_cold.fit(X_train, y_train)
    transformation_cold_plus_one = lmnn_cold.components_

    diff_warm = np.sum(np.abs(transformation_warm_plus_one -
                              transformation_warm))
    diff_cold = np.sum(np.abs(transformation_cold_plus_one -
                              transformation_cold))

    assert_true(diff_warm < 2.0,
                "Transformer changed significantly after one iteration even "
                "though it was warm-started.")

    assert_true(diff_cold > diff_warm,
                "Cold-started transformer changed less significantly than "
                "warm-started transformer after one iteration.")


def test_warm_start_diff_inputs():
    X, y = datasets.make_classification(n_samples=30, n_features=5,
                                        n_classes=4, n_redundant=0,
                                        n_informative=5, random_state=0)

    lmnn = LargeMarginNearestNeighbor(warm_start=True, max_iter=5)
    lmnn.fit(X, y)

    X_less_features, y = \
        datasets.make_classification(n_samples=30, n_features=4, n_classes=4,
                                     n_redundant=0, n_informative=4,
                                     random_state=0)
    assert_raise_message(ValueError,
                         'The new inputs dimensionality ({}) does not '
                         'match the input dimensionality of the '
                         'previously learned transformation ({}).'
                         .format(X_less_features.shape[1],
                                 lmnn.components_.shape[1]),
                         lmnn.fit, X_less_features, y)


def test_verbose():
    # assert there is proper output when verbose = 1
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, verbose=1)
    try:
        lmnn.fit(iris_data, iris_target)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # check output
    assert("[LargeMarginNearestNeighbor]" in out)
    assert("Finding principal components" in out)
    assert ("Finding the target neighbors" in out)
    assert ("Computing static part of the gradient" in out)
    assert ("Finding principal components" in out)
    assert ("Training took" in out)

    # assert by default there is no output (verbose=0)
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3)
    try:
        lmnn.fit(iris_data, iris_target)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # check output
    assert(out == '')


def test_random_state():
    """Assert that when having more than max_impostors (forcing sampling),
    the same impostors will be sampled given the same random_state and
    different impostors will be sampled given a different random_state
    leading to a different transformation"""

    X = iris_data
    y = iris_target

    # Use init='identity' to ensure reproducibility
    params = {'n_neighbors': 3, 'max_impostors': 5, 'random_state': 1,
              'max_iter': 10, 'init': 'identity'}

    lmnn = LargeMarginNearestNeighbor(**params)
    lmnn.fit(X, y)
    transformation_1 = lmnn.components_

    lmnn = LargeMarginNearestNeighbor(**params)
    lmnn.fit(X, y)
    transformation_2 = lmnn.components_

    # This assertion fails on 32bit systems if init='pca'
    assert_allclose(transformation_1, transformation_2)

    params['random_state'] = 2
    lmnn = LargeMarginNearestNeighbor(**params)
    lmnn.fit(X, y)
    transformation_3 = lmnn.components_

    assert_false(np.allclose(transformation_2, transformation_3))


def test_singleton_class():
    X = iris_data
    y = iris_target
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, stratify=y)

    # one singleton class
    singleton_class = 1
    ind_singleton, = np.where(y_tr == singleton_class)
    y_tr[ind_singleton] = 2
    y_tr[ind_singleton[0]] = singleton_class

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, max_iter=30)
    lmnn.fit(X_tr, y_tr)

    # One non-singleton class
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, stratify=y)
    ind_1, = np.where(y_tr == 1)
    ind_2, = np.where(y_tr == 2)
    y_tr[ind_1] = 0
    y_tr[ind_1[0]] = 1
    y_tr[ind_2] = 0
    y_tr[ind_2[0]] = 2

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, max_iter=30)
    assert_raise_message(ValueError,
                         'LargeMarginNearestNeighbor needs at least 2 '
                         'non-singleton classes, got 1.',
                         lmnn.fit, X_tr, y_tr)


def test_callback():

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, callback='my_cb')
    assert_raise_message(ValueError,
                         '`callback` is not callable.',
                         lmnn.fit, iris_data, iris_target)

    max_iter = 10

    def my_cb(transformation, n_iter):
        rem_iter = max_iter - n_iter
        print('{} iterations remaining...'.format(rem_iter))

    # assert that my_cb is called
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, callback=my_cb,
                                      max_iter=max_iter, verbose=1)
    try:
        lmnn.fit(iris_data, iris_target)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # check output
    assert('{} iterations remaining...'.format(max_iter-1) in out)


def test_terminate_early():
    X = iris_data
    y = iris_target
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, max_iter=5)
    lmnn.fit(X_train, y_train)


def test_store_opt_result():
    X = iris_data
    y = iris_target
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    lmnn = LargeMarginNearestNeighbor(n_neighbors=3, max_iter=5,
                                      store_opt_result=True)
    lmnn.fit(X_train, y_train)
    transformation = lmnn.opt_result_.x
    assert_equal(transformation.size, X.shape[1]**2)


def test_paired_distances_blockwise():
    n, d = 10000, 100  # 4 or 8 MiB
    X = rng.rand(n, d)
    ind_a = rng.permutation(n)
    ind_b = rng.permutation(n)

    distances = paired_euclidean_distances(X[ind_a], X[ind_b])
    distances_blockwise = _paired_distances_blockwise(
        X, ind_a, ind_b, squared=False, block_size=1)
    assert_array_equal(distances, distances_blockwise)


def test_euclidean_distances_without_checks():
    X = rng.rand(100, 20)
    Y = rng.rand(50, 20)

    # 2 matrices with no precomputed norms
    distances1 = euclidean_distances(X, Y)
    distances2 = _euclidean_distances_without_checks(X, Y)

    assert_array_equal(distances1, distances2)

    # 1 matrix with itself with squared row_norms precomputed and transposed
    XX = row_norms(X, squared=True)[np.newaxis, :]
    distances1 = euclidean_distances(X, X_norm_squared=XX)
    distances2 = _euclidean_distances_without_checks(X, X_norm_squared=XX)

    assert_array_equal(distances1, distances2)


def test_pipeline_equivalency():
    X = iris_data
    y = iris_target
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    # Use init='identity' to ensure reproducibility
    lmnn_params = dict(n_neighbors=3, max_iter=10, init='identity',
                       random_state=42)
    n_neighbors = 3

    lmnn = LargeMarginNearestNeighbor(**lmnn_params)
    lmnn.fit(X_train, y_train)

    lmnn_pipe = make_lmnn_pipeline(**lmnn_params)
    lmnn_pipe.fit(X_train, y_train)

    pipe_transformation = lmnn_pipe.named_steps.lmnn.components_
    assert_array_almost_equal(lmnn.components_, pipe_transformation)

    knn = KNeighborsClassifier(n_neighbors=n_neighbors)
    knn.fit(lmnn.transform(X_train), y_train)
    score = knn.score(lmnn.transform(X_test), y_test)

    score_pipe = lmnn_pipe.score(X_test, y_test)

    assert_equal(score, score_pipe)
