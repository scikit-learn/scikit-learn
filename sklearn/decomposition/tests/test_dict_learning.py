from __future__ import division
import pytest

import numpy as np
import itertools

from sklearn.exceptions import ConvergenceWarning

from sklearn.utils import check_array

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import TempMemmap

from sklearn.decomposition import DictionaryLearning
from sklearn.decomposition import MiniBatchDictionaryLearning
from sklearn.decomposition import SparseCoder
from sklearn.decomposition import dict_learning_online
from sklearn.decomposition import sparse_encode


rng_global = np.random.RandomState(0)
n_samples, n_features = 10, 8
X = rng_global.randn(n_samples, n_features)


def test_sparse_encode_shapes_omp():
    rng = np.random.RandomState(0)
    algorithms = ['omp', 'lasso_lars', 'lasso_cd', 'lars', 'threshold']
    for n_components, n_samples in itertools.product([1, 5], [1, 9]):
        X_ = rng.randn(n_samples, n_features)
        dictionary = rng.randn(n_components, n_features)
        for algorithm, n_jobs in itertools.product(algorithms, [1, 3]):
            code = sparse_encode(X_, dictionary, algorithm=algorithm,
                                 n_jobs=n_jobs)
            assert_equal(code.shape, (n_samples, n_components))


def test_dict_learning_shapes():
    n_components = 5
    dico = DictionaryLearning(n_components, random_state=0).fit(X)
    assert_equal(dico.components_.shape, (n_components, n_features))

    n_components = 1
    dico = DictionaryLearning(n_components, random_state=0).fit(X)
    assert_equal(dico.components_.shape, (n_components, n_features))
    assert_equal(dico.transform(X).shape, (X.shape[0], n_components))


def test_dict_learning_overcomplete():
    n_components = 12
    dico = DictionaryLearning(n_components, random_state=0).fit(X)
    assert dico.components_.shape == (n_components, n_features)


def test_max_iter():
    def ricker_function(resolution, center, width):
        """Discrete sub-sampled Ricker (Mexican hat) wavelet"""
        x = np.linspace(0, resolution - 1, resolution)
        x = ((2 / ((np.sqrt(3 * width) * np.pi ** 1 / 4)))
             * (1 - ((x - center) ** 2 / width ** 2))
             * np.exp((-(x - center) ** 2) / (2 * width ** 2)))
        return x

    def ricker_matrix(width, resolution, n_components):
        """Dictionary of Ricker (Mexican hat) wavelets"""
        centers = np.linspace(0, resolution - 1, n_components)
        D = np.empty((n_components, resolution))
        for i, center in enumerate(centers):
            D[i] = ricker_function(resolution, center, width)
        D /= np.sqrt(np.sum(D ** 2, axis=1))[:, np.newaxis]
        return D

    transform_algorithm = 'lasso_cd'
    resolution = 1024
    subsampling = 3  # subsampling factor
    n_components = resolution // subsampling

    # Compute a wavelet dictionary
    D_multi = np.r_[tuple(ricker_matrix(width=w, resolution=resolution,
                          n_components=n_components // 5)
                          for w in (10, 50, 100, 500, 1000))]

    X = np.linspace(0, resolution - 1, resolution)
    first_quarter = X < resolution / 4
    X[first_quarter] = 3.
    X[np.logical_not(first_quarter)] = -1.
    X = X.reshape(1, -1)

    # check that the underlying model fails to converge
    with pytest.warns(ConvergenceWarning):
        model = SparseCoder(D_multi, transform_algorithm=transform_algorithm,
                            transform_max_iter=1)
        model.fit_transform(X)

    # check that the underlying model converges w/o warnings
    with pytest.warns(None) as record:
        model = SparseCoder(D_multi, transform_algorithm=transform_algorithm,
                            transform_max_iter=2000)
        model.fit_transform(X)
    assert not record.list


# positive lars deprecated 0.22
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.parametrize("transform_algorithm", [
    "lasso_lars",
    "lasso_cd",
    "lars",
    "threshold",
])
@pytest.mark.parametrize("positive_code", [
    False,
    True,
])
@pytest.mark.parametrize("positive_dict", [
    False,
    True,
])
def test_dict_learning_positivity(transform_algorithm,
                                  positive_code,
                                  positive_dict):
    n_components = 5
    dico = DictionaryLearning(
        n_components, transform_algorithm=transform_algorithm, random_state=0,
        positive_code=positive_code, positive_dict=positive_dict).fit(X)
    code = dico.transform(X)
    if positive_dict:
        assert (dico.components_ >= 0).all()
    else:
        assert (dico.components_ < 0).any()
    if positive_code:
        assert (code >= 0).all()
    else:
        assert (code < 0).any()


def test_dict_learning_reconstruction():
    n_components = 12
    dico = DictionaryLearning(n_components, transform_algorithm='omp',
                              transform_alpha=0.001, random_state=0)
    code = dico.fit(X).transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X)

    dico.set_params(transform_algorithm='lasso_lars')
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)

    # used to test lars here too, but there's no guarantee the number of
    # nonzero atoms is right.


def test_dict_learning_reconstruction_parallel():
    # regression test that parallel reconstruction works with n_jobs=-1
    n_components = 12
    dico = DictionaryLearning(n_components, transform_algorithm='omp',
                              transform_alpha=0.001, random_state=0, n_jobs=-1)
    code = dico.fit(X).transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X)

    dico.set_params(transform_algorithm='lasso_lars')
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)


def test_dict_learning_lassocd_readonly_data():
    n_components = 12
    with TempMemmap(X) as X_read_only:
        dico = DictionaryLearning(n_components, transform_algorithm='lasso_cd',
                                  transform_alpha=0.001, random_state=0,
                                  n_jobs=-1)
        with ignore_warnings(category=ConvergenceWarning):
            code = dico.fit(X_read_only).transform(X_read_only)
        assert_array_almost_equal(np.dot(code, dico.components_), X_read_only,
                                  decimal=2)


def test_dict_learning_nonzero_coefs():
    n_components = 4
    dico = DictionaryLearning(n_components, transform_algorithm='lars',
                              transform_n_nonzero_coefs=3, random_state=0)
    code = dico.fit(X).transform(X[np.newaxis, 1])
    assert len(np.flatnonzero(code)) == 3

    dico.set_params(transform_algorithm='omp')
    code = dico.transform(X[np.newaxis, 1])
    assert_equal(len(np.flatnonzero(code)), 3)


def test_dict_learning_unknown_fit_algorithm():
    n_components = 5
    dico = DictionaryLearning(n_components, fit_algorithm='<unknown>')
    assert_raises(ValueError, dico.fit, X)


def test_dict_learning_split():
    n_components = 5
    dico = DictionaryLearning(n_components, transform_algorithm='threshold',
                              random_state=0)
    code = dico.fit(X).transform(X)
    dico.split_sign = True
    split_code = dico.transform(X)

    assert_array_almost_equal(split_code[:, :n_components] -
                              split_code[:, n_components:], code)


def test_dict_learning_online_shapes():
    rng = np.random.RandomState(0)
    n_components = 8
    code, dictionary = dict_learning_online(X, n_components=n_components,
                                            alpha=1, random_state=rng)
    assert_equal(code.shape, (n_samples, n_components))
    assert_equal(dictionary.shape, (n_components, n_features))
    assert_equal(np.dot(code, dictionary).shape, X.shape)


# positive lars deprecated 0.22
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.parametrize("transform_algorithm", [
    "lasso_lars",
    "lasso_cd",
    "lars",
    "threshold",
])
@pytest.mark.parametrize("positive_code", [
    False,
    True,
])
@pytest.mark.parametrize("positive_dict", [
    False,
    True,
])
def test_dict_learning_online_positivity(transform_algorithm,
                                         positive_code,
                                         positive_dict):
    rng = np.random.RandomState(0)
    n_components = 8

    dico = MiniBatchDictionaryLearning(
        n_components, transform_algorithm=transform_algorithm, random_state=0,
        positive_code=positive_code, positive_dict=positive_dict).fit(X)
    code = dico.transform(X)
    if positive_dict:
        assert (dico.components_ >= 0).all()
    else:
        assert (dico.components_ < 0).any()
    if positive_code:
        assert (code >= 0).all()
    else:
        assert (code < 0).any()

    code, dictionary = dict_learning_online(X, n_components=n_components,
                                            alpha=1, random_state=rng,
                                            positive_dict=positive_dict,
                                            positive_code=positive_code)
    if positive_dict:
        assert (dictionary >= 0).all()
    else:
        assert (dictionary < 0).any()
    if positive_code:
        assert (code >= 0).all()
    else:
        assert (code < 0).any()


def test_dict_learning_online_verbosity():
    n_components = 5
    # test verbosity
    from sklearn.externals.six.moves import cStringIO as StringIO
    import sys

    old_stdout = sys.stdout
    try:
        sys.stdout = StringIO()
        dico = MiniBatchDictionaryLearning(n_components, n_iter=20, verbose=1,
                                           random_state=0)
        dico.fit(X)
        dico = MiniBatchDictionaryLearning(n_components, n_iter=20, verbose=2,
                                           random_state=0)
        dico.fit(X)
        dict_learning_online(X, n_components=n_components, alpha=1, verbose=1,
                             random_state=0)
        dict_learning_online(X, n_components=n_components, alpha=1, verbose=2,
                             random_state=0)
    finally:
        sys.stdout = old_stdout

    assert dico.components_.shape == (n_components, n_features)


def test_dict_learning_online_estimator_shapes():
    n_components = 5
    dico = MiniBatchDictionaryLearning(n_components, n_iter=20, random_state=0)
    dico.fit(X)
    assert dico.components_.shape == (n_components, n_features)


def test_dict_learning_online_overcomplete():
    n_components = 12
    dico = MiniBatchDictionaryLearning(n_components, n_iter=20,
                                       random_state=0).fit(X)
    assert dico.components_.shape == (n_components, n_features)


def test_dict_learning_online_initialization():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)
    dico = MiniBatchDictionaryLearning(n_components, n_iter=0,
                                       dict_init=V, random_state=0).fit(X)
    assert_array_equal(dico.components_, V)


def test_dict_learning_online_readonly_initialization():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)
    V.setflags(write=False)
    MiniBatchDictionaryLearning(n_components, n_iter=1, dict_init=V,
                                random_state=0, shuffle=False).fit(X)


def test_dict_learning_online_partial_fit():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    dict1 = MiniBatchDictionaryLearning(n_components, n_iter=10 * len(X),
                                        batch_size=1,
                                        alpha=1, shuffle=False, dict_init=V,
                                        random_state=0).fit(X)
    dict2 = MiniBatchDictionaryLearning(n_components, alpha=1,
                                        n_iter=1, dict_init=V,
                                        random_state=0)
    for i in range(10):
        for sample in X:
            dict2.partial_fit(sample[np.newaxis, :])

    assert not np.all(sparse_encode(X, dict1.components_, alpha=1) == 0)
    assert_array_almost_equal(dict1.components_, dict2.components_,
                              decimal=2)


def test_sparse_encode_shapes():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    for algo in ('lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'):
        code = sparse_encode(X, V, algorithm=algo)
        assert_equal(code.shape, (n_samples, n_components))


# positive lars deprecated 0.22
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.parametrize("positive", [
    False,
    True,
])
def test_sparse_encode_positivity(positive):
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    for algo in ('lasso_lars', 'lasso_cd', 'lars', 'threshold'):
        code = sparse_encode(X, V, algorithm=algo, positive=positive)
        if positive:
            assert (code >= 0).all()
        else:
            assert (code < 0).any()

    try:
        sparse_encode(X, V, algorithm='omp', positive=positive)
    except ValueError:
        if not positive:
            raise


def test_sparse_encode_input():
    n_components = 100
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    Xf = check_array(X, order='F')
    for algo in ('lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'):
        a = sparse_encode(X, V, algorithm=algo)
        b = sparse_encode(Xf, V, algorithm=algo)
        assert_array_almost_equal(a, b)


def test_sparse_encode_error():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    code = sparse_encode(X, V, alpha=0.001)
    assert not np.all(code == 0)
    assert_less(np.sqrt(np.sum((np.dot(code, V) - X) ** 2)), 0.1)


def test_sparse_encode_error_default_sparsity():
    rng = np.random.RandomState(0)
    X = rng.randn(100, 64)
    D = rng.randn(2, 64)
    code = ignore_warnings(sparse_encode)(X, D, algorithm='omp',
                                          n_nonzero_coefs=None)
    assert_equal(code.shape, (100, 2))


def test_unknown_method():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    assert_raises(ValueError, sparse_encode, X, V, algorithm="<unknown>")


def test_sparse_coder_estimator():
    n_components = 12
    rng = np.random.RandomState(0)
    V = rng.randn(n_components, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    code = SparseCoder(dictionary=V, transform_algorithm='lasso_lars',
                       transform_alpha=0.001).transform(X)
    assert not np.all(code == 0)
    assert_less(np.sqrt(np.sum((np.dot(code, V) - X) ** 2)), 0.1)


def test_sparse_coder_parallel_mmap():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/5956
    # Test that SparseCoder does not error by passing reading only
    # arrays to child processes

    rng = np.random.RandomState(777)
    n_components, n_features = 40, 64
    init_dict = rng.rand(n_components, n_features)
    # Ensure that `data` is >2M. Joblib memory maps arrays
    # if they are larger than 1MB. The 4 accounts for float32
    # data type
    n_samples = int(2e6) // (4 * n_features)
    data = np.random.rand(n_samples, n_features).astype(np.float32)

    sc = SparseCoder(init_dict, transform_algorithm='omp', n_jobs=2)
    sc.fit_transform(data)
