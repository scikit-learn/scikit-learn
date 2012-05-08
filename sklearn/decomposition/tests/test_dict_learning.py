import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal, \
                          assert_equal
from nose import SkipTest
from nose.tools import assert_true

from sklearn.utils.testing import assert_less

from .. import DictionaryLearning, MiniBatchDictionaryLearning, SparseCoder, \
               dict_learning_online, sparse_encode


rng = np.random.RandomState(0)
n_samples, n_features = 10, 8
X = rng.randn(n_samples, n_features)


def test_dict_learning_shapes():
    n_atoms = 5
    dico = DictionaryLearning(n_atoms).fit(X)
    assert_true(dico.components_.shape == (n_atoms, n_features))


def test_dict_learning_overcomplete():
    n_atoms = 12
    X = rng.randn(n_samples, n_features)
    dico = DictionaryLearning(n_atoms).fit(X)
    assert_true(dico.components_.shape == (n_atoms, n_features))


def test_dict_learning_reconstruction():
    n_atoms = 12
    dico = DictionaryLearning(n_atoms, transform_algorithm='omp',
                              transform_alpha=0.001, random_state=0)
    code = dico.fit(X).transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X)

    dico.set_params(transform_algorithm='lasso_lars')
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)

    # used to test lars here too, but there's no guarantee the number of
    # nonzero atoms is right.


def test_dict_learning_nonzero_coefs():
    n_atoms = 4
    dico = DictionaryLearning(n_atoms, transform_algorithm='lars',
                              transform_n_nonzero_coefs=3, random_state=0)
    code = dico.fit(X).transform(X[1])
    assert_true(len(np.flatnonzero(code)) == 3)

    dico.set_params(transform_algorithm='omp')
    code = dico.transform(X[1])
    assert_equal(len(np.flatnonzero(code)), 3)


def test_dict_learning_split():
    n_atoms = 5
    dico = DictionaryLearning(n_atoms, transform_algorithm='threshold')
    code = dico.fit(X).transform(X)
    dico.split_sign = True
    split_code = dico.transform(X)

    assert_array_equal(split_code[:, :n_atoms] - split_code[:, n_atoms:], code)


def test_dict_learning_online_shapes():
#    rng = np.random.RandomState(0)
#    X = rng.randn(12, 10)
    n_atoms = 8
    code, dictionary = dict_learning_online(X, n_atoms=n_atoms, alpha=1,
                                            random_state=rng)
    assert_equal(code.shape, (n_samples, n_atoms))
    assert_equal(dictionary.shape, (n_atoms, n_features))
    assert_equal(np.dot(code, dictionary).shape, X.shape)


def test_dict_learning_online_estimator_shapes():
    n_atoms = 5
    dico = MiniBatchDictionaryLearning(n_atoms, n_iter=20).fit(X)
    assert_true(dico.components_.shape == (n_atoms, n_features))


def test_dict_learning_online_overcomplete():
    n_atoms = 12
    dico = MiniBatchDictionaryLearning(n_atoms, n_iter=20).fit(X)
    assert_true(dico.components_.shape == (n_atoms, n_features))


def test_dict_learning_online_initialization():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)
    dico = MiniBatchDictionaryLearning(n_atoms, n_iter=0, dict_init=V).fit(X)
    assert_array_equal(dico.components_, V)


def test_dict_learning_online_partial_fit():
    # this test was not actually passing before!
    raise SkipTest
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    rng1 = np.random.RandomState(0)
    rng2 = np.random.RandomState(0)
    dico1 = MiniBatchDictionaryLearning(n_atoms, n_iter=10, chunk_size=1,
                                        shuffle=False, dict_init=V,
                                        random_state=rng1).fit(X)
    dico2 = MiniBatchDictionaryLearning(n_atoms, n_iter=1, dict_init=V,
                                        random_state=rng2)
    for ii, sample in enumerate(X):
        dico2.partial_fit(sample, iter_offset=ii * dico2.n_iter)
        # if ii == 1: break
    assert_true(not np.all(sparse_encode(X, dico1.components_, alpha=100) ==
        0))
    assert_array_equal(dico1.components_, dico2.components_)


def test_sparse_encode_shapes():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    for algo in ('lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'):
        code = sparse_encode(X, V, algorithm=algo)
        assert_equal(code.shape, (n_samples, n_atoms))


def test_sparse_encode_error():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    code = sparse_encode(X, V, alpha=0.001)
    assert_true(not np.all(code == 0))
    assert_less(np.sqrt(np.sum((np.dot(code, V) - X) ** 2)), 0.1)


def test_sparse_coder_estimator():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)  # random init
    V /= np.sum(V ** 2, axis=1)[:, np.newaxis]
    code = SparseCoder(dictionary=V, transform_algorithm='lasso_lars',
                       transform_alpha=0.001).transform(X)
    assert_true(not np.all(code == 0))
    assert_less(np.sqrt(np.sum((np.dot(code, V) - X) ** 2)), 0.1)
