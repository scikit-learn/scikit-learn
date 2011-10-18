import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal, \
                          assert_equal

from .. import DictionaryLearning, MiniBatchDictionaryLearning, \
               dict_learning_online
from ..dict_learning import sparse_encode, sparse_encode_parallel

rng = np.random.RandomState(0)
n_samples, n_features = 10, 8
X = rng.randn(n_samples, n_features)


def test_dict_learning_shapes():
    n_atoms = 5
    dico = DictionaryLearning(n_atoms).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_overcomplete():
    n_atoms = 12
    X = rng.randn(n_samples, n_features)
    dico = DictionaryLearning(n_atoms).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_reconstruction():
    n_atoms = 12
    dico = DictionaryLearning(n_atoms, transform_algorithm='omp',
                              transform_alpha=0.001, random_state=0)
    code = dico.fit(X).transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X)

    dico.set_params(transform_algorithm='lasso_lars')
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)

    dico.set_params(transform_algorithm='lars')
    code = dico.transform(X)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)


def test_dict_learning_nonzero_coefs():
    n_atoms = 4
    dico = DictionaryLearning(n_atoms, transform_algorithm='lars',
                              transform_n_nonzero_coefs=3, random_state=0)
    code = dico.fit(X).transform(X[0])
    assert len(np.flatnonzero(code)) == 3

    dico.set_params(transform_algorithm='omp')
    code = dico.transform(X[0])
    assert len(np.flatnonzero(code)) == 3


def test_dict_learning_split():
    n_atoms = 5
    dico = DictionaryLearning(n_atoms, transform_algorithm='threshold')
    code = dico.fit(X).transform(X)
    dico.split_sign = True
    split_code = dico.transform(X)

    assert_array_equal(split_code[:, :n_atoms] - split_code[:, n_atoms:], code)


def test_dict_learning_online_shapes():
    rng = np.random.RandomState(0)
    X = rng.randn(12, 10)
    dictionaryT, codeT = dict_learning_online(X.T, n_atoms=8, alpha=1,
                                              random_state=rng)
    assert_equal(codeT.shape, (8, 12))
    assert_equal(dictionaryT.shape, (10, 8))
    assert_equal(np.dot(codeT.T, dictionaryT.T).shape, X.shape)


def test_dict_learning_online_estimator_shapes():
    n_atoms = 5
    dico = MiniBatchDictionaryLearning(n_atoms, n_iter=20).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_online_overcomplete():
    n_atoms = 12
    dico = MiniBatchDictionaryLearning(n_atoms, n_iter=20).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_online_initialization():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)
    dico = MiniBatchDictionaryLearning(n_atoms, n_iter=0, dict_init=V).fit(X)
    assert_array_equal(dico.components_, V)


def test_dict_learning_online_partial_fit():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)  # random init
    rng1 = np.random.RandomState(0)
    rng2 = np.random.RandomState(0)
    dico1 = MiniBatchDictionaryLearning(n_atoms, n_iter=10, chunk_size=1,
                                        shuffle=False, dict_init=V,
                                        transform_algorithm='threshold',
                                        random_state=rng1).fit(X)
    dico2 = MiniBatchDictionaryLearning(n_atoms, n_iter=1, dict_init=V,
                                        transform_algorithm='threshold',
                                        random_state=rng2)
    for ii, sample in enumerate(X):
        dico2.partial_fit(sample, iter_offset=ii * dico2.n_iter)

    assert_array_equal(dico1.components_, dico2.components_)


def test_sparse_code():
    rng = np.random.RandomState(0)
    dictionary = rng.randn(10, 3)
    real_code = np.zeros((3, 5))
    real_code.ravel()[rng.randint(15, size=6)] = 1.0
    Y = np.dot(dictionary, real_code)
    est_code_1 = sparse_encode(dictionary, Y, alpha=1.0)
    est_code_2 = sparse_encode_parallel(dictionary, Y, alpha=1.0)
    assert_equal(est_code_1.shape, real_code.shape)
    assert_equal(est_code_1, est_code_2)
    assert_equal(est_code_1.nonzero(), real_code.nonzero())
