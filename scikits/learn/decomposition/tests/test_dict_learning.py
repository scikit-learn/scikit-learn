import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from nose.plugins.skip import SkipTest

from ...datasets import make_sparse_coded_signal
from .. import DictionaryLearning, DictionaryLearningOnline

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
    dico = DictionaryLearning(n_atoms, transform_algorithm='omp')
    code = dico.fit(X).transform(X, eps=0.01)
    assert_array_almost_equal(np.dot(code, dico.components_), X)

    dico.transform_algorithm = 'lasso_lars'
    code = dico.transform(X, alpha=0.001)
    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=2)


def test_dict_learning_split():
    n_atoms = 5
    dico = DictionaryLearning(n_atoms, transform_algorithm='threshold')
    code = dico.fit(X).transform(X, alpha=1)
    dico.split_sign = True
    split_code = dico.transform(X, alpha=1)

    assert_array_equal(split_code[:, :n_atoms] - split_code[:, n_atoms:], code)


def test_dict_learning_online_shapes():
    n_atoms = 5
    dico = DictionaryLearningOnline(n_atoms, n_iter=20).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_online_overcomplete():
    n_atoms = 12
    dico = DictionaryLearningOnline(n_atoms, n_iter=20).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_online_initialization():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)
    dico = DictionaryLearningOnline(n_atoms, n_iter=0, dict_init=V).fit(X)
    assert_array_equal(dico.components_, V)


def test_dict_learning_online_partial_fit():
    n_atoms = 12
    V = rng.randn(n_atoms, n_features)  # random init
    rng1 = np.random.RandomState(0)
    rng2 = np.random.RandomState(0)
    dico1 = DictionaryLearningOnline(n_atoms, n_iter=10, chunk_size=1,
                                     shuffle=False, dict_init=V,
                                     transform_algorithm='threshold',
                                     random_state=rng1).fit(X)
    dico2 = DictionaryLearningOnline(n_atoms, n_iter=1, dict_init=V,
                                     transform_algorithm='threshold',
                                     random_state=rng2)
    for ii, sample in enumerate(X):
        dico2.partial_fit(sample, iter_offset=ii * dico2.n_iter)

    assert_array_equal(dico1.components_, dico2.components_)
