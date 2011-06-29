import numpy as np
from numpy.testing import assert_array_almost_equal
from nose.plugins.skip import SkipTest

from .. import DictionaryLearning

def test_dict_learning_shapes():
    n_samples, n_features = 10, 8
    n_atoms = 5
    X = np.random.randn(n_samples, n_features)
    dico = DictionaryLearning(n_atoms).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)


def test_dict_learning_overcomplete():
    raise SkipTest
    n_samples, n_features = 10, 8
    n_atoms = 12
    X = np.random.randn(n_samples, n_features)
    dico = DictionaryLearning(n_atoms).fit(X)
    assert dico.components_.shape == (n_atoms, n_features)

def test_dict_learning_reconstruction():
    n_samples, n_features = 10, 8
    n_atoms = 5
    U = np.zeros((n_samples, n_atoms)).ravel()
    U[np.random.randint(len(U), size=len(U) / n_samples)] = 1.0
    U = U.reshape((n_samples, n_atoms))
    V = np.random.randn(n_atoms, n_features)

    X = np.dot(U, V)
    dico = DictionaryLearning(n_atoms)
    code = dico.fit(X).transform(X, method='lars', alpha=0.01)

    assert_array_almost_equal(np.dot(code, dico.components_), X, decimal=1)
