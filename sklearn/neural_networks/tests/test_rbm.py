import numpy as np

from numpy.testing import assert_almost_equal, assert_array_equal

from sklearn.datasets import load_digits
from sklearn.neural_networks import RestrictedBolzmannMachine

np.seterr(all='warn')


def test_fit():
    X = load_digits().data

    X -= X.min()
    X /= X.max()

    rbm = RestrictedBolzmannMachine(n_components=64, learning_rate=0.1,
                                    n_particles=10, n_iter=7, random_state=9)
    rbm.fit(X)

    assert_almost_equal(rbm.pseudo_likelihood(X).mean(), -20., decimal=0)
