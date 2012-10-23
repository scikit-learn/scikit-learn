import numpy as np

from numpy.testing import assert_array_almost_equal
from unittest import TestCase

from sklearn.datasets import load_digits
from sklearn.neural_networks import RestrictedBolzmannMachine

np.seterr(all='warn')


class TestBaseRBM(TestCase):

    def setUp(self):
        self.prng = np.random.RandomState(9)

    def test_fit(self):
        X = load_digits().data
        
        X -= X.min()
        X /= X.max()
        
        rbm = RestrictedBolzmannMachine(n_components=64, learning_rate=0.1,
            n_particles=10, n_iter=10, random_state=self.prng)
        
        rbm.fit(X)
        
        assert_array_almost_equal(-20., rbm.pseudo_likelihood(X).mean(),
            decimal=0)
        
        
