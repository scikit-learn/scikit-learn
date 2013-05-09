import numpy as np

from numpy.testing import assert_almost_equal, assert_array_equal

from sklearn.base import clone
from sklearn.datasets import load_digits
from sklearn.neural_network import BernoulliRBM

np.seterr(all='warn')

Xdigits = load_digits().data
Xdigits -= Xdigits.min()
Xdigits /= Xdigits.max()


def test_fit():
    X = Xdigits.copy()

    rbm = BernoulliRBM(n_components=64, learning_rate=0.1,
                       batch_size=10, n_iter=7, random_state=9)
    rbm.fit(X)

    assert_almost_equal(rbm.pseudo_likelihood(X).mean(), -20., decimal=0)

    # in-place tricks shouldn't have modified X
    assert_array_equal(X, Xdigits)


def test_fit_transform():
    """Check proper implementation of fit_transform"""
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=16, batch_size=5,
                        n_iter=5, random_state=42)
    rbm2 = clone(rbm1)

    Xt1 = rbm1.fit(X).transform(X)
    Xt2 = rbm2.fit_transform(X)

    assert_array_equal(Xt1, Xt2)


def test_transform():
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=16, batch_size=5,
                        n_iter=5, random_state=42)
    rbm1.fit(X)
    
    Xt1 = rbm1.transform(X)
    Xt2 = rbm1.mean_hiddens(X)

    assert_array_equal(Xt1, Xt2)


def test_sample_hiddens():
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=2, batch_size=5,
                        n_iter=5, random_state=42)
    rbm1.fit(X)
    
    h = rbm1.mean_hiddens(X[0])
    hs = np.mean([rbm1.sample_hiddens(X[0]) for i in range(100)], 0)
    
    assert_almost_equal(h, hs, decimal=1)


def test_gibbs():
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=2, batch_size=5,
                        n_iter=5, random_state=42)
    rbm1.fit(X)
    
    Xt1 = np.mean([rbm1.gibbs(X[0]) for i in range(100)], 0)
    Xt2 = np.mean([rbm1.sample_visibles(rbm1.sample_hiddens(X[0]))
        for i in range(100)], 0)
    
    assert_almost_equal(Xt1, Xt2, decimal=1)
