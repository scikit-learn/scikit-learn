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

    assert_almost_equal(rbm.pseudo_likelihood(X).mean(), -21., decimal=0)

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
    Xt2 = rbm1._mean_hiddens(X)

    assert_array_equal(Xt1, Xt2)


def test_sample_hiddens():
    rng = np.random.RandomState(0)
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=2, batch_size=5,
                        n_iter=5, random_state=42)
    rbm1.fit(X)

    h = rbm1._mean_hiddens(X[0])
    hs = np.mean([rbm1._sample_hiddens(X[0], rng) for i in range(100)], 0)

    assert_almost_equal(h, hs, decimal=1)


def test_fit_gibbs():
    """
    Gibbs on the RBM hidden layer should be able to recreate [[0], [1]]
    from the same input
    """
    rng = np.random.RandomState(42)
    X = np.array([[0.], [1.]])
    rbm1 = BernoulliRBM(n_components=2, batch_size=2,
                        n_iter=42, random_state=rng)
                        # you need that much iters
    rbm1.fit(X)
    assert_almost_equal(rbm1.components_,
                        np.array([[0.02649814], [0.02009084]]), decimal=4)
    assert_almost_equal(rbm1.gibbs(X), X)


def test_gibbs():
    rng = np.random.RandomState(42)
    X = Xdigits[:100]
    rbm1 = BernoulliRBM(n_components=2, batch_size=5,
                        n_iter=5, random_state=rng)
    rbm1.fit(X)

    Xt1 = np.mean([rbm1.gibbs(X[0]) for i in range(200)], 0)
    Xt2 = np.mean([rbm1._sample_visibles(rbm1._sample_hiddens(X[0], rng), rng)
                   for i in range(1000)], 0)

    assert_almost_equal(Xt1, Xt2, decimal=1)
