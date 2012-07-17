"""
Test the fastica algorithm.
"""
import warnings
import numpy as np
from numpy.testing import assert_almost_equal
from nose.tools import assert_true

from scipy import stats
import itertools

from sklearn.utils.testing import assert_less

from .. import FastICA, fastica, PCA
from ..fastica_ import _gs_decorrelation


def center_and_norm(x, axis=-1):
    """ Centers and norms x **in place**

        Parameters
        -----------
        x: ndarray
            Array with an axis of observations (statistical units) measured on
            random variables.
        axis: int, optional
            Axis along which the mean and variance are calculated.
    """
    x = np.rollaxis(x, axis)
    x -= x.mean(axis=0)
    x /= x.std(axis=0)


def test_gs():
    """
    Test gram schmidt orthonormalization
    """
    # generate a random orthogonal  matrix
    rng = np.random.RandomState(0)
    W, _, _ = np.linalg.svd(rng.randn(10, 10))
    w = rng.randn(10)
    _gs_decorrelation(w, W, 10)
    assert_less((w ** 2).sum(), 1.e-10)
    w = rng.randn(10)
    u = _gs_decorrelation(w, W, 5)
    tmp = np.dot(u, W.T)
    assert_less((tmp[:5] ** 2).sum(), 1.e-10)


def test_fastica(add_noise=False):
    """ Test the FastICA algorithm on very simple data.
    """
    # scipy.stats uses the global RNG:
    rng = np.random.RandomState(0)
    n_samples = 1000
    # Generate two sources:
    s1 = (2 * np.sin(np.linspace(0, 100, n_samples)) > 0) - 1
    s2 = stats.t.rvs(1, size=n_samples)
    s = np.c_[s1, s2].T
    center_and_norm(s)
    s1, s2 = s

    # Mixing angle
    phi = 0.6
    mixing = np.array([[np.cos(phi),  np.sin(phi)],
                       [np.sin(phi), -np.cos(phi)]])
    m = np.dot(mixing, s)

    if add_noise:
        m += 0.1 * rng.randn(2, 1000)

    center_and_norm(m)

    algos = ['parallel', 'deflation']
    nls = ['logcosh', 'exp', 'cube']
    whitening = [True, False]
    for algo, nl, whiten in itertools.product(algos, nls, whitening):
        if whiten:
            k_, mixing_, s_ = fastica(m.T, fun=nl, algorithm=algo)
        else:
            X = PCA(n_components=2, whiten=True).fit_transform(m.T)
            k_, mixing_, s_ = fastica(X, fun=nl, algorithm=algo,
                                     whiten=False)
        s_ = s_.T
        # Check that the mixing model described in the docstring holds:
        if whiten:
            assert_almost_equal(s_, np.dot(np.dot(mixing_, k_), m))

        center_and_norm(s_)
        s1_, s2_ = s_
        # Check to see if the sources have been estimated
        # in the wrong order
        if abs(np.dot(s1_, s2)) > abs(np.dot(s1_, s1)):
            s2_, s1_ = s_
        s1_ *= np.sign(np.dot(s1_, s1))
        s2_ *= np.sign(np.dot(s2_, s2))

        # Check that we have estimated the original sources
        if add_noise == False:
            assert_almost_equal(np.dot(s1_, s1) / n_samples, 1, decimal=2)
            assert_almost_equal(np.dot(s2_, s2) / n_samples, 1, decimal=2)
        else:
            assert_almost_equal(np.dot(s1_, s1) / n_samples, 1, decimal=1)
            assert_almost_equal(np.dot(s2_, s2) / n_samples, 1, decimal=1)

    # Test FastICA class
    ica = FastICA(fun=nl, algorithm=algo, random_state=0)
    ica.fit(m.T)
    ica.get_mixing_matrix()
    assert_true(ica.components_.shape == (2, 2))
    assert_true(ica.sources_.shape == (1000, 2))


def test_fastica_nowhiten():
    m = [[0, 1], [1, 0]]
    ica = FastICA(whiten=False, random_state=0)
    ica.fit(m)
    ica.get_mixing_matrix()

    # test for issue #697
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ica = FastICA(n_components=1, whiten=False, random_state=0)
        ica.fit(m)  # should raise warning
        assert_true(len(w) == 1)  # 1 warning should be raised


def test_non_square_fastica(add_noise=False):
    """ Test the FastICA algorithm on very simple data.
    """
    rng = np.random.RandomState(0)

    n_samples = 1000
    # Generate two sources:
    t = np.linspace(0, 100, n_samples)
    s1 = np.sin(t)
    s2 = np.ceil(np.sin(np.pi * t))
    s = np.c_[s1, s2].T
    center_and_norm(s)
    s1, s2 = s

    # Mixing matrix
    mixing = rng.randn(6, 2)
    m = np.dot(mixing, s)

    if add_noise:
        m += 0.1 * rng.randn(6, n_samples)

    center_and_norm(m)

    k_, mixing_, s_ = fastica(m.T, n_components=2, random_state=rng)
    s_ = s_.T

    # Check that the mixing model described in the docstring holds:
    assert_almost_equal(s_, np.dot(np.dot(mixing_, k_), m))

    center_and_norm(s_)
    s1_, s2_ = s_
    # Check to see if the sources have been estimated
    # in the wrong order
    if abs(np.dot(s1_, s2)) > abs(np.dot(s1_, s1)):
        s2_, s1_ = s_
    s1_ *= np.sign(np.dot(s1_, s1))
    s2_ *= np.sign(np.dot(s2_, s2))

    # Check that we have estimated the original sources
    if add_noise == False:
        assert_almost_equal(np.dot(s1_, s1) / n_samples, 1, decimal=3)
        assert_almost_equal(np.dot(s2_, s2) / n_samples, 1, decimal=3)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
