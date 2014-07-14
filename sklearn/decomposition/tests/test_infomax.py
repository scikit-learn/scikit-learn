# Authors: Denis A. Engemann <denis.engemann@gmail.com>
#
# License: BSD (3-clause)

"""
Test the infomax algorithm.
Parts of this code are taken from scikit-learn
"""

import itertools
import numpy as np


from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal

from scipy import stats

from sklearn.decomposition import InfomaxICA, infomax, RandomizedPCA


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


def test_infomax_simple(add_noise=False):
    """ Test the infomax algorithm on very simple data.
    """
    rng = np.random.RandomState(0)
    # scipy.stats uses the global RNG:
    np.random.seed(0)
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

    extended = [True, False]
    whiten = [True, False]
    for extended, whiten in itertools.product(extended, whiten):
        if whiten:
            k_, mixing_, s_ = infomax(m.T, extended=extended)
        else:
            X = RandomizedPCA(n_components=2, whiten=True).fit_transform(m.T)
            k_, mixing_, s_ = infomax(X, extended=extended, whiten=False)

        center_and_norm(s_)
        s_ = s_.T
        s1_, s2_ = s_
        # Check to see if the sources have been estimated
        # in the wrong order
        if abs(np.dot(s1_, s2)) > abs(np.dot(s1_, s1)):
            s2_, s1_ = s_
        s1_ *= np.sign(np.dot(s1_, s1))
        s2_ *= np.sign(np.dot(s2_, s2))

        # Check that we have estimated the original sources
        if not add_noise:
            assert_almost_equal(np.dot(s1_, s1) / n_samples, 1, decimal=2)
            assert_almost_equal(np.dot(s2_, s2) / n_samples, 1, decimal=2)
        else:
            assert_almost_equal(np.dot(s1_, s1) / n_samples, 1, decimal=1)
            assert_almost_equal(np.dot(s2_, s2) / n_samples, 1, decimal=1)

        # Test FastICA class
        _, _, sources_fun = infomax(m.T, extneded=extended, whiten=whiten,
                                    random_state=0)
        ica = InfomaxICA(extended=extended, whiten=whiten, random_state=0)
        sources = ica.fit_transform(m.T)
        assert_equal(ica.components_.shape, (2, 2))
        assert_equal(sources.shape, (1000, 2))

        assert_array_almost_equal(sources_fun, sources)
        assert_array_almost_equal(sources, ica.transform(m.T))

        assert_equal(ica.mixing_.shape, (2, 2))


def test_non_square_infomax(add_noise=False):
    """ Test the infomax algorithm on very simple data.
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
    n_observed = 6
    mixing = rng.randn(n_observed, 2)
    m = np.dot(mixing, s)

    if add_noise:
        m += 0.1 * rng.randn(n_observed, n_samples)

    center_and_norm(m)

    k_, mixing_, s_ = infomax(m.T, n_components=2, random_state=rng)
    s = s_.T
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
    if not add_noise:
        assert_almost_equal(np.dot(s1_, s1) / n_samples, 1, decimal=2)
        assert_almost_equal(np.dot(s2_, s2) / n_samples, 1, decimal=2)


def test_fit_transform():
    """Test FastICA.fit_transform"""
    rng = np.random.RandomState(0)
    X = rng.random_sample((100, 10))
    for whiten, n_components in [[True, 5], [False, 10]]:

        ica = InfomaxICA(n_components=n_components, whiten=whiten,
                         random_state=0)
        Xt = ica.fit_transform(X)
        assert_equal(ica.components_.shape, (n_components, 10))
        assert_equal(Xt.shape, (100, n_components))

        ica = InfomaxICA(n_components=n_components, whiten=whiten,
                         random_state=0)
        ica.fit(X)
        assert_equal(ica.components_.shape, (n_components, 10))
        Xt2 = ica.transform(X)

        assert_array_almost_equal(Xt, Xt2)


def test_inverse_transform():
    """Test FastICA.inverse_transform"""
    n_features = 10
    n_samples = 100
    n1, n2 = 5, 10
    rng = np.random.RandomState(0)
    X = rng.random_sample((n_samples, n_features))
    expected = {(True, n1): (n_features, n1),
                (True, n2): (n_features, n2),
                (False, n1): (n_features, n2),
                (False, n2): (n_features, n2)}
    for whiten in [True, False]:
        for n_components in [n1, n2]:
            ica = InfomaxICA(n_components=n_components, random_state=rng,
                             whiten=whiten)
            Xt = ica.fit_transform(X)
            expected_shape = expected[(whiten, n_components)]
            assert_equal(ica.mixing_.shape, expected_shape)
            X2 = ica.inverse_transform(Xt)
            assert_equal(X.shape, X2.shape)

            # reversibility test in non-reduction case
            if n_components == X.shape[1]:
                assert_array_almost_equal(X, X2)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
