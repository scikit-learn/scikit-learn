"""Testing for GaussianMixtureIC"""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from sklearn.exceptions import NotFittedError
from sklearn.metrics import adjusted_rand_score

from sklearn.mixture import GaussianMixtureIC


def _test_inputs(X, error_type, **kws):
    with pytest.raises(error_type):
        gmIC = GaussianMixtureIC(**kws)
        gmIC.fit(X)


def test_n_components():
    # Generate random data
    X = np.random.normal(0, 1, size=(100, 3))

    # min_components < 1
    _test_inputs(X, ValueError, min_components=0)

    # min_components integer
    _test_inputs(X, TypeError, min_components="1")

    # max_components < min_components
    _test_inputs(X, ValueError, max_components=0)

    # max_components integer
    _test_inputs(X, TypeError, max_components="1")

    # max_components > n_samples
    _test_inputs(X, ValueError, max_components=101)

    # min_components > n_samples
    _test_inputs(
        X, ValueError, **{"min_components": 101, "max_components": 102}
    )


def test_input_param():
    # Generate random data
    X = np.random.normal(0, 1, size=(100, 3))

    # affinity is not an array, string or list
    _test_inputs(X, TypeError, affinity=1)

    # affinity is not in ['euclidean', 'manhattan', 'cosine', 'none']
    _test_inputs(X, ValueError, affinity="1")

    # linkage is not an array, string or list
    _test_inputs(X, TypeError, linkage=1)

    # linkage is not in ['single', 'average', 'complete', 'ward']
    _test_inputs(X, ValueError, linkage="1")

    # covariance type is not an array, string or list
    _test_inputs(X, TypeError, covariance_type=1)

    # covariance type is not in ['spherical', 'diag', 'tied', 'full']
    _test_inputs(X, ValueError, covariance_type="1")

    # euclidean is not an affinity option when ward is a linkage option
    _test_inputs(X, ValueError, **{"affinity": "manhattan", "linkage": "ward"})

    # criterion is not "aic" or "bic"
    _test_inputs(X, ValueError, criterion="cic")


def test_labels_init():
    # Generate random data
    X = np.random.normal(0, 1, size=(100, 3))

    # label_init is not a 1-D array
    _test_inputs(X, TypeError, label_init=np.zeros([100, 2]))

    # label_init is not 1-D array, a list or None.
    _test_inputs(X, TypeError, label_init="label")

    # label_init length is not equal to n_samples
    _test_inputs(X, ValueError, label_init=np.zeros([50, 1]))

    # label_init length does not match min_components and max_components
    _test_inputs(
        X,
        ValueError,
        **{
            "label_init": np.zeros([100, 1]),
            "min_components": 2,
            "max_components": 3,
        },
    )


def test_predict_without_fit():
    # Generate random data
    X = np.random.normal(0, 1, size=(100, 3))

    with pytest.raises(NotFittedError):
        gmIC = GaussianMixtureIC(min_components=2)
        gmIC.predict(X)


def test_cosine_with_0():
    # Generate data that contains a zero vector
    X = np.array(
        [
            [0, 1, 0],
            [1, 0, 1],
            [0, 0, 0],
            [1, 1, 0],
            [0, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
            [1, 0, 0],
            [0, 1, 1],
            [1, 1, 0],
            [0, 1, 0],
        ]
    )

    # warning when there are valid affinity options other than "cosine"
    with pytest.warns(UserWarning):
        gmIC = GaussianMixtureIC(min_components=2, affinity="all")
        gmIC.fit(X)

    # error when "cosine" is the only affinity option
    with pytest.raises(ValueError):
        gmIC = GaussianMixtureIC(min_components=2, affinity="cosine")
        gmIC.fit(X)


def test_two_class():
    """
    Easily separable two gaussian problem.
    """
    np.random.seed(1)

    n = 100
    d = 3

    X1 = np.random.normal(2, 0.5, size=(n, d))
    X2 = np.random.normal(-2, 0.5, size=(n, d))
    X = np.vstack((X1, X2))
    y = np.repeat([0, 1], n)

    # test BIC
    gmIC = GaussianMixtureIC(max_components=5, criterion="bic", n_jobs=2)
    gmIC.fit(X, y)
    n_components = gmIC.n_components_

    # Assert that the two cluster model is the best
    assert_equal(n_components, 2)

    # Assert that we get perfect clustering
    ari = adjusted_rand_score(y, gmIC.fit_predict(X))
    assert_allclose(ari, 1)

    # test AIC
    gmIC = GaussianMixtureIC(max_components=5, criterion="aic", n_jobs=2)
    gmIC.fit(X, y)
    n_components = gmIC.n_components_

    # AIC gets the number of components wrong
    assert_equal(n_components >= 1, True)
    assert_equal(n_components <= 5, True)


def test_five_class():
    """
    Easily separable five gaussian problem.
    """
    np.random.seed(1)

    n = 100
    mus = [[i * 5, 0] for i in range(5)]
    cov = np.eye(2)  # balls

    X = np.vstack([np.random.multivariate_normal(mu, cov, n) for mu in mus])

    # test BIC
    gmIC = GaussianMixtureIC(
        min_components=3,
        max_components=10,
        covariance_type="all",
        criterion="bic",
    )
    gmIC.fit(X)
    assert_equal(gmIC.n_components_, 5)

    # test AIC
    gmIC = GaussianMixtureIC(
        min_components=3,
        max_components=10,
        covariance_type="all",
        criterion="aic",
    )
    gmIC.fit(X)
    # AIC fails often so there is no assertion here
    assert_equal(gmIC.n_components_ >= 3, True)
    assert_equal(gmIC.n_components_ <= 10, True)


# Generate random data for each of the four covariance types
n = 100
mu1 = [-10, 0]
mu2 = [10, 0]

# Spherical
np.random.seed(1)
cov1 = 2 * np.eye(2)
cov2 = 2 * np.eye(2)
X1 = np.random.multivariate_normal(mu1, cov1, n)
X2 = np.random.multivariate_normal(mu2, cov2, n)
X_spherical = np.concatenate((X1, X2))

# Diagonal
np.random.seed(10)
cov1 = np.diag([1, 1])
cov2 = np.diag([2, 1])
X1 = np.random.multivariate_normal(mu1, cov1, n)
X2 = np.random.multivariate_normal(mu2, cov2, n)
X_diag = np.concatenate((X1, X2))

# Tied
cov1 = np.array([[2, 1], [1, 2]])
cov2 = np.array([[2, 1], [1, 2]])
X1 = np.random.multivariate_normal(mu1, cov1, n)
X2 = np.random.multivariate_normal(mu2, cov2, n)
X_tied = np.concatenate((X1, X2))

# Full
cov1 = np.array([[2, -1], [-1, 2]])
cov2 = np.array([[2, 1], [1, 2]])
X1 = np.random.multivariate_normal(mu1, cov1, n)
X2 = np.random.multivariate_normal(mu2, cov2, n)
X_full = np.concatenate((X1, X2))


@pytest.mark.parametrize(
    "test_data,expected_cov_type",
    [
        (X_spherical, "spherical"),
        (X_diag, "diag"),
        (X_tied, "tied"),
        (X_full, "full"),
    ],
)
def test_covariances(test_data, expected_cov_type):
    """
    Easily separable two gaussian problem.
    """
    gmIC = GaussianMixtureIC(min_components=2, covariance_type="all")
    gmIC.fit(test_data)
    assert_equal(gmIC.covariance_type_, expected_cov_type)
