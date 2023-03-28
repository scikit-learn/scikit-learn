import pytest
from functools import partial
import numpy as np

from sklearn.mixture._base import (
    _check_responsibilities,
    get_responsibilities,
)
from sklearn import cluster
from .test_gaussian_mixture import RandomData


@pytest.fixture(scope="module")
def random_data():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)

    n_components = rand_data.n_components
    X = rand_data.X["full"]
    n_samples = X.shape[0]

    return rng, n_components, X, n_samples


def kmeans_resp(random_data):
    rng, n_components, X, n_samples = random_data
    # old resp
    resp = np.zeros((n_samples, n_components))
    label = (
        cluster.KMeans(n_clusters=n_components, n_init=1, random_state=rng)
        .fit(X)
        .labels_
    )
    resp[np.arange(n_samples), label] = 1
    # new resp
    labels = np.argmax(resp, axis=1)
    new_resp = get_responsibilities(n_samples, n_components, labels=labels)
    return resp, new_resp


def random_from_data_resp(random_data):
    rng, n_components, X, n_samples = random_data

    # old resp
    resp = np.zeros((n_samples, n_components))
    indices = rng.choice(n_samples, size=n_components, replace=False)
    resp[indices, np.arange(n_components)] = 1

    # new resp
    new_resp = get_responsibilities(n_samples, n_components, indices=indices)
    return resp, new_resp


def kmeans_plusplus_resp(random_data):
    rng, n_components, X, n_samples = random_data

    # old resp
    resp = np.zeros((n_samples, n_components))
    _, indices = cluster.kmeans_plusplus(X, n_components, random_state=rng)
    resp[indices, np.arange(n_components)] = 1

    # new resp
    new_resp = get_responsibilities(n_samples, n_components, indices=indices)
    return resp, new_resp


INIT_METHODS = {
    "kmeans": kmeans_resp,
    "random_from_data": random_from_data_resp,
    "k-means++": kmeans_plusplus_resp,
}


@pytest.mark.parametrize("init_method", ["kmeans", "random_from_data", "k-means++"])
def test_get_responsibilities_regression(init_method, random_data):
    """
    Test that get_responsibilities returns the same result as the old method
    for existing non-callable init_params."""

    resp, new_resp = INIT_METHODS[init_method](random_data)
    assert resp == pytest.approx(new_resp)


def test_check_responsibilities(random_data):
    rng, n_components, X, n_samples = random_data
    # Use default kmeans
    resp, new_resp = kmeans_resp(random_data)

    # Check roundtrip
    resp = _check_responsibilities(resp, n_components=n_components, n_samples=n_samples)
    new_resp = _check_responsibilities(
        new_resp, n_components=n_components, n_samples=n_samples
    )
    assert resp == pytest.approx(new_resp)

    # check bad shape
    with pytest.raises(ValueError):
        _check_responsibilities(
            resp[3:], n_components=n_components, n_samples=n_samples
        )

    # check bad range
    with pytest.raises(ValueError):
        _check_responsibilities(
            resp + 2, n_components=n_components, n_samples=n_samples
        )

    # check negative
    with pytest.raises(ValueError):
        _check_responsibilities(
            resp - 2, n_components=n_components, n_samples=n_samples
        )


def test_random_init(random_data):
    """
    Test that random initialization works as expected with check
    responsibilities. This is separate because random is only default
    that uses non-binary responsibilities."""

    rng, n_components, X, n_samples = random_data

    # Check output is proper responsibilities
    resp = rng.uniform(size=(n_samples, n_components))
    resp /= resp.sum(axis=1)[:, np.newaxis]
    resp = _check_responsibilities(resp, n_components, n_samples)
    assert resp.shape == (n_samples, n_components)


def test_callable(random_data):
    """
    Test that mocked callable methods inside BaseMixture work as expected."""

    rng, n_components, X, n_samples = random_data

    def callable_init(X, n_components, n_samples, rng):
        """Callable init for GaussianMixture"""
        kmean_ = cluster.KMeans(n_clusters=n_components, n_init=1, random_state=rng)
        kmean_.fit(X)
        labels = kmean_.labels_

        return get_responsibilities(
            n_samples=n_samples, n_components=n_components, labels=labels
        )

    test_callable = partial(
        callable_init, n_components=n_components, n_samples=X.shape[0], rng=rng
    )

    # Check output is proper responsibilities
    resp = test_callable(X)
    resp = _check_responsibilities(resp, n_components, n_samples)
    assert resp.shape == (n_samples, n_components)

    # Check bad callable
    def bad_callable(a, b, c):
        return np.array([1, 1, 1])

    bad_resp = bad_callable(1, 2, 3)
    with pytest.raises(ValueError):
        _check_responsibilities(bad_resp, n_components, n_samples)
