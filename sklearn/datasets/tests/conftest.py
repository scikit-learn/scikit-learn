""" Network tests are only run, if data is already locally available,
or if download is specifically requested by environment variable."""
from os import environ
import pytest
from sklearn.datasets import (
    fetch_20newsgroups as _fetch_20newsgroups,
    fetch_20newsgroups_vectorized as _fetch_20newsgroups_vectorized,
    fetch_california_housing as _fetch_california_housing,
    fetch_covtype as _fetch_covtype,
    fetch_kddcup99 as _fetch_kddcup99,
    fetch_olivetti_faces as _fetch_olivetti_faces,
    fetch_rcv1 as _fetch_rcv1,
)


def _wrapped_fetch(f, dataset_name):
    """ Fetch dataset (download if missing and requested by environment) """
    download_if_missing = environ.get('SKLEARN_SKIP_NETWORK_TESTS', '1') == '0'

    def wrapped(*args, **kwargs):
        kwargs['download_if_missing'] = download_if_missing
        try:
            return f(*args, **kwargs)
        except IOError:
            pytest.skip("Download {} to run this test".format(dataset_name))

    return wrapped


@pytest.fixture
def fetch_20newsgroups():
    return _wrapped_fetch(_fetch_20newsgroups, dataset_name='20newsgroups')


@pytest.fixture
def fetch_20newsgroups_vectorized():
    return _wrapped_fetch(_fetch_20newsgroups_vectorized,
                          dataset_name='20newsgroups_vectorized')


@pytest.fixture
def fetch_california_housing():
    return _wrapped_fetch(_fetch_california_housing,
                          dataset_name='california_housing')


@pytest.fixture
def fetch_covtype():
    return _wrapped_fetch(_fetch_covtype, dataset_name='covtype')


@pytest.fixture
def fetch_kddcup99():
    return _wrapped_fetch(_fetch_kddcup99, dataset_name='kddcup99')


@pytest.fixture
def fetch_olivetti_faces():
    return _wrapped_fetch(_fetch_olivetti_faces, dataset_name='olivetti_faces')


@pytest.fixture
def fetch_rcv1():
    return _wrapped_fetch(_fetch_rcv1, dataset_name='rcv1')
