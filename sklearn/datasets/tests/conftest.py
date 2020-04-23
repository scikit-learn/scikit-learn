""" Network tests are only run, if data is already locally available,
or if download is specifically requested by environment variable."""
import builtins
from os import environ
import pytest
from sklearn.datasets import fetch_20newsgroups
from sklearn.datasets import fetch_20newsgroups_vectorized
from sklearn.datasets import fetch_california_housing
from sklearn.datasets import fetch_covtype
from sklearn.datasets import fetch_kddcup99
from sklearn.datasets import fetch_olivetti_faces
from sklearn.datasets import fetch_rcv1


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
def fetch_20newsgroups_fxt():
    return _wrapped_fetch(fetch_20newsgroups, dataset_name='20newsgroups')


@pytest.fixture
def fetch_20newsgroups_vectorized_fxt():
    return _wrapped_fetch(fetch_20newsgroups_vectorized,
                          dataset_name='20newsgroups_vectorized')


@pytest.fixture
def fetch_california_housing_fxt():
    return _wrapped_fetch(fetch_california_housing,
                          dataset_name='california_housing')


@pytest.fixture
def fetch_covtype_fxt():
    return _wrapped_fetch(fetch_covtype, dataset_name='covtype')


@pytest.fixture
def fetch_kddcup99_fxt():
    return _wrapped_fetch(fetch_kddcup99, dataset_name='kddcup99')


@pytest.fixture
def fetch_olivetti_faces_fxt():
    return _wrapped_fetch(fetch_olivetti_faces, dataset_name='olivetti_faces')


@pytest.fixture
def fetch_rcv1_fxt():
    return _wrapped_fetch(fetch_rcv1, dataset_name='rcv1')


@pytest.fixture
def hide_available_pandas(monkeypatch):
    """ Pretend pandas was not installed. """
    import_orig = builtins.__import__

    def mocked_import(name, *args, **kwargs):
        if name == 'pandas':
            raise ImportError()
        return import_orig(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', mocked_import)
