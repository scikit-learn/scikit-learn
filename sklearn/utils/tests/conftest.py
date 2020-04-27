import pytest

import sklearn


@pytest.fixture(scope='function')
def with_changed_only_true():
    sklearn.set_config(print_changed_only=True)
    yield
    sklearn.set_config(print_changed_only=None)  # reset to default


@pytest.fixture(scope='function')
def with_changed_only_false():
    sklearn.set_config(print_changed_only=False)
    yield
    sklearn.set_config(print_changed_only=None)  # reset to default
