import pytest

import sklearn


@pytest.fixture(scope='function')
def print_changed_only_setter(value):
    sklearn.set_config(print_changed_only=value)
    yield
    sklearn.set_config(print_changed_only=None)  # reset to default
