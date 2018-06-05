import warnings

import pytest

from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_warns_message


def warns():
    warnings.warn("some warning")
    return 1


@ignore_warnings()
def test_1():
    print('before:', warnings.filters)
    assert_warns_message(UserWarning, 'some warning', warns)
    print('after:', warnings.filters)
    # This warning is visible because assert_warns_message resets
    # warnings.filters.
    warns()


def test_12():
    print('test12:', warnings.filters)


ignore_common_warnings = pytest.mark.filterwarnings('ignore::UserWarning')


@ignore_common_warnings
def test_2():
    warns()


@ignore_common_warnings
def test_3():
    assert_warns_message(UserWarning, 'some warning', warns)
    # This warning is visible
    warns()
