# Authors: Raghav RV <rvraghav93@gmail.com>
# License: BSD 3 clause


import sys

from sklearn.utils.deprecation import _is_deprecated
from sklearn.utils.deprecation import deprecated
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import SkipTest


@deprecated('qwerty')
class MockClass1:
    pass


class MockClass2:
    @deprecated('mockclass2_method')
    def method(self):
        pass


class MockClass3:
    @deprecated()
    def __init__(self):
        pass


class MockClass4:
    pass


@deprecated()
def mock_function():
    return 10


def test_deprecated():
    assert_warns_message(DeprecationWarning, 'qwerty', MockClass1)
    assert_warns_message(DeprecationWarning, 'mockclass2_method',
                         MockClass2().method)
    assert_warns_message(DeprecationWarning, 'deprecated', MockClass3)
    val = assert_warns_message(DeprecationWarning, 'deprecated', mock_function)
    assert val == 10


def test_is_deprecated():
    if sys.version_info < (3, 5):
        raise SkipTest("This test will run only on python3.5 and above")
    # Test if _is_deprecated helper identifies wrapping via deprecated
    # NOTE it works only for class methods and functions
    assert _is_deprecated(MockClass1.__init__)
    assert _is_deprecated(MockClass2().method)
    assert _is_deprecated(MockClass3.__init__)
    assert not _is_deprecated(MockClass4.__init__)
    assert _is_deprecated(mock_function)
