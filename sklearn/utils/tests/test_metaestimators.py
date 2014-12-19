from sklearn.utils.metaestimators import if_delegate_has_method
from nose.tools import assert_true


class Prefix(object):
    def func():
        pass


class MockMetaEstimator(object):
    """This is a mock meta estimator"""
    a_prefix = Prefix()

    @if_delegate_has_method(delegate="a_prefix")
    def func(self):
        """This is a mock delegated function"""
        pass


def test_delegated_docstring():
    assert_true("This is a mock delegated function"
                in str(MockMetaEstimator.__dict__['func'].__doc__))
    assert_true("This is a mock delegated function"
                in str(MockMetaEstimator.func.__doc__))
    assert_true("This is a mock delegated function"
                in str(MockMetaEstimator().func.__doc__))
