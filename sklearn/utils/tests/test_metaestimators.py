from sklearn.utils.testing import assert_true, assert_false
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils.metaestimators import _Router


class Prefix(object):
    def func(self):
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


class MetaEst(object):
    """A mock meta estimator"""
    def __init__(self, sub_est, better_sub_est=None):
        self.sub_est = sub_est
        self.better_sub_est = better_sub_est

    @if_delegate_has_method(delegate='sub_est')
    def predict(self):
        pass


class MetaEstTestTuple(MetaEst):
    """A mock meta estimator to test passing a tuple of delegates"""

    @if_delegate_has_method(delegate=('sub_est', 'better_sub_est'))
    def predict(self):
        pass


class MetaEstTestList(MetaEst):
    """A mock meta estimator to test passing a list of delegates"""

    @if_delegate_has_method(delegate=['sub_est', 'better_sub_est'])
    def predict(self):
        pass


class HasPredict(object):
    """A mock sub-estimator with predict method"""

    def predict(self):
        pass


class HasNoPredict(object):
    """A mock sub-estimator with no predict method"""
    pass


def test_if_delegate_has_method():
    assert_true(hasattr(MetaEst(HasPredict()), 'predict'))
    assert_false(hasattr(MetaEst(HasNoPredict()), 'predict'))
    assert_false(
        hasattr(MetaEstTestTuple(HasNoPredict(), HasNoPredict()), 'predict'))
    assert_true(
        hasattr(MetaEstTestTuple(HasPredict(), HasNoPredict()), 'predict'))
    assert_false(
        hasattr(MetaEstTestTuple(HasNoPredict(), HasPredict()), 'predict'))
    assert_false(
        hasattr(MetaEstTestList(HasNoPredict(), HasPredict()), 'predict'))
    assert_true(
        hasattr(MetaEstTestList(HasPredict(), HasPredict()), 'predict'))


FOO = [1, 2, 3]
BAR = [4, 5, 6]
TEST_ROUTER_BASIC_PARAMS = [
    ({'*': [('*', '*')]},
     {'bar': BAR, 'foo': FOO}, {'bar': BAR, 'foo': FOO}, set()),
    ({'foo': [('*', '*')]},
     {'foo': FOO}, {'foo': FOO}, {'bar'}),
    ({'foo': ('*', '*')},
     {'foo': FOO}, {'foo': FOO}, {'bar'}),
    ({'foo': [('d*1', '*')]},
     {'foo': FOO}, {}, {'bar'}),
    ({'foo': [('d*1', 'goo')]},
     {'goo': FOO}, {}, {'bar'}),
    ({'foo': [('dest1', 'goo')]},
     {'goo': FOO}, {}, {'bar'}),
    ({'foo': [('dest11', 'goo')]},
     {}, {}, {'foo', 'bar'}),
]

# pytest.mark.parametrize(['routing', 'expected_dest1',
#                           'expected_dest2', 'expected_remainder'],
#                         TEST_ROUTER_BASIC_PARAMS)


def test_router_basic():
    for params in TEST_ROUTER_BASIC_PARAMS:
        print(params)
        routing, expected_dest1, expected_dest2, expected_remainder = params
        router = _Router(routing)
        (dest1, dest2), remainder = router({'foo': FOO, 'bar': BAR},
                                           ['dest1', 'dest2'])
        assert dest1 == expected_dest1
        assert dest2 == expected_dest2
        assert remainder == expected_remainder
