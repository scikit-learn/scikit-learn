import sys
import pytest

from sklearn.utils.metaestimators import if_delegate_has_method


class Prefix:
    def func(self):
        pass


class MockMetaEstimator:
    """This is a mock meta estimator"""
    a_prefix = Prefix()

    @if_delegate_has_method(delegate="a_prefix")
    def func(self):
        """This is a mock delegated function"""
        pass


def test_delegated_docstring():
    assert "This is a mock delegated function" \
                in str(MockMetaEstimator.__dict__['func'].__doc__)
    assert "This is a mock delegated function" \
           in str(MockMetaEstimator.func.__doc__)
    assert "This is a mock delegated function" \
           in str(MockMetaEstimator().func.__doc__)


class MetaEst:
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


class MetaEstMultipleDelegations(MetaEst):
    """A mock meta estimator to test multiple delegations"""

    @if_delegate_has_method(delegate='sub_est')
    def fit(self):
        pass

    @if_delegate_has_method(delegate='better_sub_est')
    def fit_predict(self):
        pass

    @if_delegate_has_method(delegate='sub_est')
    def does_not_exist(self):
        pass


class HasPredict:
    """A mock sub-estimator with predict method"""

    def predict(self):
        pass


class HasNoPredict:
    """A mock sub-estimator with no predict method"""
    pass


class HasFitPredict(HasPredict):
    """A mock sub-estimator with fit and predict methods"""

    def fit(self):
        pass

    def fit_predict(self):
        pass


def test_if_delegate_has_method():
    assert hasattr(MetaEst(HasPredict()), 'predict')
    assert not hasattr(MetaEst(HasNoPredict()), 'predict')
    assert not hasattr(MetaEstTestTuple(HasNoPredict(), HasNoPredict()),
                       'predict')
    assert hasattr(MetaEstTestTuple(HasPredict(), HasNoPredict()), 'predict')
    assert not hasattr(MetaEstTestTuple(HasNoPredict(), HasPredict()),
                       'predict')
    assert not hasattr(MetaEstTestList(HasNoPredict(), HasPredict()),
                       'predict')
    assert hasattr(MetaEstTestList(HasPredict(), HasPredict()), 'predict')

    obj = MetaEstMultipleDelegations(HasFitPredict())
    assert hasattr(obj, 'predict')
    assert hasattr(obj, 'fit')
    assert not hasattr(obj, 'does_not_exist')
    assert not hasattr(obj, 'fit_predict')


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires __set_name__")
def test_if_delegate_has_method_dir():
    assert 'predict' in dir(MetaEst(HasPredict()))
    assert 'predict' not in dir(MetaEst(HasNoPredict()))
    assert 'predict' not in dir(MetaEstTestTuple(HasNoPredict(), HasNoPredict())) # noqa
    assert 'predict' in dir(MetaEstTestTuple(HasPredict(), HasNoPredict()))
    assert 'predict' not in dir(MetaEstTestTuple(HasNoPredict(), HasPredict()))
    assert 'predict' not in dir(MetaEstTestList(HasNoPredict(), HasPredict()))
    assert 'predict' in dir(MetaEstTestList(HasPredict(), HasPredict()))

    attrs = dir(MetaEstMultipleDelegations(HasFitPredict()))
    assert 'predict' in attrs
    assert 'fit' in attrs
    assert 'does_not_exist' not in attrs
    assert 'fit_predict' not in attrs
