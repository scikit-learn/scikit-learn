from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils.metaestimators import available_if


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
    assert "This is a mock delegated function" in str(
        MockMetaEstimator.__dict__["func"].__doc__
    )
    assert "This is a mock delegated function" in str(MockMetaEstimator.func.__doc__)
    assert "This is a mock delegated function" in str(MockMetaEstimator().func.__doc__)


class MetaEst:
    """A mock meta estimator"""

    def __init__(self, sub_est, better_sub_est=None):
        self.sub_est = sub_est
        self.better_sub_est = better_sub_est

    @if_delegate_has_method(delegate="sub_est")
    def predict(self):
        pass


class MetaEstTestTuple(MetaEst):
    """A mock meta estimator to test passing a tuple of delegates"""

    @if_delegate_has_method(delegate=("sub_est", "better_sub_est"))
    def predict(self):
        pass


class MetaEstTestList(MetaEst):
    """A mock meta estimator to test passing a list of delegates"""

    @if_delegate_has_method(delegate=["sub_est", "better_sub_est"])
    def predict(self):
        pass


class HasPredict:
    """A mock sub-estimator with predict method"""

    def predict(self):
        pass


class HasNoPredict:
    """A mock sub-estimator with no predict method"""

    pass


def test_if_delegate_has_method():
    assert hasattr(MetaEst(HasPredict()), "predict")
    assert not hasattr(MetaEst(HasNoPredict()), "predict")
    assert not hasattr(MetaEstTestTuple(HasNoPredict(), HasNoPredict()), "predict")
    assert hasattr(MetaEstTestTuple(HasPredict(), HasNoPredict()), "predict")
    assert not hasattr(MetaEstTestTuple(HasNoPredict(), HasPredict()), "predict")
    assert not hasattr(MetaEstTestList(HasNoPredict(), HasPredict()), "predict")
    assert hasattr(MetaEstTestList(HasPredict(), HasPredict()), "predict")


class AvailableParameterEstimator:
    """This estimator's `available` parameter toggles the presence of a method"""

    def __init__(self, available=True):
        self.available = available

    @available_if(lambda est: est.available)
    def available_func(self):
        """This is a mock available_if function"""
        pass


def test_available_if_docstring():
    assert "This is a mock available_if function" in str(
        AvailableParameterEstimator.__dict__["available_func"].__doc__
    )
    assert "This is a mock available_if function" in str(
        AvailableParameterEstimator.available_func.__doc__
    )
    assert "This is a mock available_if function" in str(
        AvailableParameterEstimator().available_func.__doc__
    )


def test_available_if():
    assert hasattr(AvailableParameterEstimator(), "available_func")
    assert not hasattr(AvailableParameterEstimator(available=False), "available_func")
