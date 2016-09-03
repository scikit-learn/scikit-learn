import warnings
import numpy as np
from nose.tools import assert_true, assert_false, assert_raises
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.linear_model import SGDClassifier
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from sklearn.grid_search import GridSearchCV


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


def test_stochastic_gradient_loss_param():
    # Make sure the predict_proba works when loss is specified
    # as one of the parameters in the param_grid.
    param_grid = {
        'loss': ['log'],
    }
    X = np.arange(20).reshape(5, -1)
    y = [0, 0, 1, 1, 1]
    clf = GridSearchCV(estimator=SGDClassifier(loss='hinge'),
                       param_grid=param_grid)

    # When the estimator is not fitted, `predict_proba` is not available as the
    # loss is 'hinge'.
    assert_false(hasattr(clf, "predict_proba"))
    clf.fit(X, y)
    clf.predict_proba(X)
    clf.predict_log_proba(X)

    # Make sure `predict_proba` is not available when setting loss=['hinge']
    # in param_grid
    param_grid = {
        'loss': ['hinge'],
    }
    clf = GridSearchCV(estimator=SGDClassifier(loss='hinge'),
                       param_grid=param_grid)
    assert_false(hasattr(clf, "predict_proba"))
    clf.fit(X, y)
    assert_false(hasattr(clf, "predict_proba"))


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
