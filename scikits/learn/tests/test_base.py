from nose.tools import assert_raises

from scikits.learn.base import BaseEstimator, BaseClassifier, BaseRegressor


class MyEstimator(BaseEstimator):

    def __init__(self, l1=0):
        self.l1 = l1


class MyRegressor(BaseRegressor):

    def __init__(self, l1=0):
        self.l1 = l1


class MyClassifier(BaseClassifier):

    def __init__(self, l1=0):
        self.l1 = l1


def test_repr():
    """ Smoke test the repr of the 
    """
    my_estimator = MyEstimator()
    repr(my_estimator)

    my_regressor = MyRegressor()
    repr(my_regressor)

    my_classifier = MyClassifier()
    repr(my_classifier)

