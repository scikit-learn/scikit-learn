from nose.tools import assert_raises

from scikits.learn.base_estimator import BaseEstimator

class IncorrectEstimator(BaseEstimator):
    # A class to check that it is impossible to instanciate an estimator
    # without specifying its parameters
    pass


class MyEstimator(BaseEstimator):
    _params = {'l1': float}

    def __init__(self, l1=0):
        self._set_params(l1=l1)


def test_incorrect_estimator():
    assert_raises(AssertionError, IncorrectEstimator)


def test_repr():
    """ Smoke test the repr of the 
    """
    my_estimator = MyEstimator()
    repr(my_estimator)

