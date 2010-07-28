from nose.tools import assert_raises

from scikits.learn.base import BaseEstimator


class MyEstimator(BaseEstimator):

    def __init__(self, l1=0):
        self.l1 = l1

def test_repr():
    """ Smoke test the repr of the 
    """
    my_estimator = MyEstimator()
    repr(my_estimator)

