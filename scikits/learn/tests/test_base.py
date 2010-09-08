from nose.tools import assert_true, assert_false, assert_equal
from ..base import BaseEstimator

class MyEstimator(BaseEstimator):

    def __init__(self, l1=0):
        self.l1 = l1

def test_reinit():
    """Tests that BaseEstimator._new() creates a correct deep copy.

    We create an estimator, make a copy of its original state
    (which, in this case, is the current state of the setimator),
    and check that the obtained copy is a correct deep copy.

    """
    from scikits.learn.feature_selection import SelectFpr, f_classif

    selector = SelectFpr(f_classif, alpha=0.1)
    new_selector = selector._reinit()
    assert_true(selector is not new_selector)
    assert_equal(selector._get_params(), new_selector._get_params())
    
def test_reinit_2():
    """Tests that BaseEstimator._new() doesn't copy everything.

    We first create an estimator, give it an own attribute, and
    make a copy of its original state. Then we check that the copy doesn't have
    the specific attribute we manually added to the initial estimator.
    
    """
    from scikits.learn.feature_selection import SelectFpr, f_classif

    selector = SelectFpr(f_classif, alpha=0.1)
    selector.own_attribute = "test"
    new_selector = selector._reinit()
    assert_false(hasattr(new_selector, "own_attribute"))
    

def test_repr():
    """ Smoke test the repr of the 
    """
    my_estimator = MyEstimator()
    repr(my_estimator)

