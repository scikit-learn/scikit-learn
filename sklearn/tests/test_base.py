
# Author: Gael Varoquaux
# License: BSD

from nose.tools import assert_true, assert_false, assert_equal, \
    assert_raises
from ..base import BaseEstimator, clone, is_classifier


#############################################################################
# A few test classes
class MyEstimator(BaseEstimator):

    def __init__(self, l1=0):
        self.l1 = l1


class K(BaseEstimator):
    def __init__(self, c=None, d=None):
        self.c = c
        self.d = d


class T(BaseEstimator):
    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b


class Buggy(BaseEstimator):
    " A buggy estimator that does not set its parameters right. "

    def __init__(self, a=None):
        self.a = 1


#############################################################################
# The tests

def test_clone():
    """Tests that clone creates a correct deep copy.

    We create an estimator, make a copy of its original state
    (which, in this case, is the current state of the setimator),
    and check that the obtained copy is a correct deep copy.

    """
    from sklearn.feature_selection import SelectFpr, f_classif

    selector = SelectFpr(f_classif, alpha=0.1)
    new_selector = clone(selector)
    assert_true(selector is not new_selector)
    assert_equal(selector._get_params(), new_selector._get_params())


def test_clone_2():
    """Tests that clone doesn't copy everything.

    We first create an estimator, give it an own attribute, and
    make a copy of its original state. Then we check that the copy doesn't have
    the specific attribute we manually added to the initial estimator.

    """
    from sklearn.feature_selection import SelectFpr, f_classif

    selector = SelectFpr(f_classif, alpha=0.1)
    selector.own_attribute = "test"
    new_selector = clone(selector)
    assert_false(hasattr(new_selector, "own_attribute"))


def test_clone_buggy():
    """ Check that clone raises an error on buggy estimators """
    buggy = Buggy()
    buggy.a = 2
    assert_raises(AssertionError, clone, buggy)


def test_repr():
    """ Smoke test the repr of the
    """
    my_estimator = MyEstimator()
    repr(my_estimator)
    test = T(K(), K())
    assert_equal(repr(test),
                "T(a=K(c=None, d=None), b=K(c=None, d=None))"
                )


def test_str():
    """ Smoke test the str of the
    """
    my_estimator = MyEstimator()
    str(my_estimator)


def test_get_params():
    test = T(K(), K())

    assert_true('a__d' in test._get_params(deep=True))
    assert_true('a__d' not in test._get_params(deep=False))

    test.set_params(a__d=2)
    assert test.a.d == 2
    assert_raises(AssertionError, test.set_params, a__a=2)


def test_is_classifier():
    from ..svm import SVC
    from ..pipeline import Pipeline
    from ..grid_search import GridSearchCV
    svc = SVC()
    assert_true(is_classifier(svc))
    assert_true(is_classifier(GridSearchCV(svc, {'C': [0.1, 1]})))
    assert_true(is_classifier(Pipeline([('svc', svc)])))
    assert_true(is_classifier(Pipeline([('svc_cv',
                              GridSearchCV(svc, {'C': [0.1, 1]}))])))
