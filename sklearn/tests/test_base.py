# Author: Gael Varoquaux
# License: BSD 3 clause

import sys

import numpy as np
import scipy.sparse as sp

import sklearn
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import assert_warns_message

from sklearn.base import BaseEstimator, clone, is_classifier
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn import datasets
from sklearn.utils import deprecated

from sklearn.base import TransformerMixin
from sklearn.utils.mocking import MockDataFrame
import pickle


#############################################################################
# A few test classes
class MyEstimator(BaseEstimator):

    def __init__(self, l1=0, empty=None):
        self.l1 = l1
        self.empty = empty


class K(BaseEstimator):
    def __init__(self, c=None, d=None):
        self.c = c
        self.d = d


class T(BaseEstimator):
    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b


class ModifyInitParams(BaseEstimator):
    """Deprecated behavior.
    Equal parameters but with a type cast.
    Doesn't fulfill a is a
    """
    def __init__(self, a=np.array([0])):
        self.a = a.copy()


class DeprecatedAttributeEstimator(BaseEstimator):
    def __init__(self, a=None, b=None):
        self.a = a
        if b is not None:
            DeprecationWarning("b is deprecated and renamed 'a'")
            self.a = b

    @property
    @deprecated("Parameter 'b' is deprecated and renamed to 'a'")
    def b(self):
        return self._b


class Buggy(BaseEstimator):
    " A buggy estimator that does not set its parameters right. "

    def __init__(self, a=None):
        self.a = 1


class NoEstimator(object):
    def __init__(self):
        pass

    def fit(self, X=None, y=None):
        return self

    def predict(self, X=None):
        return None


class VargEstimator(BaseEstimator):
    """scikit-learn estimators shouldn't have vargs."""
    def __init__(self, *vargs):
        pass


#############################################################################
# The tests

def test_clone():
    # Tests that clone creates a correct deep copy.
    # We create an estimator, make a copy of its original state
    # (which, in this case, is the current state of the estimator),
    # and check that the obtained copy is a correct deep copy.

    from sklearn.feature_selection import SelectFpr, f_classif

    selector = SelectFpr(f_classif, alpha=0.1)
    new_selector = clone(selector)
    assert_true(selector is not new_selector)
    assert_equal(selector.get_params(), new_selector.get_params())

    selector = SelectFpr(f_classif, alpha=np.zeros((10, 2)))
    new_selector = clone(selector)
    assert_true(selector is not new_selector)


def test_clone_2():
    # Tests that clone doesn't copy everything.
    # We first create an estimator, give it an own attribute, and
    # make a copy of its original state. Then we check that the copy doesn't
    # have the specific attribute we manually added to the initial estimator.

    from sklearn.feature_selection import SelectFpr, f_classif

    selector = SelectFpr(f_classif, alpha=0.1)
    selector.own_attribute = "test"
    new_selector = clone(selector)
    assert_false(hasattr(new_selector, "own_attribute"))


def test_clone_buggy():
    # Check that clone raises an error on buggy estimators.
    buggy = Buggy()
    buggy.a = 2
    assert_raises(RuntimeError, clone, buggy)

    no_estimator = NoEstimator()
    assert_raises(TypeError, clone, no_estimator)

    varg_est = VargEstimator()
    assert_raises(RuntimeError, clone, varg_est)


def test_clone_empty_array():
    # Regression test for cloning estimators with empty arrays
    clf = MyEstimator(empty=np.array([]))
    clf2 = clone(clf)
    assert_array_equal(clf.empty, clf2.empty)

    clf = MyEstimator(empty=sp.csr_matrix(np.array([[0]])))
    clf2 = clone(clf)
    assert_array_equal(clf.empty.data, clf2.empty.data)


def test_clone_nan():
    # Regression test for cloning estimators with default parameter as np.nan
    clf = MyEstimator(empty=np.nan)
    clf2 = clone(clf)

    assert_true(clf.empty is clf2.empty)


def test_clone_copy_init_params():
    # test for deprecation warning when copying or casting an init parameter
    est = ModifyInitParams()
    message = ("Estimator ModifyInitParams modifies parameters in __init__. "
               "This behavior is deprecated as of 0.18 and support "
               "for this behavior will be removed in 0.20.")

    assert_warns_message(DeprecationWarning, message, clone, est)


def test_clone_sparse_matrices():
    sparse_matrix_classes = [
        getattr(sp, name)
        for name in dir(sp) if name.endswith('_matrix')]

    PY26 = sys.version_info[:2] == (2, 6)
    if PY26:
        # sp.dok_matrix can not be deepcopied in Python 2.6
        sparse_matrix_classes.remove(sp.dok_matrix)

    for cls in sparse_matrix_classes:
        sparse_matrix = cls(np.eye(5))
        clf = MyEstimator(empty=sparse_matrix)
        clf_cloned = clone(clf)
        assert_true(clf.empty.__class__ is clf_cloned.empty.__class__)
        assert_array_equal(clf.empty.toarray(), clf_cloned.empty.toarray())


def test_repr():
    # Smoke test the repr of the base estimator.
    my_estimator = MyEstimator()
    repr(my_estimator)
    test = T(K(), K())
    assert_equal(
        repr(test),
        "T(a=K(c=None, d=None), b=K(c=None, d=None))"
    )

    some_est = T(a=["long_params"] * 1000)
    assert_equal(len(repr(some_est)), 415)


def test_str():
    # Smoke test the str of the base estimator
    my_estimator = MyEstimator()
    str(my_estimator)


def test_get_params():
    test = T(K(), K())

    assert_true('a__d' in test.get_params(deep=True))
    assert_true('a__d' not in test.get_params(deep=False))

    test.set_params(a__d=2)
    assert_true(test.a.d == 2)
    assert_raises(ValueError, test.set_params, a__a=2)


def test_get_params_deprecated():
    # deprecated attribute should not show up as params
    est = DeprecatedAttributeEstimator(a=1)

    assert_true('a' in est.get_params())
    assert_true('a' in est.get_params(deep=True))
    assert_true('a' in est.get_params(deep=False))

    assert_true('b' not in est.get_params())
    assert_true('b' not in est.get_params(deep=True))
    assert_true('b' not in est.get_params(deep=False))


def test_is_classifier():
    svc = SVC()
    assert_true(is_classifier(svc))
    assert_true(is_classifier(GridSearchCV(svc, {'C': [0.1, 1]})))
    assert_true(is_classifier(Pipeline([('svc', svc)])))
    assert_true(is_classifier(Pipeline(
        [('svc_cv', GridSearchCV(svc, {'C': [0.1, 1]}))])))


def test_set_params():
    # test nested estimator parameter setting
    clf = Pipeline([("svc", SVC())])
    # non-existing parameter in svc
    assert_raises(ValueError, clf.set_params, svc__stupid_param=True)
    # non-existing parameter of pipeline
    assert_raises(ValueError, clf.set_params, svm__stupid_param=True)
    # we don't currently catch if the things in pipeline are estimators
    # bad_pipeline = Pipeline([("bad", NoEstimator())])
    # assert_raises(AttributeError, bad_pipeline.set_params,
    #               bad__stupid_param=True)


def test_score_sample_weight():

    rng = np.random.RandomState(0)

    # test both ClassifierMixin and RegressorMixin
    estimators = [DecisionTreeClassifier(max_depth=2),
                  DecisionTreeRegressor(max_depth=2)]
    sets = [datasets.load_iris(),
            datasets.load_boston()]

    for est, ds in zip(estimators, sets):
        est.fit(ds.data, ds.target)
        # generate random sample weights
        sample_weight = rng.randint(1, 10, size=len(ds.target))
        # check that the score with and without sample weights are different
        assert_not_equal(est.score(ds.data, ds.target),
                         est.score(ds.data, ds.target,
                                   sample_weight=sample_weight),
                         msg="Unweighted and weighted scores "
                             "are unexpectedly equal")


def test_clone_pandas_dataframe():

    class DummyEstimator(BaseEstimator, TransformerMixin):
        """This is a dummy class for generating numerical features

        This feature extractor extracts numerical features from pandas data
        frame.

        Parameters
        ----------

        df: pandas data frame
            The pandas data frame parameter.

        Notes
        -----
        """
        def __init__(self, df=None, scalar_param=1):
            self.df = df
            self.scalar_param = scalar_param

        def fit(self, X, y=None):
            pass

        def transform(self, X, y=None):
            pass

    # build and clone estimator
    d = np.arange(10)
    df = MockDataFrame(d)
    e = DummyEstimator(df, scalar_param=1)
    cloned_e = clone(e)

    # the test
    assert_true((e.df == cloned_e.df).values.all())
    assert_equal(e.scalar_param, cloned_e.scalar_param)


class TreeNoVersion(DecisionTreeClassifier):
    def __getstate__(self):
        return self.__dict__


class TreeBadVersion(DecisionTreeClassifier):
    def __getstate__(self):
        return dict(self.__dict__.items(), _sklearn_version="something")


def test_pickle_version_warning():
    # check that warnings are raised when unpickling in a different version

    # first, check no warning when in the same version:
    iris = datasets.load_iris()
    tree = DecisionTreeClassifier().fit(iris.data, iris.target)
    tree_pickle = pickle.dumps(tree)
    assert_true(b"version" in tree_pickle)
    assert_no_warnings(pickle.loads, tree_pickle)

    # check that warning is raised on different version
    tree = TreeBadVersion().fit(iris.data, iris.target)
    tree_pickle_other = pickle.dumps(tree)
    message = ("Trying to unpickle estimator TreeBadVersion from "
               "version {0} when using version {1}. This might lead to "
               "breaking code or invalid results. "
               "Use at your own risk.".format("something",
                                              sklearn.__version__))
    assert_warns_message(UserWarning, message, pickle.loads, tree_pickle_other)

    # check that not including any version also works:
    # TreeNoVersion has no getstate, like pre-0.18
    tree = TreeNoVersion().fit(iris.data, iris.target)

    tree_pickle_noversion = pickle.dumps(tree)
    assert_false(b"version" in tree_pickle_noversion)
    message = message.replace("something", "pre-0.18")
    message = message.replace("TreeBadVersion", "TreeNoVersion")
    # check we got the warning about using pre-0.18 pickle
    assert_warns_message(UserWarning, message, pickle.loads,
                         tree_pickle_noversion)

    # check that no warning is raised for external estimators
    TreeNoVersion.__module__ = "notsklearn"
    assert_no_warnings(pickle.loads, tree_pickle_noversion)
