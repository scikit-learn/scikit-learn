# Author: Gael Varoquaux
# License: BSD 3 clause

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
from sklearn.utils.testing import assert_dict_equal
from sklearn.utils.testing import ignore_warnings

from sklearn.base import BaseEstimator, clone, is_classifier
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn import datasets

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

    est = ModifyInitParams()
    assert_raises(RuntimeError, clone, est)


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


def test_clone_sparse_matrices():
    sparse_matrix_classes = [
        getattr(sp, name)
        for name in dir(sp) if name.endswith('_matrix')]

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


def test_set_params_passes_all_parameters():
    # Make sure all parameters are passed together to set_params
    # of nested estimator. Regression test for #9944

    class TestDecisionTree(DecisionTreeClassifier):
        def set_params(self, **kwargs):
            super(TestDecisionTree, self).set_params(**kwargs)
            # expected_kwargs is in test scope
            assert kwargs == expected_kwargs
            return self

    expected_kwargs = {'max_depth': 5, 'min_samples_leaf': 2}
    for est in [Pipeline([('estimator', TestDecisionTree())]),
                GridSearchCV(TestDecisionTree(), {})]:
        est.set_params(estimator__max_depth=5,
                       estimator__min_samples_leaf=2)


def test_set_params_updates_valid_params():
    # Check that set_params tries to set SVC().C, not
    # DecisionTreeClassifier().C
    gscv = GridSearchCV(DecisionTreeClassifier(), {})
    gscv.set_params(estimator=SVC(), estimator__C=42.0)
    assert gscv.estimator.C == 42.0


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

        def transform(self, X):
            pass

    # build and clone estimator
    d = np.arange(10)
    df = MockDataFrame(d)
    e = DummyEstimator(df, scalar_param=1)
    cloned_e = clone(e)

    # the test
    assert_true((e.df == cloned_e.df).values.all())
    assert_equal(e.scalar_param, cloned_e.scalar_param)


def test_pickle_version_warning_is_not_raised_with_matching_version():
    iris = datasets.load_iris()
    tree = DecisionTreeClassifier().fit(iris.data, iris.target)
    tree_pickle = pickle.dumps(tree)
    assert_true(b"version" in tree_pickle)
    tree_restored = assert_no_warnings(pickle.loads, tree_pickle)

    # test that we can predict with the restored decision tree classifier
    score_of_original = tree.score(iris.data, iris.target)
    score_of_restored = tree_restored.score(iris.data, iris.target)
    assert_equal(score_of_original, score_of_restored)


class TreeBadVersion(DecisionTreeClassifier):
    def __getstate__(self):
        return dict(self.__dict__.items(), _sklearn_version="something")


pickle_error_message = (
    "Trying to unpickle estimator {estimator} from "
    "version {old_version} when using version "
    "{current_version}. This might "
    "lead to breaking code or invalid results. "
    "Use at your own risk.")


def test_pickle_version_warning_is_issued_upon_different_version():
    iris = datasets.load_iris()
    tree = TreeBadVersion().fit(iris.data, iris.target)
    tree_pickle_other = pickle.dumps(tree)
    message = pickle_error_message.format(estimator="TreeBadVersion",
                                          old_version="something",
                                          current_version=sklearn.__version__)
    assert_warns_message(UserWarning, message, pickle.loads, tree_pickle_other)


class TreeNoVersion(DecisionTreeClassifier):
    def __getstate__(self):
        return self.__dict__


def test_pickle_version_warning_is_issued_when_no_version_info_in_pickle():
    iris = datasets.load_iris()
    # TreeNoVersion has no getstate, like pre-0.18
    tree = TreeNoVersion().fit(iris.data, iris.target)

    tree_pickle_noversion = pickle.dumps(tree)
    assert_false(b"version" in tree_pickle_noversion)
    message = pickle_error_message.format(estimator="TreeNoVersion",
                                          old_version="pre-0.18",
                                          current_version=sklearn.__version__)
    # check we got the warning about using pre-0.18 pickle
    assert_warns_message(UserWarning, message, pickle.loads,
                         tree_pickle_noversion)


def test_pickle_version_no_warning_is_issued_with_non_sklearn_estimator():
    iris = datasets.load_iris()
    tree = TreeNoVersion().fit(iris.data, iris.target)
    tree_pickle_noversion = pickle.dumps(tree)
    try:
        module_backup = TreeNoVersion.__module__
        TreeNoVersion.__module__ = "notsklearn"
        assert_no_warnings(pickle.loads, tree_pickle_noversion)
    finally:
        TreeNoVersion.__module__ = module_backup


class DontPickleAttributeMixin(object):
    def __getstate__(self):
        data = self.__dict__.copy()
        data["_attribute_not_pickled"] = None
        return data

    def __setstate__(self, state):
        state["_restored"] = True
        self.__dict__.update(state)


class MultiInheritanceEstimator(BaseEstimator, DontPickleAttributeMixin):
    def __init__(self, attribute_pickled=5):
        self.attribute_pickled = attribute_pickled
        self._attribute_not_pickled = None


def test_pickling_when_getstate_is_overwritten_by_mixin():
    estimator = MultiInheritanceEstimator()
    estimator._attribute_not_pickled = "this attribute should not be pickled"

    serialized = pickle.dumps(estimator)
    estimator_restored = pickle.loads(serialized)
    assert_equal(estimator_restored.attribute_pickled, 5)
    assert_equal(estimator_restored._attribute_not_pickled, None)
    assert_true(estimator_restored._restored)


def test_pickling_when_getstate_is_overwritten_by_mixin_outside_of_sklearn():
    try:
        estimator = MultiInheritanceEstimator()
        text = "this attribute should not be pickled"
        estimator._attribute_not_pickled = text
        old_mod = type(estimator).__module__
        type(estimator).__module__ = "notsklearn"

        serialized = estimator.__getstate__()
        assert_dict_equal(serialized, {'_attribute_not_pickled': None,
                                       'attribute_pickled': 5})

        serialized['attribute_pickled'] = 4
        estimator.__setstate__(serialized)
        assert_equal(estimator.attribute_pickled, 4)
        assert_true(estimator._restored)
    finally:
        type(estimator).__module__ = old_mod


class SingleInheritanceEstimator(BaseEstimator):
    def __init__(self, attribute_pickled=5):
        self.attribute_pickled = attribute_pickled
        self._attribute_not_pickled = None

    def __getstate__(self):
        data = self.__dict__.copy()
        data["_attribute_not_pickled"] = None
        return data


@ignore_warnings(category=(UserWarning))
def test_pickling_works_when_getstate_is_overwritten_in_the_child_class():
    estimator = SingleInheritanceEstimator()
    estimator._attribute_not_pickled = "this attribute should not be pickled"

    serialized = pickle.dumps(estimator)
    estimator_restored = pickle.loads(serialized)
    assert_equal(estimator_restored.attribute_pickled, 5)
    assert_equal(estimator_restored._attribute_not_pickled, None)
