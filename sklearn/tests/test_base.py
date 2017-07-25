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
from sklearn.utils import deprecated

from sklearn.base import TransformerMixin, frozen_fit
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

    for cls in sparse_matrix_classes:
        sparse_matrix = cls(np.eye(5))
        clf = MyEstimator(empty=sparse_matrix)
        clf_cloned = clone(clf)
        assert_true(clf.empty.__class__ is clf_cloned.empty.__class__)
        assert_array_equal(clf.empty.toarray(), clf_cloned.empty.toarray())


def test_clone_frozen():
    est = DecisionTreeClassifier()
    est.fit([[0], [1]], [0, 1])
    assert clone(est) is not est
    est.frozen = True
    assert clone(est) is est
    # is still fitted
    assert est.predict([[0]]) == 0
    # freezing works recursively
    seq = [est]
    assert clone(seq) is not seq
    assert clone(seq)[0] is est

    # freezing first then fitting works too
    est = DecisionTreeClassifier()
    est.frozen = True
    est.fit([[0], [1]], [0, 1])
    assert clone(est) is est
    assert est.predict([[0]]) == 0

    # can freeze an unfitted estimator (not sure why, but worth testing)
    est = DecisionTreeClassifier()
    est.frozen = True
    assert clone(est) is est


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


def test_frozen_fit():
    class DummyEstimator(BaseEstimator):
        def fit(self, X, y, **kwargs):
            self.X_ = X
            self.y_ = y
            self.kwargs_ = kwargs
            self._last_call = 'fit'
            return self

        def fit_transform(self, X, y, **kwargs):
            self.X_ = X
            self.y_ = y
            self.kwargs_ = kwargs
            self._last_call = 'fit_transform'
            return np.array(self.X_[0]) + X

        def transform(self, X):
            self._last_call = 'transform'
            return np.array(self.X_[0]) + X

        def fit_wobble(self, X, y, **kwargs):
            self.X_ = X
            self.y_ = y
            self.kwargs_ = kwargs
            self._last_call = 'fit_wobble'
            return self.y_[0] + X[0][0]

        def wobble(self, X):
            self._last_call = 'wobble'
            return self.y_[0] + X[0][0]

    X_freeze = [[5]]
    y_freeze = [-1]
    z_freeze = [0]
    X_train = [[10]]
    y_train = [1]
    z_train = [0]

    for fit_method, method in [('fit', None),
                               ('fit_transform', 'transform'),
                               ('fit_wobble', 'wobble')]:

        # est is not frozen
        est = DummyEstimator().fit(X_freeze, y_freeze, z=z_freeze)

        result = frozen_fit(est, fit_method, X_train, y_train, z=z_train)
        # check it called .fit_transform(), not .fit().transform(), for example
        assert est._last_call == fit_method
        # check model was re-fit
        assert est.X_ == X_train
        assert est.y_ == y_train
        assert est.kwargs_ == {'z': z_train}
        if fit_method == 'fit':
            assert result is est
        else:
            assert_array_equal(result, getattr(est, method)(X_train))

        # est is not frozen
        est = DummyEstimator().fit(X_freeze, y_freeze, z=z_freeze)
        est.frozen = True
        result = frozen_fit(est, fit_method, X_train, y_train, z=z_train)
        # check model was not re-fit
        assert est.X_ == X_freeze
        assert est.y_ == y_freeze
        assert est.kwargs_ == {'z': z_freeze}
        if fit_method == 'fit':
            assert result is est
        else:
            assert est._last_call == method
            assert_array_equal(result, getattr(est, method)(X_train))
