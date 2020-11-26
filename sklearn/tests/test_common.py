"""
General tests for all estimators in sklearn.
"""

# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD 3 clause

import os
import warnings
import sys
import re
import pkgutil
from inspect import isgenerator
from functools import partial

import numpy as np
import pytest

from sklearn.utils import all_estimators
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.estimator_checks import check_estimator

import sklearn
from sklearn.base import BiclusterMixin

from sklearn.decomposition import PCA
from sklearn.linear_model._base import LinearClassifierMixin
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.pipeline import make_pipeline

from sklearn.utils import IS_PYPY
from sklearn.utils._testing import SkipTest
from sklearn.utils.estimator_checks import (
    _construct_instance,
    _set_checking_parameters,
    _get_check_estimator_ids,
    check_class_weight_balanced_linear_classifier,
    parametrize_with_checks,
    check_n_features_in_after_fitting,
)
from sklearn.multiclass import check_classification_targets
from sklearn.utils.validation import (
    check_array,
    check_is_fitted,
    check_X_y,
)


def test_all_estimator_no_base_class():
    # test that all_estimators doesn't find abstract classes.
    for name, Estimator in all_estimators():
        msg = ("Base estimators such as {0} should not be included"
               " in all_estimators").format(name)
        assert not name.lower().startswith('base'), msg


def _sample_func(x, y=1):
    pass


@pytest.mark.parametrize("val, expected", [
    (partial(_sample_func, y=1), "_sample_func(y=1)"),
    (_sample_func, "_sample_func"),
    (partial(_sample_func, 'world'), "_sample_func"),
    (LogisticRegression(C=2.0), "LogisticRegression(C=2.0)"),
    (LogisticRegression(random_state=1, solver='newton-cg',
                        class_weight='balanced', warm_start=True),
     "LogisticRegression(class_weight='balanced',random_state=1,"
     "solver='newton-cg',warm_start=True)")
])
def test_get_check_estimator_ids(val, expected):
    assert _get_check_estimator_ids(val) == expected


def _tested_estimators():
    for name, Estimator in all_estimators():
        if issubclass(Estimator, BiclusterMixin):
            continue
        try:
            estimator = _construct_instance(Estimator)
        except SkipTest:
            continue

        yield estimator


@parametrize_with_checks(list(_tested_estimators()))
def test_estimators(estimator, check, request):
    # Common tests for estimator instances
    with ignore_warnings(category=(FutureWarning,
                                   ConvergenceWarning,
                                   UserWarning, FutureWarning)):
        _set_checking_parameters(estimator)
        check(estimator)


def test_check_estimator_generate_only():
    all_instance_gen_checks = check_estimator(LogisticRegression(),
                                              generate_only=True)
    assert isgenerator(all_instance_gen_checks)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
# ignore deprecated open(.., 'U') in numpy distutils
def test_configure():
    # Smoke test the 'configure' step of setup, this tests all the
    # 'configure' functions in the setup.pys in scikit-learn
    # This test requires Cython which is not necessarily there when running
    # the tests of an installed version of scikit-learn or when scikit-learn
    # is installed in editable mode by pip build isolation enabled.
    pytest.importorskip("Cython")
    cwd = os.getcwd()
    setup_path = os.path.abspath(os.path.join(sklearn.__path__[0], '..'))
    setup_filename = os.path.join(setup_path, 'setup.py')
    if not os.path.exists(setup_filename):
        pytest.skip('setup.py not available')
    # XXX unreached code as of v0.22
    try:
        os.chdir(setup_path)
        old_argv = sys.argv
        sys.argv = ['setup.py', 'config']

        with warnings.catch_warnings():
            # The configuration spits out warnings when not finding
            # Blas/Atlas development headers
            warnings.simplefilter('ignore', UserWarning)
            with open('setup.py') as f:
                exec(f.read(), dict(__name__='__main__'))
    finally:
        sys.argv = old_argv
        os.chdir(cwd)


def _tested_linear_classifiers():
    classifiers = all_estimators(type_filter='classifier')

    with warnings.catch_warnings(record=True):
        for name, clazz in classifiers:
            required_parameters = getattr(clazz, "_required_parameters", [])
            if len(required_parameters):
                # FIXME
                continue

            if ('class_weight' in clazz().get_params().keys() and
                    issubclass(clazz, LinearClassifierMixin)):
                yield name, clazz


@pytest.mark.parametrize("name, Classifier",
                         _tested_linear_classifiers())
def test_class_weight_balanced_linear_classifiers(name, Classifier):
    check_class_weight_balanced_linear_classifier(name, Classifier)


@ignore_warnings
def test_import_all_consistency():
    # Smoke test to check that any name in a __all__ list is actually defined
    # in the namespace of the module or package.
    pkgs = pkgutil.walk_packages(path=sklearn.__path__, prefix='sklearn.',
                                 onerror=lambda _: None)
    submods = [modname for _, modname, _ in pkgs]
    for modname in submods + ['sklearn']:
        if ".tests." in modname:
            continue
        if IS_PYPY and ('_svmlight_format_io' in modname or
                        'feature_extraction._hashing_fast' in modname):
            continue
        package = __import__(modname, fromlist="dummy")
        for name in getattr(package, '__all__', ()):
            assert hasattr(package, name),\
                "Module '{0}' has no attribute '{1}'".format(modname, name)


def test_root_import_all_completeness():
    EXCEPTIONS = ('utils', 'tests', 'base', 'setup', 'conftest')
    for _, modname, _ in pkgutil.walk_packages(path=sklearn.__path__,
                                               onerror=lambda _: None):
        if '.' in modname or modname.startswith('_') or modname in EXCEPTIONS:
            continue
        assert modname in sklearn.__all__


def test_all_tests_are_importable():
    # Ensure that for each contentful subpackage, there is a test directory
    # within it that is also a subpackage (i.e. a directory with __init__.py)

    HAS_TESTS_EXCEPTIONS = re.compile(r'''(?x)
                                      \.externals(\.|$)|
                                      \.tests(\.|$)|
                                      \._
                                      ''')
    lookup = {name: ispkg
              for _, name, ispkg
              in pkgutil.walk_packages(sklearn.__path__, prefix='sklearn.')}
    missing_tests = [name for name, ispkg in lookup.items()
                     if ispkg
                     and not HAS_TESTS_EXCEPTIONS.search(name)
                     and name + '.tests' not in lookup]
    assert missing_tests == [], ('{0} do not have `tests` subpackages. '
                                 'Perhaps they require '
                                 '__init__.py or an add_subpackage directive '
                                 'in the parent '
                                 'setup.py'.format(missing_tests))


def test_class_support_removed():
    # Make sure passing classes to check_estimator or parametrize_with_checks
    # raises an error

    msg = "Passing a class was deprecated.* isn't supported anymore"
    with pytest.raises(TypeError, match=msg):
        check_estimator(LogisticRegression)

    with pytest.raises(TypeError, match=msg):
        parametrize_with_checks([LogisticRegression])


def _generate_search_cv_instances():
    for SearchCV, (Estimator, param_grid) in zip(
        [GridSearchCV, RandomizedSearchCV],
        [
            (Ridge, {"alpha": [0.1, 1.0]}),
            (LogisticRegression, {"C": [0.1, 1.0]}),
        ],
    ):
        yield SearchCV(Estimator(), param_grid)

    for SearchCV, (Estimator, param_grid) in zip(
        [GridSearchCV, RandomizedSearchCV],
        [
            (Ridge, {"ridge__alpha": [0.1, 1.0]}),
            (LogisticRegression, {"logisticregression__C": [0.1, 1.0]}),
        ],
    ):
        yield SearchCV(
            make_pipeline(PCA(), Estimator()), param_grid
        ).set_params(error_score="raise")


@parametrize_with_checks(list(_generate_search_cv_instances()))
def test_search_cv(estimator, check, request):
    # Common tests for SearchCV instances
    # We have a separate test because those meta-estimators can accept a
    # wide range of base estimators (classifiers, regressors, pipelines)
    with ignore_warnings(
        category=(
            FutureWarning,
            ConvergenceWarning,
            UserWarning,
            FutureWarning,
        )
    ):
        check(estimator)


# TODO: When more modules get added, we can remove it from this list to make
# sure it gets tested. After we finish each module we can move the checks
# into sklearn.utils.estimator_checks.check_n_features_in.
#
# check_estimators_partial_fit_n_features can either be removed or updated
# with the two more assertions:
# 1. `n_features_in_` is set during the first call to `partial_fit`.
# 2. More strict when it comes to the error message.
#
# check_classifiers_train would need to be updated with the error message
N_FEATURES_IN_AFTER_FIT_MODULES_TO_IGNORE = {
    'calibration',
    'compose',
    'covariance',
    'cross_decomposition',
    'discriminant_analysis',
    'ensemble',
    'feature_extraction',
    'feature_selection',
    'gaussian_process',
    'isotonic',
    'linear_model',
    'manifold',
    'mixture',
    'model_selection',
    'multiclass',
    'multioutput',
    'naive_bayes',
    'neighbors',
    'pipeline',
    'random_projection',
    'semi_supervised',
    'svm',
}

N_FEATURES_IN_AFTER_FIT_ESTIMATORS = [
    est for est in _tested_estimators() if est.__module__.split('.')[1] not in
    N_FEATURES_IN_AFTER_FIT_MODULES_TO_IGNORE
]


@pytest.mark.parametrize("estimator", N_FEATURES_IN_AFTER_FIT_ESTIMATORS,
                         ids=_get_check_estimator_ids)
def test_check_n_features_in_after_fitting(estimator):
    _set_checking_parameters(estimator)
    check_n_features_in_after_fitting(estimator.__class__.__name__, estimator)


class MinimalClassifier:
    """Minimal classifier implementation with inheriting from BaseEstimator."""
    _estimator_type = "classifier"

    def __init__(self, param=None):
        self.param = param

    def __repr__(self):
        # Only required when using pytest-xdist to get an id not associated
        # with the memory location
        return self.__class__.__name__

    def get_params(self, deep=True):
        return {"param": self.param}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def __getstate__(self):
        return self.__dict__.copy()

    def __setstate__(self, state):
        self.__dict__.update(state)

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        check_classification_targets(y)
        self.n_features_in_ = X.shape[1]
        self.classes_, counts = np.unique(y, return_counts=True)
        self._most_frequent_class_idx = counts.argmax()
        return self

    def predict_proba(self, X):
        check_is_fitted(self)
        X = check_array(X)
        if X.shape[1] != self.n_features_in_:
            raise ValueError
        proba_shape = (X.shape[0], self.classes_.size)
        y_proba = np.zeros(shape=proba_shape, dtype=np.float64)
        y_proba[:, self._most_frequent_class_idx] = 1.0
        return y_proba

    def predict(self, X):
        y_proba = self.predict_proba(X)
        y_pred = y_proba.argmax(axis=1)
        return self.classes_[y_pred]

    def score(self, X, y):
        from sklearn.metrics import accuracy_score
        return accuracy_score(y, self.predict(X))


class MinimalRegressor:
    """Minimal regressor implementation with inheriting from BaseEstimator."""
    _estimator_type = "regressor"

    def __init__(self, param=None):
        self.param = param

    def __repr__(self):
        # Only required when using pytest-xdist to get an id not associated
        # with the memory location
        return self.__class__.__name__

    def get_params(self, deep=True):
        return {"param": self.param}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def __getstate__(self):
        return self.__dict__.copy()

    def __setstate__(self, state):
        self.__dict__.update(state)

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        self.n_features_in_ = X.shape[1]
        self._mean = np.mean(y)
        return self

    def predict(self, X):
        check_is_fitted(self)
        X = check_array(X)
        if X.shape[1] != self.n_features_in_:
            raise ValueError
        return np.ones(shape=(X.shape[0],)) * self._mean

    def score(self, X, y):
        from sklearn.metrics import r2_score
        return r2_score(y, self.predict(X))


class MinimalTransformer:
    """Minimal transformer implementation with inheriting from
    BaseEstimator."""

    def __init__(self, param=None):
        self.param = param

    def __repr__(self):
        # Only required when using pytest-xdist to get an id not associated
        # with the memory location
        return self.__class__.__name__

    def get_params(self, deep=True):
        return {"param": self.param}

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        return self

    def __getstate__(self):
        return self.__dict__.copy()

    def __setstate__(self, state):
        self.__dict__.update(state)

    def fit(self, X, y=None):
        X = check_array(X)
        self.n_features_in_ = X.shape[1]
        return self

    def transform(self, X, y=None):
        check_is_fitted(self)
        X = check_array(X)
        if X.shape[1] != self.n_features_in_:
            raise ValueError
        return X

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X, y)


def _generate_minimal_compatible_instances():
    """Generate instance containing estimators from minimal class compatible
    implementation."""
    for SearchCV, (Estimator, param_grid) in zip(
        [GridSearchCV, RandomizedSearchCV],
        [
            (MinimalRegressor, {"param": [1, 10]}),
            (MinimalClassifier, {"param": [1, 10]}),
        ],
    ):
        yield SearchCV(Estimator(), param_grid)

    for SearchCV, (Estimator, param_grid) in zip(
        [GridSearchCV, RandomizedSearchCV],
        [
            (MinimalRegressor, {
                "minimaltransformer__param": [1, 10],
                "minimalregressor__param": [1, 10],
            }),
            (MinimalClassifier, {
                "minimaltransformer__param": [1, 10],
                "minimalclassifier__param": [1, 10],
            }),
        ],
    ):
        yield SearchCV(
            make_pipeline(MinimalTransformer(), Estimator()), param_grid
        ).set_params(error_score="raise")


# FIXME: hopefully in 0.25
@pytest.mark.xfail(
    reason=("This test is currently failing because checks are granular "
            "enough. Once checks are split with some kind of only API tests, "
            "this test should enabled.")
)
@parametrize_with_checks(list(_generate_minimal_compatible_instances()))
def test_minimal_class_implementation_checks(estimator, check):
    # Check that third-party library can run tests without inheriting from
    # BaseEstimator.
    check(estimator)
