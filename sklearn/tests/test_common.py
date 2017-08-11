"""
General tests for all estimators in sklearn.
"""

# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD 3 clause
from __future__ import print_function

import os
import warnings
import sys
import re
import pkgutil
import numpy as np

from sklearn.externals.six import PY3
from sklearn.utils.testing import assert_false, clean_warning_registry
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import _named_check
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raise_message

import sklearn
from sklearn.cluster.bicluster import BiclusterMixin

from sklearn.linear_model.base import LinearClassifierMixin
from sklearn.utils.estimator_checks import (
    _yield_all_checks,
    set_checking_parameters,
    check_parameters_default_constructible,
    check_no_fit_attributes_set_in_init,
    check_class_weight_balanced_linear_classifier)


def test_all_estimator_no_base_class():
    # test that all_estimators doesn't find abstract classes.
    for name, Estimator in all_estimators():
        msg = ("Base estimators such as {0} should not be included"
               " in all_estimators").format(name)
        assert_false(name.lower().startswith('base'), msg=msg)


def test_all_estimators():
    # Test that estimators are default-constructible, cloneable
    # and have working repr.
    estimators = all_estimators(include_meta_estimators=True)

    # Meta sanity-check to make sure that the estimator introspection runs
    # properly
    assert_greater(len(estimators), 0)

    for name, Estimator in estimators:
        # some can just not be sensibly default constructed
        yield (_named_check(check_parameters_default_constructible, name),
               name, Estimator)


def test_non_meta_estimators():
    # input validation etc for non-meta estimators
    estimators = all_estimators()
    for name, Estimator in estimators:
        if issubclass(Estimator, BiclusterMixin):
            continue
        if name.startswith("_"):
            continue
        estimator = Estimator()
        # check this on class
        yield _named_check(
            check_no_fit_attributes_set_in_init, name), name, Estimator

        for check in _yield_all_checks(name, estimator):
            set_checking_parameters(estimator)
            yield _named_check(check, name), name, estimator


def test_configure():
    # Smoke test the 'configure' step of setup, this tests all the
    # 'configure' functions in the setup.pys in the scikit
    cwd = os.getcwd()
    setup_path = os.path.abspath(os.path.join(sklearn.__path__[0], '..'))
    setup_filename = os.path.join(setup_path, 'setup.py')
    if not os.path.exists(setup_filename):
        return
    try:
        os.chdir(setup_path)
        old_argv = sys.argv
        sys.argv = ['setup.py', 'config']
        clean_warning_registry()
        with warnings.catch_warnings():
            # The configuration spits out warnings when not finding
            # Blas/Atlas development headers
            warnings.simplefilter('ignore', UserWarning)
            if PY3:
                with open('setup.py') as f:
                    exec(f.read(), dict(__name__='__main__'))
            else:
                execfile('setup.py', dict(__name__='__main__'))
    finally:
        sys.argv = old_argv
        os.chdir(cwd)


def test_class_weight_balanced_linear_classifiers():
    classifiers = all_estimators(type_filter='classifier')

    clean_warning_registry()
    with warnings.catch_warnings(record=True):
        linear_classifiers = [
            (name, clazz)
            for name, clazz in classifiers
            if ('class_weight' in clazz().get_params().keys() and
                issubclass(clazz, LinearClassifierMixin))]

    for name, Classifier in linear_classifiers:
        yield _named_check(check_class_weight_balanced_linear_classifier,
                           name), name, Classifier


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
        package = __import__(modname, fromlist="dummy")
        for name in getattr(package, '__all__', ()):
            if getattr(package, name, None) is None:
                raise AttributeError(
                    "Module '{0}' has no attribute '{1}'".format(
                        modname, name))


def test_root_import_all_completeness():
    EXCEPTIONS = ('utils', 'tests', 'base', 'setup')
    for _, modname, _ in pkgutil.walk_packages(path=sklearn.__path__,
                                               onerror=lambda _: None):
        if '.' in modname or modname.startswith('_') or modname in EXCEPTIONS:
            continue
        assert_in(modname, sklearn.__all__)


def test_all_tests_are_importable():
    # Ensure that for each contentful subpackage, there is a test directory
    # within it that is also a subpackage (i.e. a directory with __init__.py)

    HAS_TESTS_EXCEPTIONS = re.compile(r'''(?x)
                                      \.externals(\.|$)|
                                      \.tests(\.|$)|
                                      \._
                                      ''')
    lookup = dict((name, ispkg)
                  for _, name, ispkg
                  in pkgutil.walk_packages(sklearn.__path__,
                                           prefix='sklearn.'))
    missing_tests = [name for name, ispkg in lookup.items()
                     if ispkg
                     and not HAS_TESTS_EXCEPTIONS.search(name)
                     and name + '.tests' not in lookup]
    assert_equal(missing_tests, [],
                 '{0} do not have `tests` subpackages. Perhaps they require '
                 '__init__.py or an add_subpackage directive in the parent '
                 'setup.py'.format(missing_tests))


@ignore_warnings(category=(DeprecationWarning))
def test_classes_parameter_extra_classes():
    # Test whether adding extra classes doesn't change the prediction of the
    # existing classes just adds 0 value columns for new classes
    rng = np.random.RandomState(0)
    X = rng.randint(5, size=(6, 100))
    y = np.array([1, 1, 1, 2, 2, 2])
    estimators_to_check = set(['GaussianNB',
                               'BernoulliNB',
                               'MultinomialNB',
                               'DecisionTreeClassifier',
                               'ExtraTreeClassifier',
                               'RandomForestClassifier',
                               'ExtraTreesClassifier',
                               'GradientBoostingClassifier',
                               'LinearSVC',
                               # 'SVC',
                               # 'NuSVC'
                               ])
    num_predict_proba_checks = 8
    num_decision_function_checks = 2

    for name, estimator in all_estimators():
        if name in estimators_to_check:
            # remove the estimator from set:
            estimators_to_check.remove(name)
            params = {}
            if hasattr(estimator(), 'random_state'):
                params['random_state'] = 0
            if hasattr(estimator(), 'probability'):
                params['probability'] = True

            if hasattr(estimator(), 'predict_proba'):
                num_predict_proba_checks -= 1
                clf = estimator()

                clf.set_params(**params)
                pred0 = clf.fit(X, y).predict_proba(X)

                params['classes'] = [1, 2]
                clf.set_params(**params)
                pred1 = clf.fit(X, y).predict_proba(X)

                params['classes'] = [1, 2, 3, 4]
                clf.set_params(**params)
                pred2 = clf.fit(X, y).predict_proba(X)

                params['classes'] = [0, 1, 2, 3]
                clf.set_params(**params)
                pred3 = clf.fit(X, y).predict_proba(X)

                # check shapes:
                assert_equal(pred0.shape, (X.shape[0], 2))
                assert_equal(pred1.shape, (X.shape[0], 2))
                assert_equal(pred2.shape, (X.shape[0], 4))
                assert_equal(pred3.shape, (X.shape[0], 4))

                # check same columns are equal. since these are probabilities,
                # checking upto 3 decimal places should be fine.
                assert_array_almost_equal(pred0, pred1, 3)
                assert_array_almost_equal(pred1[:, 0:2], pred2[:, 0:2], 3)
                assert_array_almost_equal(pred1[:, 0:2], pred3[:, 1:3], 3)
                assert_array_almost_equal(pred2[:, 2:4], 0, 3)
                assert_array_almost_equal(pred3[:, [0, 3]], 0, 3)

            if hasattr(estimator, 'decision_function'):
                num_decision_function_checks -= 1
                clf = estimator()
                params.pop("classes", None)

                clf.set_params(**params)
                pred0 = clf.fit(X, y).decision_function(X)
                pred0 = np.atleast_2d(pred0).T

                params['classes'] = [1, 2]
                clf.set_params(**params)
                pred1 = clf.fit(X, y).decision_function(X)
                pred1 = np.atleast_2d(pred1).T

                params['classes'] = [1, 2, 3, 4]
                clf.set_params(**params)
                pred2 = clf.fit(X, y).decision_function(X)

                params['classes'] = [0, 1, 2, 3]
                clf.set_params(**params)
                pred3 = clf.fit(X, y).decision_function(X)

                # check shapes:
                assert_equal(pred0.shape, (X.shape[0], 1))
                assert_equal(pred1.shape, (X.shape[0], 1))
                assert_equal(pred2.shape, (X.shape[0], 4))
                assert_equal(pred3.shape, (X.shape[0], 4))

                # check same columns have same sign
                assert_array_almost_equal(pred0, pred1, 1)
                assert_array_almost_equal(pred2[:, 0:2],
                                          pred3[:, 1:3], 1)
                assert_array_almost_equal(pred2[:, 2:4],
                                          pred3[:, [0, 3]], 1)

            # Test whether duplicate classes result in error
            expected_msg = ("Classses parameter should contain all unique"
                            " values, duplicates found in [1 1 2]")
            expected_msg2 = ("Classses parameter should contain sorted values"
                             ", unsorted values found in [1 3 2]")
            expected_msg3 = ("The target label(s) [1] in y do not exist in the"
                             "initial classes [2, 3]")

            assert_raise_message(ValueError, expected_msg,
                                 estimator(classes=[1, 1, 2]).fit, X, y)
            assert_raise_message(ValueError, expected_msg2,
                                 estimator(classes=[1, 3, 2]).fit, X, y)
            assert_raise_message(ValueError, expected_msg3,
                                 estimator(classes=[2, 3]).fit, X, y)
            # Run only if partial_fit is present:
            if hasattr(estimator, 'partial_fit'):
                assert_raise_message(ValueError, expected_msg,
                                     estimator().partial_fit, X, y,
                                     classes=[1, 1, 2])
                assert_raise_message(ValueError, expected_msg2,
                                     estimator().partial_fit, X, y,
                                     classes=[1, 3, 2])
                assert_raise_message(ValueError, expected_msg3,
                                     estimator().partial_fit, X, y,
                                     classes=[2, 3])

    # check that all estimators in estimator_to_check were tested and right
    # number of predict_proba and decision_function checks were done
    assert_equal(estimators_to_check, set(),
                 "All estimators specified in list estimators_to_check were"
                 "not checked")
    assert_equal(num_predict_proba_checks, 0,
                 "The required number of predict_proba checks not done, "
                 "remaining: %d" % num_predict_proba_checks)
    assert_equal(num_decision_function_checks, 0,
                 "The required number of decision_function checks not done, "
                 "remaining: %d" % num_decision_function_checks)


@ignore_warnings(category=(DeprecationWarning))
def test_classes_parameter_extra_classes_multilabel():
    # Test same as extra_classes but for multilabel
    rng = np.random.RandomState(0)
    X = rng.randint(5, size=(6, 100))
    y = np.array([1, 1, 1, 2, 2, 2])
    y = np.hstack([y, y]).reshape(len(y), 2)
    estimators_to_check = set(['DecisionTreeClassifier',
                               'ExtraTreeClassifier',
                               'RandomForestClassifier',
                               'ExtraTreesClassifier'])
    num_predict_proba_checks = 4

    for name, estimator in all_estimators():
        if name in estimators_to_check:
            # remove the estimator from set:
            estimators_to_check.remove(name)

            clf = estimator()
            if hasattr(clf, 'random_state'):
                params = {'random_state': 0}
            else:
                params = {}

            if hasattr(clf, 'predict_proba'):
                num_predict_proba_checks -= 1

                clf0 = estimator().fit(X, y)
                pred0 = clf0.predict_proba(X)

                clf.set_params(**params)
                pred0 = clf.fit(X, y).predict_proba(X)

                params['classes'] = [[1, 2], [1, 2]]
                clf.set_params(**params)
                pred1 = clf.fit(X, y).predict_proba(X)

                params['classes'] = [[1, 2, 3, 4], [1, 2, 3, 4]]
                clf.set_params(**params)
                pred2 = clf.fit(X, y).predict_proba(X)

                params['classes'] = [[0, 1, 2, 3], [0, 1, 2, 3]]
                clf.set_params(**params)
                pred3 = clf.fit(X, y).predict_proba(X)

                for i in range(2):
                    pred0_i = pred0[i]
                    pred1_i = pred1[i]
                    pred2_i = pred2[i]
                    pred3_i = pred3[i]
                    # check shapes:
                    assert_equal(pred0_i.shape, (X.shape[0], 2))
                    assert_equal(pred1_i.shape, (X.shape[0], 2))
                    assert_equal(pred2_i.shape, (X.shape[0], 4))
                    assert_equal(pred3_i.shape, (X.shape[0], 4))

                    # check same columns are equal.
                    assert_array_almost_equal(pred0_i, pred1_i, 3)
                    assert_array_almost_equal(pred1_i[:, 0:2],
                                              pred2_i[:, 0:2], 3)
                    assert_array_almost_equal(pred1_i[:, 0:2],
                                              pred3_i[:, 1:3], 3)
                    assert_array_almost_equal(pred2_i[:, 2:4], 0, 3)
                    assert_array_almost_equal(pred3_i[:, [0, 3]], 0, 3)

    # check that all estimators in estimator_to_check were tested
    assert_equal(estimators_to_check, set(),
                 "All estimators specified in list estimators_to_check were"
                 "not checked")
    assert_equal(num_predict_proba_checks, 0,
                 "The required number of predict_proba checks not done, "
                 "remaining: %d" % num_predict_proba_checks)
