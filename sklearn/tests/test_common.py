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
    estimators_to_check = ['GaussianNB', 'BernoulliNB', 'MultinomialNB',
                           'DecisionTreeClassifier']
    for name, estimator in all_estimators():
        if estimator in estimators_to_check:
            clf1 = estimator(classes=[1, 2]).fit(X, y)
            pred1 = clf1.predict_proba(X)

            clf2 = estimator(classes=[1, 2, 3, 4]).fit(X, y)
            pred2 = clf2.predict_proba(X)

            clf3 = estimator(classes=[0, 1, 2, 3]).fit(X, y)
            pred3 = clf3.predict_proba(X)

            # check shapes:
            assert_equal(pred1.shape, (X.shape[0], 2))
            assert_equal(pred2.shape, (X.shape[0], 4))
            assert_equal(pred3.shape, (X.shape[0], 4))

            # check same columns are equal:
            assert_array_almost_equal(pred1[:, 0], pred2[:, 0])
            assert_array_almost_equal(pred1[:, 0], pred3[:, 1])
            assert_array_almost_equal(pred1[:, 1], pred2[:, 1])
            assert_array_almost_equal(pred1[:, 1], pred3[:, 2])
            assert_array_almost_equal(pred2[:, 2], 0)
            assert_array_almost_equal(pred2[:, 3], 0)
            assert_array_almost_equal(pred3[:, 0], 0)
            assert_array_almost_equal(pred3[:, 3], 0)

    # Test whether duplicate classes result in error
            expected_msg = ("Classses parameter should contain all unique"
                            " values, duplicates found in [1, 1, 2]")
            expected_msg2 = ("Classses parameter should contain sorted values"
                             ", unsorted values found in [1, 3, 2]")

            assert_raise_message(ValueError, expected_msg,
                                 estimator(classes=[1, 1, 2]).fit, X, y)
            assert_raise_message(ValueError, expected_msg,
                                 estimator().partial_fit, X, y,
                                 classes=[1, 1, 2])
            assert_raise_message(ValueError, expected_msg2,
                                 estimator().partial_fit, X, y,
                                 classes=[1, 3, 2])
