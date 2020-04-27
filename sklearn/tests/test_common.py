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

import pytest


from sklearn.utils import all_estimators
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.estimator_checks import check_estimator

import sklearn
from sklearn.base import BiclusterMixin

from sklearn.linear_model._base import LinearClassifierMixin
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.utils import IS_PYPY
from sklearn.utils._testing import SkipTest
from sklearn.utils.estimator_checks import (
    _mark_xfail_checks,
    _construct_instance,
    _set_checking_parameters,
    _set_check_estimator_ids,
    check_parameters_default_constructible,
    check_class_weight_balanced_linear_classifier,
    parametrize_with_checks)


def test_all_estimator_no_base_class():
    # test that all_estimators doesn't find abstract classes.
    for name, Estimator in all_estimators():
        msg = ("Base estimators such as {0} should not be included"
               " in all_estimators").format(name)
        assert not name.lower().startswith('base'), msg


@ignore_warnings("Passing a class is depr", category=FutureWarning)  # 0.24
def test_estimator_cls_parameterize_with_checks():
    # TODO: remove test in 0.24
    # Non-regression test for #16707 to ensure that parametrize_with_checks
    # works with estimator classes
    param_checks = parametrize_with_checks([LogisticRegression])
    # Using the generator does not raise
    list(param_checks.args[1])


def test_mark_xfail_checks_with_unconsructable_estimator():
    class MyEstimator:
        def __init__(self):
            raise ValueError("This is bad")

    estimator, check = _mark_xfail_checks(MyEstimator, 42, None)
    assert estimator == MyEstimator
    assert check == 42


@pytest.mark.parametrize(
        'name, Estimator',
        all_estimators()
)
def test_parameters_default_constructible(name, Estimator):
    # Test that estimators are default-constructible
    check_parameters_default_constructible(name, Estimator)


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
def test_set_check_estimator_ids(val, expected):
    assert _set_check_estimator_ids(val) == expected


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


@ignore_warnings("Passing a class is depr", category=FutureWarning)  # 0.24
def test_check_estimator_generate_only():
    # TODO in 0.24: remove checks on passing a class
    estimator_cls_gen_checks = check_estimator(LogisticRegression,
                                               generate_only=True)
    all_instance_gen_checks = check_estimator(LogisticRegression(),
                                              generate_only=True)
    assert isgenerator(estimator_cls_gen_checks)
    assert isgenerator(all_instance_gen_checks)

    estimator_cls_checks = list(estimator_cls_gen_checks)
    all_instance_checks = list(all_instance_gen_checks)

    # all classes checks include check_parameters_default_constructible
    assert len(estimator_cls_checks) == len(all_instance_checks) + 1

    # TODO: meta-estimators like GridSearchCV has required parameters
    # that do not have default values. This is expected to change in the future
    with pytest.raises(SkipTest):
        for estimator, check in check_estimator(GridSearchCV,
                                                generate_only=True):
            check(estimator)


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


# TODO: remove in 0.24
def test_class_support_deprecated():
    # Make sure passing classes to check_estimator or parametrize_with_checks
    # is deprecated

    msg = "Passing a class is deprecated"
    with pytest.warns(FutureWarning, match=msg):
        check_estimator(LogisticRegression)

    with pytest.warns(FutureWarning, match=msg):
        parametrize_with_checks([LogisticRegression])

    # Make sure check_parameters_default_constructible accepts instances now
    check_parameters_default_constructible('name', LogisticRegression())
