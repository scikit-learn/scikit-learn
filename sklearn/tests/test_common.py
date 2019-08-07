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
import functools

import pytest

from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning, SkipTestWarning

import sklearn
from sklearn.base import RegressorMixin
from sklearn.cluster.bicluster import BiclusterMixin

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model.base import LinearClassifierMixin
from sklearn.linear_model import Ridge
from sklearn.utils import IS_PYPY
from sklearn.utils.estimator_checks import (
    _yield_all_checks,
    _safe_tags,
    set_checking_parameters,
    check_parameters_default_constructible,
    check_class_weight_balanced_linear_classifier)


def test_all_estimator_no_base_class():
    # test that all_estimators doesn't find abstract classes.
    for name, Estimator in all_estimators():
        msg = ("Base estimators such as {0} should not be included"
               " in all_estimators").format(name)
        assert not name.lower().startswith('base'), msg


@pytest.mark.parametrize(
        'name, Estimator',
        all_estimators()
)
def test_parameters_default_constructible(name, Estimator):
    # Test that estimators are default-constructible
    check_parameters_default_constructible(name, Estimator)


def _tested_estimators():
    for name, Estimator in all_estimators():
        if issubclass(Estimator, BiclusterMixin):
            continue
        if name.startswith("_"):
            continue
        # FIXME _skip_test should be used here (if we could)

        required_parameters = getattr(Estimator, "_required_parameters", [])
        if len(required_parameters):
            if required_parameters in (["estimator"], ["base_estimator"]):
                if issubclass(Estimator, RegressorMixin):
                    estimator = Estimator(Ridge())
                else:
                    estimator = Estimator(LinearDiscriminantAnalysis())
            else:
                warnings.warn("Can't instantiate estimator {} which requires "
                              "parameters {}".format(name,
                                                     required_parameters),
                              SkipTestWarning)
                continue
        else:
            estimator = Estimator()
        yield name, estimator


def _generate_checks_per_estimator(check_generator, estimators):
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        for name, estimator in estimators:
            for check in check_generator(name, estimator):
                yield estimator, check


def _rename_partial(val):
    if isinstance(val, functools.partial):
        kwstring = "".join(["{}={}".format(k, v)
                            for k, v in val.keywords.items()])
        return "{}({})".format(val.func.__name__, kwstring)
    # FIXME once we have short reprs we can use them here!
    if hasattr(val, "get_params") and not isinstance(val, type):
        return type(val).__name__


@pytest.mark.parametrize(
        "estimator, check",
        _generate_checks_per_estimator(_yield_all_checks,
                                       _tested_estimators()),
        ids=_rename_partial
)
def test_estimators(estimator, check):
    # Common tests for estimator instances
    with ignore_warnings(category=(DeprecationWarning, ConvergenceWarning,
                                   UserWarning, FutureWarning)):
        set_checking_parameters(estimator)
        name = estimator.__class__.__name__
        check(name, estimator)


@ignore_warnings(category=DeprecationWarning)
# ignore deprecated open(.., 'U') in numpy distutils
def test_configure():
    # Smoke test the 'configure' step of setup, this tests all the
    # 'configure' functions in the setup.pys in scikit-learn
    cwd = os.getcwd()
    setup_path = os.path.abspath(os.path.join(sklearn.__path__[0], '..'))
    setup_filename = os.path.join(setup_path, 'setup.py')
    if not os.path.exists(setup_filename):
        return
    try:
        os.chdir(setup_path)
        old_argv = sys.argv
        sys.argv = ['setup.py', 'config']

        # This test will run every setup.py and eventually call
        # check_openmp_support(), which tries to compile a C file that uses
        # OpenMP, unless SKLEARN_NO_OPENMP is set. Some users might want to run
        # the tests without having build-support for OpenMP. In particular, mac
        # users need to set some environment variables to build with openmp
        # support, and these might not be set anymore at test time. We thus
        # temporarily set SKLEARN_NO_OPENMP, so that this test runs smoothly.
        old_env = os.getenv('SKLEARN_NO_OPENMP')
        os.environ['SKLEARN_NO_OPENMP'] = "True"

        with warnings.catch_warnings():
            # The configuration spits out warnings when not finding
            # Blas/Atlas development headers
            warnings.simplefilter('ignore', UserWarning)
            with open('setup.py') as f:
                exec(f.read(), dict(__name__='__main__'))
    finally:
        sys.argv = old_argv
        if old_env is not None:
            os.environ['SKLEARN_NO_OPENMP'] = old_env
        else:
            del os.environ['SKLEARN_NO_OPENMP']
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
        if IS_PYPY and ('_svmlight_format' in modname or
                        'feature_extraction._hashing' in modname):
            continue
        package = __import__(modname, fromlist="dummy")
        for name in getattr(package, '__all__', ()):
            if getattr(package, name, None) is None:
                raise AttributeError(
                    "Module '{0}' has no attribute '{1}'".format(
                        modname, name))


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
