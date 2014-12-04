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

from sklearn.externals.six import PY3
from sklearn.utils.testing import assert_false, clean_warning_registry
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import _named_check

import sklearn
from sklearn.cluster.bicluster import BiclusterMixin
from sklearn.decomposition import ProjectedGradientNMF

from sklearn.linear_model.base import LinearClassifierMixin
from sklearn.utils.estimator_checks import (
    _yield_all_checks,
    CROSS_DECOMPOSITION,
    check_parameters_default_constructible,
    check_class_weight_balanced_linear_classifier,
    check_no_fit_attributes_set_in_init,
    check_transformer_n_iter,
    check_non_transformer_estimators_n_iter,
    check_get_params_invariance)


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
        for check in _yield_all_checks(name, Estimator):
            if issubclass(Estimator, ProjectedGradientNMF):
                # The ProjectedGradientNMF class is deprecated
                with ignore_warnings():
                    yield _named_check(check, name), name, Estimator
            else:
                yield _named_check(check, name), name, Estimator


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


# These are the known trailing-_ properties on current estimators.
#
# Such properties cause hasattr(estimator, "param_") to return True even if
# fit has not yet been called, making it harder to check whether the estimator
# has been properly fitted to data, but we allow these for backwards compat.
#
# If you're sure you need to add such a property, list it here.
_KNOWN_PROPERTIES = {
    'AdaBoostClassifier': ['feature_importances_'],
    'AdaBoostRegressor': ['feature_importances_'],
    'BernoulliNB': ['coef_', 'intercept_'],
    'DecisionTreeClassifier': ['feature_importances_'],
    'DecisionTreeRegressor': ['feature_importances_'],
    'ElasticNet': ['sparse_coef_'],
    'ExtraTreeClassifier': ['feature_importances_'],
    'ExtraTreeRegressor': ['feature_importances_'],
    'ExtraTreesClassifier': ['feature_importances_'],
    'ExtraTreesRegressor': ['feature_importances_'],
    'GradientBoostingClassifier': ['feature_importances_'],
    'GradientBoostingRegressor': ['feature_importances_'],
    'Lasso': ['sparse_coef_'],
    'MultinomialNB': ['coef_', 'intercept_'],
    'MultiTaskLasso': ['sparse_coef_'],
    'MultiTaskElasticNet': ['sparse_coef_'],
    'NuSVC': ['coef_'],
    'NuSVR': ['coef_'],
    'OneClassSVM': ['coef_'],
    'RandomForestClassifier': ['feature_importances_'],
    'RandomForestRegressor': ['feature_importances_'],
    'RidgeClassifier': ['classes_'],
    'RidgeClassifierCV': ['classes_'],
    'SpectralBiclustering': ['biclusters_'],
    'SpectralCoclustering': ['biclusters_'],
    'SVC': ['coef_'],
    'SVR': ['coef_'],
    'TfidfVectorizer': ['idf_'],
}


def test_no_fit_attributes_set_in_init():
    # Assert that estimators do not set fit attributes (ending in "_")
    # at initialization time.
    for name, Estimator in all_estimators():
        if name not in ['_BaseHMM', "GaussianHMM", "GMMHMM", 'MultinomialHMM']:
            yield (check_no_fit_attributes_set_in_init,
                   name, Estimator, _KNOWN_PROPERTIES.get(name, []))


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


def test_non_transformer_estimators_n_iter():
    # Test that all estimators of type which are non-transformer
    # and which have an attribute of max_iter, return the attribute
    # of n_iter atleast 1.
    for est_type in ['regressor', 'classifier', 'cluster']:
        regressors = all_estimators(type_filter=est_type)
        for name, Estimator in regressors:
            # LassoLars stops early for the default alpha=1.0 for
            # the iris dataset.
            if name == 'LassoLars':
                estimator = Estimator(alpha=0.)
            else:
                estimator = Estimator()
            if hasattr(estimator, "max_iter"):
                # These models are dependent on external solvers like
                # libsvm and accessing the iter parameter is non-trivial.
                if name in (['Ridge', 'SVR', 'NuSVR', 'NuSVC',
                             'RidgeClassifier', 'SVC', 'RandomizedLasso',
                             'LogisticRegressionCV']):
                    continue

                # Tested in test_transformer_n_iter below
                elif (name in CROSS_DECOMPOSITION or
                      name in ['LinearSVC', 'LogisticRegression']):
                    continue

                else:
                    # Multitask models related to ENet cannot handle
                    # if y is mono-output.
                    yield (_named_check(
                        check_non_transformer_estimators_n_iter, name),
                        name, estimator, 'Multi' in name)


def test_transformer_n_iter():
    transformers = all_estimators(type_filter='transformer')
    for name, Estimator in transformers:
        if issubclass(Estimator, ProjectedGradientNMF):
            # The ProjectedGradientNMF class is deprecated
            with ignore_warnings():
                estimator = Estimator()
        else:
            estimator = Estimator()
        # Dependent on external solvers and hence accessing the iter
        # param is non-trivial.
        external_solver = ['Isomap', 'KernelPCA', 'LocallyLinearEmbedding',
                           'RandomizedLasso', 'LogisticRegressionCV']

        if hasattr(estimator, "max_iter") and name not in external_solver:
            if isinstance(estimator, ProjectedGradientNMF):
                # The ProjectedGradientNMF class is deprecated
                with ignore_warnings():
                    yield _named_check(
                        check_transformer_n_iter, name), name, estimator
            else:
                yield _named_check(
                    check_transformer_n_iter, name), name, estimator


def test_get_params_invariance():
    # Test for estimators that support get_params, that
    # get_params(deep=False) is a subset of get_params(deep=True)
    # Related to issue #4465

    estimators = all_estimators(include_meta_estimators=False,
                                include_other=True)
    for name, Estimator in estimators:
        if hasattr(Estimator, 'get_params'):
            # If class is deprecated, ignore deprecated warnings
            if hasattr(Estimator.__init__, "deprecated_original"):
                with ignore_warnings():
                    yield _named_check(
                        check_get_params_invariance, name), name, Estimator
            else:
                yield _named_check(
                    check_get_params_invariance, name), name, Estimator
