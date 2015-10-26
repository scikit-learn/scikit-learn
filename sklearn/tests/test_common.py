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
import pkgutil
import numpy as np

from sklearn import datasets
from sklearn.base import is_classifier, is_regressor
from sklearn.cross_validation import train_test_split
from sklearn.externals.six import PY3
from sklearn.externals.six.moves import zip
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_false, clean_warning_registry
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import set_random_state
from sklearn.utils.validation import has_fit_parameter

from numpy.testing import assert_array_almost_equal, assert_equal

import sklearn
from sklearn.cluster.bicluster import BiclusterMixin
from sklearn.decomposition import ProjectedGradientNMF

from sklearn.linear_model.base import LinearClassifierMixin
from sklearn.utils.estimator_checks import (
    _yield_all_checks,
    CROSS_DECOMPOSITION,
    check_parameters_default_constructible,
    check_class_weight_balanced_linear_classifier,
    check_transformer_n_iter,
    check_non_transformer_estimators_n_iter,
    check_get_params_invariance,
    check_fit2d_predict1d,
    check_fit1d_1sample)


def test_all_estimator_no_base_class():
    # test that all_estimators doesn't find abstract classes.
    for name, Estimator in all_estimators():
        msg = ("Base estimators such as {0} should not be included"
               " in all_estimators").format(name)
        assert_false(name.lower().startswith('base'), msg=msg)


def test_all_estimators():
    # Test that estimators are default-constructible, clonable
    # and have working repr.
    estimators = all_estimators(include_meta_estimators=True)

    # Meta sanity-check to make sure that the estimator introspection runs
    # properly
    assert_greater(len(estimators), 0)

    for name, Estimator in estimators:
        # some can just not be sensibly default constructed
        yield check_parameters_default_constructible, name, Estimator


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
                    yield check, name, Estimator
            else:
                yield check, name, Estimator

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
            if 'class_weight' in clazz().get_params().keys()
               and issubclass(clazz, LinearClassifierMixin)]

    for name, Classifier in linear_classifiers:
        if name == "LogisticRegressionCV":
            # Contrary to RidgeClassifierCV, LogisticRegressionCV use actual
            # CV folds and fit a model for each CV iteration before averaging
            # the coef. Therefore it is expected to not behave exactly as the
            # other linear model.
            continue
        yield check_class_weight_balanced_linear_classifier, name, Classifier


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
                    yield (check_non_transformer_estimators_n_iter,
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
                    yield check_transformer_n_iter, name, estimator
            else:
                yield check_transformer_n_iter, name, estimator


def test_get_params_invariance():
    # Test for estimators that support get_params, that
    # get_params(deep=False) is a subset of get_params(deep=True)
    # Related to issue #4465

    estimators = all_estimators(include_meta_estimators=False, include_other=True)
    for name, Estimator in estimators:
        if hasattr(Estimator, 'get_params'):
            # The ProjectedGradientNMF class is deprecated
            if issubclass(Estimator, ProjectedGradientNMF):
                with ignore_warnings():
                    yield check_get_params_invariance, name, Estimator
            else:
                yield check_get_params_invariance, name, Estimator
            #yield check_transformer_n_iter, name, Estimator()


def test_sample_weight_consistency(random_state=42):
    exclude = [
               'AdaBoostClassifier', 'AdaBoostRegressor', 
               'BaggingClassifier', 'BaggingRegressor',
               'GradientBoostingClassifier', 'GradientBoostingRegressor', 
               'ExtraTreeClassifier', 'ExtraTreeRegressor',
               'ExtraTreesClassifier', 'ExtraTreesRegressor',
               'DecisionTreeClassifier','DecisionTreeRegressor',
               'LogisticRegression', 'LogisticRegressionCV', 'LinearSVC',
               'MultinomialNB',  # Requires positive samples
               'CalibratedClassifierCV',  # This is a meta-estimator using LinearSVC
               'SGDClassifier',  # Doesn't work. Probably more data needed
               'SGDRegressor',  # Doesn't work. Probably more data needed
               'Perceptron',  # Uses SGD too. Doesn't work. Probably more data needed
               'RidgeClassifierCV', 'RidgeCV', 
               'RandomForestClassifier', 'RandomForestRegressor',
               ]
    estimators = all_estimators()

    n_samples, n_features = 20, 5
    rng = check_random_state(random_state)

    sample_weight = rng.randint(1, 4, (n_samples,))

    X_clf, y_clf = datasets.make_classification(
        n_samples=n_samples, n_features=n_features,
        random_state=random_state)
    X_reg, y_reg = datasets.make_regression(
        n_samples=n_samples, n_features=n_features,
        n_informative=2, random_state=random_state)

    def aug(data, sample_weight):
        # raise all samples to multiplicity of the corresponding sampleweight
        aug_data = []
        for samples, weight in zip(zip(*data), sample_weight):
            for _ in range(weight):
                aug_data.append(samples)
        aug_data = map(np.array, zip(*aug_data))
        return aug_data

    train, test = train_test_split(range(n_samples))

    for name, Estimator in estimators:
        print ("%s is being analysed" % name)
        if name in exclude:
            print ("%s is being excluded" % name)
            continue
        if not has_fit_parameter(Estimator, 'sample_weight'):
            continue
        if is_classifier(Estimator):
            X, y = X_clf, y_clf
        elif is_regressor(Estimator):
            X, y = X_reg, y_reg
        else:
            print ("%s is neither classifier nor regressor" % name)
            continue

        estimator_sw = Estimator()
        set_random_state(estimator_sw, random_state=random_state)
        estimator_sw.fit(X[train], y[train], sample_weight=sample_weight[train])
        X_aug_train, y_aug_train = aug((X[train], y[train]),
                                   sample_weight[train])
        assert_equal(X_aug_train.shape[0], np.sum(sample_weight[train]))

        estimator_aug = Estimator()
        set_random_state(estimator_aug, random_state=random_state)
        estimator_aug.fit(X_aug_train, y_aug_train)

        precision = 6
        # if estimator has `coef_` attribute, then compare the two
        if hasattr(estimator_sw, 'coef_'):
            yield (assert_array_almost_equal,
                   estimator_sw.coef_, estimator_aug.coef_, precision, name+' coef_ not equal')

        pred_sw = estimator_sw.predict(X[test])
        pred_aug = estimator_aug.predict(X[test])

        yield assert_array_almost_equal, pred_sw, pred_aug, precision, name+' prediction not equal'


def test_sample_weight_0(random_state=42):
    exclude = [
               'AdaBoostClassifier', 'AdaBoostRegressor',
               'BaggingClassifier', 'BaggingRegressor',
               'GradientBoostingClassifier', 'GradientBoostingRegressor',
               'ExtraTreeClassifier', 'ExtraTreeRegressor',
               'ExtraTreesClassifier', 'ExtraTreesRegressor',
               'DecisionTreeClassifier','DecisionTreeRegressor',
               'LogisticRegression', 'LogisticRegressionCV', 'LinearSVC',
               'MultinomialNB',  # Requires positive samples
               'CalibratedClassifierCV',  # This is a meta-estimator using LinearSVC
               'SGDClassifier',  # Doesn't work. Probably more data needed
               'SGDRegressor',  # Doesn't work. Probably more data needed
               'Perceptron',  # Uses SGD too. Doesn't work. Probably more data needed
               'RidgeClassifierCV', 'RidgeCV',
               'RandomForestClassifier', 'RandomForestRegressor',
               ]
    estimators = all_estimators()

    n_samples, n_features = 20, 5
    rng = check_random_state(random_state)

    sample_weight = (rng.permutation(n_samples) < (n_samples / 2.)).astype(np.float)

    X_clf, y_clf = datasets.make_classification(
        n_samples=n_samples, n_features=n_features,
        random_state=random_state)
    X_reg, y_reg = datasets.make_regression(
        n_samples=n_samples, n_features=n_features,
        n_informative=2, random_state=random_state)

    train, test = train_test_split(range(n_samples))

    for name, Estimator in estimators:
        print ("%s is being analysed" % name)
        if name in exclude:
            print ("%s is being excluded" % name)
            continue
        if not has_fit_parameter(Estimator, 'sample_weight'):
            continue
        if is_classifier(Estimator):
            X, y = X_clf, y_clf
        elif is_regressor(Estimator):
            X, y = X_reg, y_reg
        else:
            print ("%s is neither classifier nor regressor" % name)
            continue

        estimator_sw = Estimator()
        set_random_state(estimator_sw, random_state=random_state)
        estimator_sw.fit(X[train], y[train], sample_weight=sample_weight[train])
        X_aug_train, y_aug_train = X[train][sample_weight[train]==1], \
                                   y[train][sample_weight[train]==1]
        estimator_aug = Estimator()
        set_random_state(estimator_aug, random_state=random_state)
        estimator_aug.fit(X_aug_train, y_aug_train)

        precision = 6
        # if estimator has `coef_` attribute, then compare the two
        if hasattr(estimator_sw, 'coef_'):
            yield (assert_array_almost_equal,
                   estimator_sw.coef_, estimator_aug.coef_, precision, name+' coef_ not equal')

        pred_sw = estimator_sw.predict(X[test])
        pred_aug = estimator_aug.predict(X[test])

        yield assert_array_almost_equal, pred_sw, pred_aug, precision, name+' prediction not equal'


