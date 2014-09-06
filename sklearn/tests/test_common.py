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

from sklearn.externals.six import PY3
from sklearn.externals.six.moves import zip
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import ignore_warnings

import sklearn
from sklearn.base import (ClassifierMixin, RegressorMixin,
                          TransformerMixin, ClusterMixin)
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_classification

from sklearn.cross_validation import train_test_split
from sklearn.linear_model.base import LinearClassifierMixin
from sklearn.utils.estimator_checks import (
    check_parameters_default_constructible,
    check_regressors_classifiers_sparse_data,
    check_transformer,
    check_clustering,
    check_regressors_int,
    check_regressors_train,
    check_regressors_pickle,
    check_transformer_sparse_data,
    check_transformer_pickle,
    check_estimators_nan_inf,
    check_classifiers_one_label,
    check_classifiers_train,
    check_classifiers_classes,
    check_classifiers_input_shapes,
    check_classifiers_pickle,
    check_class_weight_classifiers,
    check_class_weight_auto_classifiers,
    check_class_weight_auto_linear_classifier,
    check_estimators_overwrite_params,
    check_cluster_overwrite_params,
    check_sparsify_binary_classifier,
    check_sparsify_multiclass_classifier,
    check_classifier_data_not_an_array,
    check_regressor_data_not_an_array,
    check_transformer_data_not_an_array,
    check_transformer_n_iter,
    check_non_transformer_estimators_n_iter,
    CROSS_DECOMPOSITION)


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


def test_estimators_sparse_data():
    # All estimators should either deal with sparse data or raise an
    # exception with type TypeError and an intelligible error message
    estimators = all_estimators()
    estimators = [(name, Estimator) for name, Estimator in estimators
                  if issubclass(Estimator, (ClassifierMixin, RegressorMixin))]
    for name, Estimator in estimators:
        yield check_regressors_classifiers_sparse_data, name, Estimator


def test_transformers():
    # test if transformers do something sensible on training set
    # also test all shapes / shape errors
    transformers = all_estimators(type_filter='transformer')
    for name, Transformer in transformers:
        # All transformers should either deal with sparse data or raise an
        # exception with type TypeError and an intelligible error message
        yield check_transformer_sparse_data, name, Transformer
        yield check_transformer_pickle, name, Transformer
        if name not in ['AdditiveChi2Sampler', 'Binarizer', 'Normalizer',
                        'PLSCanonical', 'PLSRegression', 'CCA', 'PLSSVD']:
            yield check_transformer_data_not_an_array, name, Transformer
        # these don't actually fit the data, so don't raise errors
        if name not in ['AdditiveChi2Sampler', 'Binarizer', 'Normalizer']:
            # basic tests
            yield check_transformer, name, Transformer


def test_estimators_nan_inf():
    # Test that all estimators check their input for NaN's and infs
    estimators = all_estimators()
    estimators = [(name, E) for name, E in estimators
                  if (issubclass(E, ClassifierMixin) or
                      issubclass(E, RegressorMixin) or
                      issubclass(E, TransformerMixin) or
                      issubclass(E, ClusterMixin))]
    for name, Estimator in estimators:
        if name not in CROSS_DECOMPOSITION + ['Imputer']:
            yield check_estimators_nan_inf, name, Estimator


def test_clustering():
    # test if clustering algorithms do something sensible
    # also test all shapes / shape errors
    clustering = all_estimators(type_filter='cluster')
    for name, Alg in clustering:
        # test whether any classifier overwrites his init parameters during fit
        yield check_cluster_overwrite_params, name, Alg
        if name not in ('WardAgglomeration', "FeatureAgglomeration"):
            # this is clustering on the features
            # let's not test that here.
            yield check_clustering, name, Alg


def test_classifiers():
    # test if classifiers can cope with non-consecutive classes
    classifiers = all_estimators(type_filter='classifier')
    for name, Classifier in classifiers:
        # test classfiers can handle non-array data
        yield check_classifier_data_not_an_array, name, Classifier
        # test classifiers trained on a single label always return this label
        yield check_classifiers_one_label, name, Classifier
        yield check_classifiers_classes, name, Classifier
        yield check_classifiers_pickle, name, Classifier
        # basic consistency testing
        yield check_classifiers_train, name, Classifier
        if (name not in ["MultinomialNB", "LabelPropagation", "LabelSpreading"]
            # TODO some complication with -1 label
                and name not in ["DecisionTreeClassifier", "ExtraTreeClassifier"]):
                # We don't raise a warning in these classifiers, as
                # the column y interface is used by the forests.

            # test if classifiers can cope with y.shape = (n_samples, 1)
            yield check_classifiers_input_shapes, name, Classifier


def test_regressors():
    regressors = all_estimators(type_filter='regressor')
    # TODO: test with intercept
    # TODO: test with multiple responses
    for name, Regressor in regressors:
        # basic testing
        yield check_regressors_train, name, Regressor
        yield check_regressor_data_not_an_array, name, Regressor
        # Test that estimators can be pickled, and once pickled
        # give the same answer as before.
        yield check_regressors_pickle, name, Regressor
        if name != 'CCA':
            # check that the regressor handles int input
            yield check_regressors_int, name, Regressor


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


def test_class_weight_classifiers():
    # test that class_weight works and that the semantics are consistent
    classifiers = all_estimators(type_filter='classifier')

    with warnings.catch_warnings(record=True):
        classifiers = [c for c in classifiers
                       if 'class_weight' in c[1]().get_params().keys()]

    for name, Classifier in classifiers:
        if name == "NuSVC":
            # the sparse version has a parameter that doesn't do anything
            continue
        if name.endswith("NB"):
            # NaiveBayes classifiers have a somewhat different interface.
            # FIXME SOON!
            continue
        yield check_class_weight_classifiers, name, Classifier


def test_class_weight_auto_classifiers():
    """Test that class_weight="auto" improves f1-score"""

    # This test is broken; its success depends on:
    # * a rare fortuitous RNG seed for make_classification; and
    # * the use of binary F1 over a seemingly arbitrary positive class for two
    #   datasets, and weighted average F1 for the third.
    # Its expectations need to be clarified and reimplemented.
    raise SkipTest('This test requires redefinition')

    classifiers = all_estimators(type_filter='classifier')

    with warnings.catch_warnings(record=True):
        classifiers = [c for c in classifiers
                       if 'class_weight' in c[1]().get_params().keys()]

    for n_classes, weights in zip([2, 3], [[.8, .2], [.8, .1, .1]]):
        # create unbalanced dataset
        X, y = make_classification(n_classes=n_classes, n_samples=200,
                                   n_features=10, weights=weights,
                                   random_state=0, n_informative=n_classes)
        X = StandardScaler().fit_transform(X)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                            random_state=0)
        for name, Classifier in classifiers:
            if (name != "NuSVC"
                # the sparse version has a parameter that doesn't do anything
                    and not name.startswith("RidgeClassifier")
                    # RidgeClassifier behaves unexpected
                    # FIXME!
                    and not name.endswith("NB")):
                # NaiveBayes classifiers have a somewhat different interface.
                # FIXME SOON!
                yield (check_class_weight_auto_classifiers, name, Classifier,
                       X_train, y_train, X_test, y_test, weights)


def test_class_weight_auto_linear_classifiers():
    classifiers = all_estimators(type_filter='classifier')

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
        yield check_class_weight_auto_linear_classifier, name, Classifier


def test_estimators_overwrite_params():
    # test whether any classifier overwrites his init parameters during fit
    for est_type in ["classifier", "regressor", "transformer"]:
        estimators = all_estimators(type_filter=est_type)
        for name, Estimator in estimators:
            if (name not in ['CCA', '_CCA', 'PLSCanonical', 'PLSRegression',
                             'PLSSVD', 'GaussianProcess']):
                # FIXME!
                # in particular GaussianProcess!
                yield check_estimators_overwrite_params, name, Estimator


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


def test_sparsify_estimators():
    #Test if predict with sparsified estimators works.
    #Tests regression, binary classification, and multi-class classification.
    estimators = all_estimators()

    # test regression and binary classification
    for name, Estimator in estimators:
        try:
            Estimator.sparsify
            yield check_sparsify_binary_classifier, name, Estimator
        except:
            pass

    # test multiclass classification
    classifiers = all_estimators(type_filter='classifier')
    for name, Classifier in classifiers:
        try:
            Classifier.sparsify
            yield check_sparsify_multiclass_classifier, name, Classifier
        except:
            pass


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
                elif name in CROSS_DECOMPOSITION or (
                    name in ['LinearSVC', 'LogisticRegression']
                    ):
                    continue

                else:
                    # Multitask models related to ENet cannot handle
                    # if y is mono-output.
                    yield (check_non_transformer_estimators_n_iter,
                           name, estimator, 'Multi' in name)


def test_transformer_n_iter():
    transformers = all_estimators(type_filter='transformer')
    for name, Estimator in transformers:
        estimator = Estimator()
        # Dependent on external solvers and hence accessing the iter
        # param is non-trivial.
        external_solver = ['Isomap', 'KernelPCA', 'LocallyLinearEmbedding',
                           'RandomizedLasso', 'LogisticRegressionCV']

        if hasattr(estimator, "max_iter") and name not in external_solver:
            yield check_transformer_n_iter, name, estimator
