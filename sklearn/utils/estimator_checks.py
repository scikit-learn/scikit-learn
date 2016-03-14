from __future__ import print_function

import types
import warnings
import sys
import traceback
import pickle
from copy import deepcopy

import numpy as np
from scipy import sparse
import struct

from sklearn.externals.six.moves import zip
from sklearn.externals.joblib import hash, Memory
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import META_ESTIMATORS
from sklearn.utils.testing import set_random_state
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_warns


from sklearn.base import (clone, ClassifierMixin, RegressorMixin,
                          TransformerMixin, ClusterMixin, BaseEstimator)
from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.random_projection import BaseRandomProjection
from sklearn.feature_selection import SelectKBest
from sklearn.svm.base import BaseLibSVM
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import NMF, ProjectedGradientNMF
from sklearn.exceptions import ConvergenceWarning
from sklearn.exceptions import DataConversionWarning
from sklearn.model_selection import train_test_split

from sklearn.utils import shuffle
from sklearn.utils.fixes import signature
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris, load_boston, make_blobs


BOSTON = None
CROSS_DECOMPOSITION = ['PLSCanonical', 'PLSRegression', 'CCA', 'PLSSVD']
MULTI_OUTPUT = ['CCA', 'DecisionTreeRegressor', 'ElasticNet',
                'ExtraTreeRegressor', 'ExtraTreesRegressor', 'GaussianProcess',
                'GaussianProcessRegressor',
                'KNeighborsRegressor', 'KernelRidge', 'Lars', 'Lasso',
                'LassoLars', 'LinearRegression', 'MultiTaskElasticNet',
                'MultiTaskElasticNetCV', 'MultiTaskLasso', 'MultiTaskLassoCV',
                'OrthogonalMatchingPursuit', 'PLSCanonical', 'PLSRegression',
                'RANSACRegressor', 'RadiusNeighborsRegressor',
                'RandomForestRegressor', 'Ridge', 'RidgeCV']

# Estimators with deprecated transform methods. Should be removed in 0.19 when
# _LearntSelectorMixin is removed.
DEPRECATED_TRANSFORM = [
    "RandomForestClassifier", "RandomForestRegressor", "ExtraTreesClassifier",
    "ExtraTreesRegressor", "DecisionTreeClassifier",
    "DecisionTreeRegressor", "ExtraTreeClassifier", "ExtraTreeRegressor",
    "LinearSVC", "SGDClassifier", "SGDRegressor", "Perceptron",
    "LogisticRegression", "LogisticRegressionCV",
    "GradientBoostingClassifier", "GradientBoostingRegressor"]


def _yield_non_meta_checks(name, Estimator):
    yield check_estimators_dtypes
    yield check_fit_score_takes_y
    yield check_dtype_object
    yield check_estimators_fit_returns_self

    # Check that all estimator yield informative messages when
    # trained on empty datasets
    yield check_estimators_empty_data_messages

    if name not in CROSS_DECOMPOSITION + ['SpectralEmbedding']:
        # SpectralEmbedding is non-deterministic,
        # see issue #4236
        # cross-decomposition's "transform" returns X and Y
        yield check_pipeline_consistency

    if name not in ['Imputer']:
        # Test that all estimators check their input for NaN's and infs
        yield check_estimators_nan_inf

    if name not in ['GaussianProcess']:
        # FIXME!
        # in particular GaussianProcess!
        yield check_estimators_overwrite_params
    if hasattr(Estimator, 'sparsify'):
        yield check_sparsify_coefficients

    yield check_estimator_sparse_data

    # Test that estimators can be pickled, and once pickled
    # give the same answer as before.
    yield check_estimators_pickle


def _yield_classifier_checks(name, Classifier):
    # test classifiers can handle non-array data
    yield check_classifier_data_not_an_array
    # test classifiers trained on a single label always return this label
    yield check_classifiers_one_label
    yield check_classifiers_classes
    yield check_estimators_partial_fit_n_features
    # basic consistency testing
    yield check_classifiers_train
    yield check_classifiers_regression_target
    if (name not in ["MultinomialNB", "LabelPropagation", "LabelSpreading"]
        # TODO some complication with -1 label
            and name not in ["DecisionTreeClassifier",
                             "ExtraTreeClassifier"]):
            # We don't raise a warning in these classifiers, as
            # the column y interface is used by the forests.

        yield check_supervised_y_2d
    # test if NotFittedError is raised
    yield check_estimators_unfitted
    if 'class_weight' in Classifier().get_params().keys():
        yield check_class_weight_classifiers

def check_supervised_y_no_nan(name, Estimator):
    # Checks that the Estimator targets are not NaN.

    rng = np.random.RandomState(888)
    X = rng.randn(10, 5)
    y = np.ones(10) * np.inf
    y = multioutput_estimator_convert_y_2d(name, y)

    errmsg = "Input contains NaN, infinity or a value too large for " \
             "dtype('float64')."
    try:
        Estimator().fit(X, y)
    except ValueError as e:
        if str(e) != errmsg:
            raise ValueError("Estimator {0} raised warning as expected, but "
                             "does not match expected error message" \
                             .format(name))
    else:
        raise ValueError("Estimator {0} should have raised error on fitting "
                         "array y with NaN value.".format(name))

def _yield_regressor_checks(name, Regressor):
    # TODO: test with intercept
    # TODO: test with multiple responses
    # basic testing
    yield check_regressors_train
    yield check_regressor_data_not_an_array
    yield check_estimators_partial_fit_n_features
    yield check_regressors_no_decision_function
    yield check_supervised_y_2d
    yield check_supervised_y_no_nan
    if name != 'CCA':
        # check that the regressor handles int input
        yield check_regressors_int
    if name != "GaussianProcessRegressor":
        # Test if NotFittedError is raised
        yield check_estimators_unfitted


def _yield_transformer_checks(name, Transformer):
    # All transformers should either deal with sparse data or raise an
    # exception with type TypeError and an intelligible error message
    if name not in ['AdditiveChi2Sampler', 'Binarizer', 'Normalizer',
                    'PLSCanonical', 'PLSRegression', 'CCA', 'PLSSVD']:
        yield check_transformer_data_not_an_array
    # these don't actually fit the data, so don't raise errors
    if name not in ['AdditiveChi2Sampler', 'Binarizer',
                    'FunctionTransformer', 'Normalizer']:
        # basic tests
        yield check_transformer_general
        yield check_transformers_unfitted


def _yield_clustering_checks(name, Clusterer):
    yield check_clusterer_compute_labels_predict
    if name not in ('WardAgglomeration', "FeatureAgglomeration"):
        # this is clustering on the features
        # let's not test that here.
        yield check_clustering
        yield check_estimators_partial_fit_n_features


def _yield_all_checks(name, Estimator):
    for check in _yield_non_meta_checks(name, Estimator):
        yield check
    if issubclass(Estimator, ClassifierMixin):
        for check in _yield_classifier_checks(name, Estimator):
            yield check
    if issubclass(Estimator, RegressorMixin):
        for check in _yield_regressor_checks(name, Estimator):
            yield check
    if issubclass(Estimator, TransformerMixin):
        if name not in DEPRECATED_TRANSFORM:
            for check in _yield_transformer_checks(name, Estimator):
                yield check
    if issubclass(Estimator, ClusterMixin):
        for check in _yield_clustering_checks(name, Estimator):
            yield check
    yield check_fit2d_predict1d
    yield check_fit2d_1sample
    yield check_fit2d_1feature
    yield check_fit1d_1feature
    yield check_fit1d_1sample


def check_estimator(Estimator):
    """Check if estimator adheres to sklearn conventions.

    This estimator will run an extensive test-suite for input validation,
    shapes, etc.
    Additional tests for classifiers, regressors, clustering or transformers
    will be run if the Estimator class inherits from the corresponding mixin
    from sklearn.base.

    Parameters
    ----------
    Estimator : class
        Class to check. Estimator is a class object (not an instance).

    """
    name = Estimator.__name__
    check_parameters_default_constructible(name, Estimator)
    for check in _yield_all_checks(name, Estimator):
        check(name, Estimator)


def _boston_subset(n_samples=200):
    global BOSTON
    if BOSTON is None:
        boston = load_boston()
        X, y = boston.data, boston.target
        X, y = shuffle(X, y, random_state=0)
        X, y = X[:n_samples], y[:n_samples]
        X = StandardScaler().fit_transform(X)
        BOSTON = X, y
    return BOSTON


def set_testing_parameters(estimator):
    # set parameters to speed up some estimators and
    # avoid deprecated behaviour
    params = estimator.get_params()
    if ("n_iter" in params
            and estimator.__class__.__name__ != "TSNE"):
        estimator.set_params(n_iter=5)
    if "max_iter" in params:
        warnings.simplefilter("ignore", ConvergenceWarning)
        if estimator.max_iter is not None:
            estimator.set_params(max_iter=min(5, estimator.max_iter))
        # LinearSVR
        if estimator.__class__.__name__ == 'LinearSVR':
            estimator.set_params(max_iter=20)
        # NMF
        if estimator.__class__.__name__ == 'NMF':
            estimator.set_params(max_iter=100)
        # MLP
        if estimator.__class__.__name__ in ['MLPClassifier', 'MLPRegressor']:
            estimator.set_params(max_iter=100)
    if "n_resampling" in params:
        # randomized lasso
        estimator.set_params(n_resampling=5)
    if "n_estimators" in params:
        # especially gradient boosting with default 100
        estimator.set_params(n_estimators=min(5, estimator.n_estimators))
    if "max_trials" in params:
        # RANSAC
        estimator.set_params(max_trials=10)
    if "n_init" in params:
        # K-Means
        estimator.set_params(n_init=2)
    if "decision_function_shape" in params:
        # SVC
        estimator.set_params(decision_function_shape='ovo')

    if estimator.__class__.__name__ == "SelectFdr":
        # be tolerant of noisy datasets (not actually speed)
        estimator.set_params(alpha=.5)

    if estimator.__class__.__name__ == "TheilSenRegressor":
        estimator.max_subpopulation = 100

    if isinstance(estimator, BaseRandomProjection):
        # Due to the jl lemma and often very few samples, the number
        # of components of the random matrix projection will be probably
        # greater than the number of features.
        # So we impose a smaller number (avoid "auto" mode)
        estimator.set_params(n_components=1)

    if isinstance(estimator, SelectKBest):
        # SelectKBest has a default of k=10
        # which is more feature than we have in most case.
        estimator.set_params(k=1)

    if isinstance(estimator, NMF):
        if not isinstance(estimator, ProjectedGradientNMF):
            estimator.set_params(solver='cd')


class NotAnArray(object):
    " An object that is convertable to an array"

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype=None):
        return self.data


def _is_32bit():
    """Detect if process is 32bit Python."""
    return struct.calcsize('P') * 8 == 32


def check_estimator_sparse_data(name, Estimator):
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X_csr = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    for sparse_format in ['csr', 'csc', 'dok', 'lil', 'coo', 'dia', 'bsr']:
        X = X_csr.asformat(sparse_format)
        # catch deprecation warnings
        with warnings.catch_warnings():
            if name in ['Scaler', 'StandardScaler']:
                estimator = Estimator(with_mean=False)
            else:
                estimator = Estimator()
        set_testing_parameters(estimator)
        # fit and predict
        try:
            estimator.fit(X, y)
            if hasattr(estimator, "predict"):
                pred = estimator.predict(X)
                assert_equal(pred.shape, (X.shape[0],))
            if hasattr(estimator, 'predict_proba'):
                probs = estimator.predict_proba(X)
                assert_equal(probs.shape, (X.shape[0], 4))
        except TypeError as e:
            if 'sparse' not in repr(e):
                print("Estimator %s doesn't seem to fail gracefully on "
                      "sparse data: error message state explicitly that "
                      "sparse input is not supported if this is not the case."
                      % name)
                raise
        except Exception:
            print("Estimator %s doesn't seem to fail gracefully on "
                  "sparse data: it should raise a TypeError if sparse input "
                  "is explicitly not supported." % name)
            raise


def check_dtype_object(name, Estimator):
    # check that estimators treat dtype object as numeric if possible
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10).astype(object)
    y = (X[:, 0] * 4).astype(np.int)
    y = multioutput_estimator_convert_y_2d(name, y)
    with warnings.catch_warnings():
        estimator = Estimator()
    set_testing_parameters(estimator)

    estimator.fit(X, y)
    if hasattr(estimator, "predict"):
        estimator.predict(X)

    if (hasattr(estimator, "transform") and
            name not in DEPRECATED_TRANSFORM):
        estimator.transform(X)

    try:
        estimator.fit(X, y.astype(object))
    except Exception as e:
        if "Unknown label type" not in str(e):
            raise

    X[0, 0] = {'foo': 'bar'}
    msg = "argument must be a string or a number"
    assert_raises_regex(TypeError, msg, estimator.fit, X, y)


@ignore_warnings
def check_fit2d_predict1d(name, Estimator):
    # check by fitting a 2d array and prediting with a 1d array
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    y = X[:, 0].astype(np.int)
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    estimator.fit(X, y)

    for method in ["predict", "transform", "decision_function",
                   "predict_proba"]:
        if hasattr(estimator, method):
            try:
                assert_warns(DeprecationWarning,
                             getattr(estimator, method), X[0])
            except ValueError:
                pass


@ignore_warnings
def check_fit2d_1sample(name, Estimator):
    # check by fitting a 2d array and prediting with a 1d array
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(1, 10))
    y = X[:, 0].astype(np.int)
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    try:
        estimator.fit(X, y)
    except ValueError:
        pass


@ignore_warnings
def check_fit2d_1feature(name, Estimator):
    # check by fitting a 2d array and prediting with a 1d array
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(10, 1))
    y = X[:, 0].astype(np.int)
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    try:
        estimator.fit(X, y)
    except ValueError:
        pass


@ignore_warnings
def check_fit1d_1feature(name, Estimator):
    # check fitting 1d array with 1 feature
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20))
    y = X.astype(np.int)
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)

    try:
        estimator.fit(X, y)
    except ValueError:
        pass


@ignore_warnings
def check_fit1d_1sample(name, Estimator):
    # check fitting 1d array with 1 feature
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20))
    y = np.array([1])
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)

    try:
        estimator.fit(X, y)
    except ValueError:
        pass


def check_transformer_general(name, Transformer):
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X = StandardScaler().fit_transform(X)
    X -= X.min()
    _check_transformer(name, Transformer, X, y)
    _check_transformer(name, Transformer, X.tolist(), y.tolist())


def check_transformer_data_not_an_array(name, Transformer):
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X = StandardScaler().fit_transform(X)
    # We need to make sure that we have non negative data, for things
    # like NMF
    X -= X.min() - .1
    this_X = NotAnArray(X)
    this_y = NotAnArray(np.asarray(y))
    _check_transformer(name, Transformer, this_X, this_y)


def check_transformers_unfitted(name, Transformer):
    X, y = _boston_subset()

    with warnings.catch_warnings(record=True):
        transformer = Transformer()

    assert_raises((AttributeError, ValueError), transformer.transform, X)


def _check_transformer(name, Transformer, X, y):
    if name in ('CCA', 'LocallyLinearEmbedding', 'KernelPCA') and _is_32bit():
        # Those transformers yield non-deterministic output when executed on
        # a 32bit Python. The same transformers are stable on 64bit Python.
        # FIXME: try to isolate a minimalistic reproduction case only depending
        # on numpy & scipy and/or maybe generate a test dataset that does not
        # cause such unstable behaviors.
        msg = name + ' is non deterministic on 32bit Python'
        raise SkipTest(msg)
    n_samples, n_features = np.asarray(X).shape
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        transformer = Transformer()
    set_random_state(transformer)
    set_testing_parameters(transformer)

    # fit

    if name in CROSS_DECOMPOSITION:
        y_ = np.c_[y, y]
        y_[::2, 1] *= 2
    else:
        y_ = y

    transformer.fit(X, y_)
    # fit_transform method should work on non fitted estimator
    transformer_clone = clone(transformer)
    X_pred = transformer_clone.fit_transform(X, y=y_)

    if isinstance(X_pred, tuple):
        for x_pred in X_pred:
            assert_equal(x_pred.shape[0], n_samples)
    else:
        # check for consistent n_samples
        assert_equal(X_pred.shape[0], n_samples)

    if hasattr(transformer, 'transform'):
        if name in CROSS_DECOMPOSITION:
            X_pred2 = transformer.transform(X, y_)
            X_pred3 = transformer.fit_transform(X, y=y_)
        else:
            X_pred2 = transformer.transform(X)
            X_pred3 = transformer.fit_transform(X, y=y_)
        if isinstance(X_pred, tuple) and isinstance(X_pred2, tuple):
            for x_pred, x_pred2, x_pred3 in zip(X_pred, X_pred2, X_pred3):
                assert_array_almost_equal(
                    x_pred, x_pred2, 2,
                    "fit_transform and transform outcomes not consistent in %s"
                    % Transformer)
                assert_array_almost_equal(
                    x_pred, x_pred3, 2,
                    "consecutive fit_transform outcomes not consistent in %s"
                    % Transformer)
        else:
            assert_array_almost_equal(
                X_pred, X_pred2, 2,
                "fit_transform and transform outcomes not consistent in %s"
                % Transformer)
            assert_array_almost_equal(
                X_pred, X_pred3, 2,
                "consecutive fit_transform outcomes not consistent in %s"
                % Transformer)
            assert_equal(len(X_pred2), n_samples)
            assert_equal(len(X_pred3), n_samples)

        # raises error on malformed input for transform
        if hasattr(X, 'T'):
            # If it's not an array, it does not have a 'T' property
            assert_raises(ValueError, transformer.transform, X.T)


@ignore_warnings
def check_pipeline_consistency(name, Estimator):
    if name in ('CCA', 'LocallyLinearEmbedding', 'KernelPCA') and _is_32bit():
        # Those transformers yield non-deterministic output when executed on
        # a 32bit Python. The same transformers are stable on 64bit Python.
        # FIXME: try to isolate a minimalistic reproduction case only depending
        # scipy and/or maybe generate a test dataset that does not
        # cause such unstable behaviors.
        msg = name + ' is non deterministic on 32bit Python'
        raise SkipTest(msg)

    # check that make_pipeline(est) gives same score as est
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X -= X.min()
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)
    set_random_state(estimator)
    pipeline = make_pipeline(estimator)
    estimator.fit(X, y)
    pipeline.fit(X, y)

    if name in DEPRECATED_TRANSFORM:
        funcs = ["score"]
    else:
        funcs = ["score", "fit_transform"]

    for func_name in funcs:
        func = getattr(estimator, func_name, None)
        if func is not None:
            func_pipeline = getattr(pipeline, func_name)
            result = func(X, y)
            result_pipe = func_pipeline(X, y)
            assert_array_almost_equal(result, result_pipe)


@ignore_warnings
def check_fit_score_takes_y(name, Estimator):
    # check that all estimators accept an optional y
    # in fit and score so they can be used in pipelines
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 3))
    y = np.arange(10) % 3
    y = multioutput_estimator_convert_y_2d(name, y)
    estimator = Estimator()
    set_testing_parameters(estimator)
    set_random_state(estimator)

    if name in DEPRECATED_TRANSFORM:
        funcs = ["fit", "score", "partial_fit", "fit_predict"]
    else:
        funcs = [
            "fit", "score", "partial_fit", "fit_predict", "fit_transform"]
    for func_name in funcs:
        func = getattr(estimator, func_name, None)
        if func is not None:
            func(X, y)
            args = [p.name for p in signature(func).parameters.values()]
            assert_true(args[1] in ["y", "Y"],
                        "Expected y or Y as second argument for method "
                        "%s of %s. Got arguments: %r."
                        % (func_name, Estimator.__name__, args))


@ignore_warnings
def check_estimators_dtypes(name, Estimator):
    rnd = np.random.RandomState(0)
    X_train_32 = 3 * rnd.uniform(size=(20, 5)).astype(np.float32)
    X_train_64 = X_train_32.astype(np.float64)
    X_train_int_64 = X_train_32.astype(np.int64)
    X_train_int_32 = X_train_32.astype(np.int32)
    y = X_train_int_64[:, 0]
    y = multioutput_estimator_convert_y_2d(name, y)

    if name in DEPRECATED_TRANSFORM:
        methods = ["predict", "decision_function", "predict_proba"]
    else:
        methods = [
            "predict", "transform", "decision_function", "predict_proba"]

    for X_train in [X_train_32, X_train_64, X_train_int_64, X_train_int_32]:
        with warnings.catch_warnings(record=True):
            estimator = Estimator()
        set_testing_parameters(estimator)
        set_random_state(estimator, 1)
        estimator.fit(X_train, y)

        for method in methods:
            if hasattr(estimator, method):
                getattr(estimator, method)(X_train)


def check_estimators_empty_data_messages(name, Estimator):
    e = Estimator()
    set_testing_parameters(e)
    set_random_state(e, 1)

    X_zero_samples = np.empty(0).reshape(0, 3)
    # The precise message can change depending on whether X or y is
    # validated first. Let us test the type of exception only:
    assert_raises(ValueError, e.fit, X_zero_samples, [])

    X_zero_features = np.empty(0).reshape(3, 0)
    # the following y should be accepted by both classifiers and regressors
    # and ignored by unsupervised models
    y = multioutput_estimator_convert_y_2d(name, np.array([1, 0, 1]))
    msg = "0 feature\(s\) \(shape=\(3, 0\)\) while a minimum of \d* is required."
    assert_raises_regex(ValueError, msg, e.fit, X_zero_features, y)


def check_estimators_nan_inf(name, Estimator):
    # Checks that Estimator X's do not contain NaN or inf.
    rnd = np.random.RandomState(0)
    X_train_finite = rnd.uniform(size=(10, 3))
    X_train_nan = rnd.uniform(size=(10, 3))
    X_train_nan[0, 0] = np.nan
    X_train_inf = rnd.uniform(size=(10, 3))
    X_train_inf[0, 0] = np.inf
    y = np.ones(10)
    y[:5] = 0
    y = multioutput_estimator_convert_y_2d(name, y)
    error_string_fit = "Estimator doesn't check for NaN and inf in fit."
    error_string_predict = ("Estimator doesn't check for NaN and inf in"
                            " predict.")
    error_string_transform = ("Estimator doesn't check for NaN and inf in"
                              " transform.")
    for X_train in [X_train_nan, X_train_inf]:
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            estimator = Estimator()
            set_testing_parameters(estimator)
            set_random_state(estimator, 1)
            # try to fit
            try:
                estimator.fit(X_train, y)
            except ValueError as e:
                if 'inf' not in repr(e) and 'NaN' not in repr(e):
                    print(error_string_fit, Estimator, e)
                    traceback.print_exc(file=sys.stdout)
                    raise e
            except Exception as exc:
                print(error_string_fit, Estimator, exc)
                traceback.print_exc(file=sys.stdout)
                raise exc
            else:
                raise AssertionError(error_string_fit, Estimator)
            # actually fit
            estimator.fit(X_train_finite, y)

            # predict
            if hasattr(estimator, "predict"):
                try:
                    estimator.predict(X_train)
                except ValueError as e:
                    if 'inf' not in repr(e) and 'NaN' not in repr(e):
                        print(error_string_predict, Estimator, e)
                        traceback.print_exc(file=sys.stdout)
                        raise e
                except Exception as exc:
                    print(error_string_predict, Estimator, exc)
                    traceback.print_exc(file=sys.stdout)
                else:
                    raise AssertionError(error_string_predict, Estimator)

            # transform
            if (hasattr(estimator, "transform") and
                    name not in DEPRECATED_TRANSFORM):
                try:
                    estimator.transform(X_train)
                except ValueError as e:
                    if 'inf' not in repr(e) and 'NaN' not in repr(e):
                        print(error_string_transform, Estimator, e)
                        traceback.print_exc(file=sys.stdout)
                        raise e
                except Exception as exc:
                    print(error_string_transform, Estimator, exc)
                    traceback.print_exc(file=sys.stdout)
                else:
                    raise AssertionError(error_string_transform, Estimator)


@ignore_warnings
def check_estimators_pickle(name, Estimator):
    """Test that we can pickle all estimators"""
    if name in DEPRECATED_TRANSFORM:
        check_methods = ["predict", "decision_function", "predict_proba"]
    else:
        check_methods = ["predict", "transform", "decision_function",
                         "predict_proba"]

    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)

    # some estimators can't do features less than 0
    X -= X.min()

    # some estimators only take multioutputs
    y = multioutput_estimator_convert_y_2d(name, y)

    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        estimator = Estimator()

    set_random_state(estimator)
    set_testing_parameters(estimator)
    estimator.fit(X, y)

    result = dict()
    for method in check_methods:
        if hasattr(estimator, method):
            result[method] = getattr(estimator, method)(X)

    # pickle and unpickle!
    pickled_estimator = pickle.dumps(estimator)
    unpickled_estimator = pickle.loads(pickled_estimator)

    for method in result:
        unpickled_result = getattr(unpickled_estimator, method)(X)
        assert_array_almost_equal(result[method], unpickled_result)


def check_estimators_partial_fit_n_features(name, Alg):
    # check if number of features changes between calls to partial_fit.
    if not hasattr(Alg, 'partial_fit'):
        return
    X, y = make_blobs(n_samples=50, random_state=1)
    X -= X.min()
    with warnings.catch_warnings(record=True):
        alg = Alg()
    if not hasattr(alg, 'partial_fit'):
        # check again as for mlp this depends on algorithm
        return

    set_testing_parameters(alg)
    try:
        if isinstance(alg, ClassifierMixin):
            classes = np.unique(y)
            alg.partial_fit(X, y, classes=classes)
        else:
            alg.partial_fit(X, y)
    except NotImplementedError:
        return

    assert_raises(ValueError, alg.partial_fit, X[:, :-1], y)


def check_clustering(name, Alg):
    X, y = make_blobs(n_samples=50, random_state=1)
    X, y = shuffle(X, y, random_state=7)
    X = StandardScaler().fit_transform(X)
    n_samples, n_features = X.shape
    # catch deprecation and neighbors warnings
    with warnings.catch_warnings(record=True):
        alg = Alg()
    set_testing_parameters(alg)
    if hasattr(alg, "n_clusters"):
        alg.set_params(n_clusters=3)
    set_random_state(alg)
    if name == 'AffinityPropagation':
        alg.set_params(preference=-100)
        alg.set_params(max_iter=100)

    # fit
    alg.fit(X)
    # with lists
    alg.fit(X.tolist())

    assert_equal(alg.labels_.shape, (n_samples,))
    pred = alg.labels_
    assert_greater(adjusted_rand_score(pred, y), 0.4)
    # fit another time with ``fit_predict`` and compare results
    if name is 'SpectralClustering':
        # there is no way to make Spectral clustering deterministic :(
        return
    set_random_state(alg)
    with warnings.catch_warnings(record=True):
        pred2 = alg.fit_predict(X)
    assert_array_equal(pred, pred2)


def check_clusterer_compute_labels_predict(name, Clusterer):
    """Check that predict is invariant of compute_labels"""
    X, y = make_blobs(n_samples=20, random_state=0)
    clusterer = Clusterer()

    if hasattr(clusterer, "compute_labels"):
        # MiniBatchKMeans
        if hasattr(clusterer, "random_state"):
            clusterer.set_params(random_state=0)

        X_pred1 = clusterer.fit(X).predict(X)
        clusterer.set_params(compute_labels=False)
        X_pred2 = clusterer.fit(X).predict(X)
        assert_array_equal(X_pred1, X_pred2)


def check_classifiers_one_label(name, Classifier):
    error_string_fit = "Classifier can't train when only one class is present."
    error_string_predict = ("Classifier can't predict when only one class is "
                            "present.")
    rnd = np.random.RandomState(0)
    X_train = rnd.uniform(size=(10, 3))
    X_test = rnd.uniform(size=(10, 3))
    y = np.ones(10)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
        set_testing_parameters(classifier)
        # try to fit
        try:
            classifier.fit(X_train, y)
        except ValueError as e:
            if 'class' not in repr(e):
                print(error_string_fit, Classifier, e)
                traceback.print_exc(file=sys.stdout)
                raise e
            else:
                return
        except Exception as exc:
            print(error_string_fit, Classifier, exc)
            traceback.print_exc(file=sys.stdout)
            raise exc
        # predict
        try:
            assert_array_equal(classifier.predict(X_test), y)
        except Exception as exc:
            print(error_string_predict, Classifier, exc)
            raise exc


@ignore_warnings  # Warnings are raised by decision function
def check_classifiers_train(name, Classifier):
    X_m, y_m = make_blobs(n_samples=300, random_state=0)
    X_m, y_m = shuffle(X_m, y_m, random_state=7)
    X_m = StandardScaler().fit_transform(X_m)
    # generate binary problem from multi-class one
    y_b = y_m[y_m != 2]
    X_b = X_m[y_m != 2]
    for (X, y) in [(X_m, y_m), (X_b, y_b)]:
        # catch deprecation warnings
        classes = np.unique(y)
        n_classes = len(classes)
        n_samples, n_features = X.shape
        with warnings.catch_warnings(record=True):
            classifier = Classifier()
        if name in ['BernoulliNB', 'MultinomialNB']:
            X -= X.min()
        set_testing_parameters(classifier)
        set_random_state(classifier)
        # raises error on malformed input for fit
        assert_raises(ValueError, classifier.fit, X, y[:-1])

        # fit
        classifier.fit(X, y)
        # with lists
        classifier.fit(X.tolist(), y.tolist())
        assert_true(hasattr(classifier, "classes_"))
        y_pred = classifier.predict(X)
        assert_equal(y_pred.shape, (n_samples,))
        # training set performance
        if name not in ['BernoulliNB', 'MultinomialNB']:
            assert_greater(accuracy_score(y, y_pred), 0.83)

        # raises error on malformed input for predict
        assert_raises(ValueError, classifier.predict, X.T)
        if hasattr(classifier, "decision_function"):
            try:
                # decision_function agrees with predict
                decision = classifier.decision_function(X)
                if n_classes is 2:
                    assert_equal(decision.shape, (n_samples,))
                    dec_pred = (decision.ravel() > 0).astype(np.int)
                    assert_array_equal(dec_pred, y_pred)
                if (n_classes is 3
                        and not isinstance(classifier, BaseLibSVM)):
                    # 1on1 of LibSVM works differently
                    assert_equal(decision.shape, (n_samples, n_classes))
                    assert_array_equal(np.argmax(decision, axis=1), y_pred)

                # raises error on malformed input
                assert_raises(ValueError,
                              classifier.decision_function, X.T)
                # raises error on malformed input for decision_function
                assert_raises(ValueError,
                              classifier.decision_function, X.T)
            except NotImplementedError:
                pass
        if hasattr(classifier, "predict_proba"):
            # predict_proba agrees with predict
            y_prob = classifier.predict_proba(X)
            assert_equal(y_prob.shape, (n_samples, n_classes))
            assert_array_equal(np.argmax(y_prob, axis=1), y_pred)
            # check that probas for all classes sum to one
            assert_array_almost_equal(np.sum(y_prob, axis=1),
                                      np.ones(n_samples))
            # raises error on malformed input
            assert_raises(ValueError, classifier.predict_proba, X.T)
            # raises error on malformed input for predict_proba
            assert_raises(ValueError, classifier.predict_proba, X.T)


def check_estimators_fit_returns_self(name, Estimator):
    """Check if self is returned when calling fit"""
    X, y = make_blobs(random_state=0, n_samples=9, n_features=4)
    y = multioutput_estimator_convert_y_2d(name, y)
    # some want non-negative input
    X -= X.min()

    estimator = Estimator()

    set_testing_parameters(estimator)
    set_random_state(estimator)

    assert_true(estimator.fit(X, y) is estimator)


@ignore_warnings
def check_estimators_unfitted(name, Estimator):
    """Check that predict raises an exception in an unfitted estimator.

    Unfitted estimators should raise either AttributeError or ValueError.
    The specific exception type NotFittedError inherits from both and can
    therefore be adequately raised for that purpose.
    """

    # Common test for Regressors as well as Classifiers
    X, y = _boston_subset()

    with warnings.catch_warnings(record=True):
        est = Estimator()

    msg = "fit"
    if hasattr(est, 'predict'):
        assert_raise_message((AttributeError, ValueError), msg,
                             est.predict, X)

    if hasattr(est, 'decision_function'):
        assert_raise_message((AttributeError, ValueError), msg,
                             est.decision_function, X)

    if hasattr(est, 'predict_proba'):
        assert_raise_message((AttributeError, ValueError), msg,
                             est.predict_proba, X)

    if hasattr(est, 'predict_log_proba'):
        assert_raise_message((AttributeError, ValueError), msg,
                             est.predict_log_proba, X)


def check_supervised_y_2d(name, Estimator):
    if "MultiTask" in name:
        # These only work on 2d, so this test makes no sense
        return
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 3))
    y = np.arange(10) % 3
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        estimator = Estimator()
    set_testing_parameters(estimator)
    set_random_state(estimator)
    # fit
    estimator.fit(X, y)
    y_pred = estimator.predict(X)

    set_random_state(estimator)
    # Check that when a 2D y is given, a DataConversionWarning is
    # raised
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", DataConversionWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        estimator.fit(X, y[:, np.newaxis])
    y_pred_2d = estimator.predict(X)
    msg = "expected 1 DataConversionWarning, got: %s" % (
        ", ".join([str(w_x) for w_x in w]))
    if name not in MULTI_OUTPUT:
        # check that we warned if we don't support multi-output
        assert_greater(len(w), 0, msg)
        assert_true("DataConversionWarning('A column-vector y"
                    " was passed when a 1d array was expected" in msg)
    assert_array_almost_equal(y_pred.ravel(), y_pred_2d.ravel())


def check_classifiers_classes(name, Classifier):
    X, y = make_blobs(n_samples=30, random_state=0, cluster_std=0.1)
    X, y = shuffle(X, y, random_state=7)
    X = StandardScaler().fit_transform(X)
    # We need to make sure that we have non negative data, for things
    # like NMF
    X -= X.min() - .1
    y_names = np.array(["one", "two", "three"])[y]

    for y_names in [y_names, y_names.astype('O')]:
        if name in ["LabelPropagation", "LabelSpreading"]:
            # TODO some complication with -1 label
            y_ = y
        else:
            y_ = y_names

        classes = np.unique(y_)
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            classifier = Classifier()
        if name == 'BernoulliNB':
            classifier.set_params(binarize=X.mean())
        set_testing_parameters(classifier)
        set_random_state(classifier)
        # fit
        classifier.fit(X, y_)

        y_pred = classifier.predict(X)
        # training set performance
        assert_array_equal(np.unique(y_), np.unique(y_pred))
        if np.any(classifier.classes_ != classes):
            print("Unexpected classes_ attribute for %r: "
                  "expected %s, got %s" %
                  (classifier, classes, classifier.classes_))


def check_regressors_int(name, Regressor):
    X, _ = _boston_subset()
    X = X[:50]
    rnd = np.random.RandomState(0)
    y = rnd.randint(3, size=X.shape[0])
    y = multioutput_estimator_convert_y_2d(name, y)
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        # separate estimators to control random seeds
        regressor_1 = Regressor()
        regressor_2 = Regressor()
    set_testing_parameters(regressor_1)
    set_testing_parameters(regressor_2)
    set_random_state(regressor_1)
    set_random_state(regressor_2)

    if name in CROSS_DECOMPOSITION:
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y

    # fit
    regressor_1.fit(X, y_)
    pred1 = regressor_1.predict(X)
    regressor_2.fit(X, y_.astype(np.float))
    pred2 = regressor_2.predict(X)
    assert_array_almost_equal(pred1, pred2, 2, name)


def check_regressors_train(name, Regressor):
    X, y = _boston_subset()
    y = StandardScaler().fit_transform(y.reshape(-1, 1))  # X is already scaled
    y = y.ravel()
    y = multioutput_estimator_convert_y_2d(name, y)
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        regressor = Regressor()
    set_testing_parameters(regressor)
    if not hasattr(regressor, 'alphas') and hasattr(regressor, 'alpha'):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01
    if name == 'PassiveAggressiveRegressor':
        regressor.C = 0.01

    # raises error on malformed input for fit
    assert_raises(ValueError, regressor.fit, X, y[:-1])
    # fit
    if name in CROSS_DECOMPOSITION:
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y
    set_random_state(regressor)
    regressor.fit(X, y_)
    regressor.fit(X.tolist(), y_.tolist())
    y_pred = regressor.predict(X)
    assert_equal(y_pred.shape, y_.shape)

    # TODO: find out why PLS and CCA fail. RANSAC is random
    # and furthermore assumes the presence of outliers, hence
    # skipped
    if name not in ('PLSCanonical', 'CCA', 'RANSACRegressor'):
        assert_greater(regressor.score(X, y_), 0.5)


@ignore_warnings
def check_regressors_no_decision_function(name, Regressor):
    # checks whether regressors have decision_function or predict_proba
    rng = np.random.RandomState(0)
    X = rng.normal(size=(10, 4))
    y = multioutput_estimator_convert_y_2d(name, X[:, 0])
    regressor = Regressor()

    set_testing_parameters(regressor)
    if hasattr(regressor, "n_components"):
        # FIXME CCA, PLS is not robust to rank 1 effects
        regressor.n_components = 1

    regressor.fit(X, y)
    funcs = ["decision_function", "predict_proba", "predict_log_proba"]
    for func_name in funcs:
        func = getattr(regressor, func_name, None)
        if func is None:
            # doesn't have function
            continue
        # has function. Should raise deprecation warning
        msg = func_name
        assert_warns_message(DeprecationWarning, msg, func, X)


def check_class_weight_classifiers(name, Classifier):
    if name == "NuSVC":
        # the sparse version has a parameter that doesn't do anything
        raise SkipTest
    if name.endswith("NB"):
        # NaiveBayes classifiers have a somewhat different interface.
        # FIXME SOON!
        raise SkipTest

    for n_centers in [2, 3]:
        # create a very noisy dataset
        X, y = make_blobs(centers=n_centers, random_state=0, cluster_std=20)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                            random_state=0)
        n_centers = len(np.unique(y_train))

        if n_centers == 2:
            class_weight = {0: 1000, 1: 0.0001}
        else:
            class_weight = {0: 1000, 1: 0.0001, 2: 0.0001}

        with warnings.catch_warnings(record=True):
            classifier = Classifier(class_weight=class_weight)
        if hasattr(classifier, "n_iter"):
            classifier.set_params(n_iter=100)
        if hasattr(classifier, "min_weight_fraction_leaf"):
            classifier.set_params(min_weight_fraction_leaf=0.01)

        set_random_state(classifier)
        classifier.fit(X_train, y_train)
        y_pred = classifier.predict(X_test)
        assert_greater(np.mean(y_pred == 0), 0.89)


def check_class_weight_balanced_classifiers(name, Classifier, X_train, y_train,
                                            X_test, y_test, weights):
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    if hasattr(classifier, "n_iter"):
        classifier.set_params(n_iter=100)

    set_random_state(classifier)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)

    classifier.set_params(class_weight='balanced')
    classifier.fit(X_train, y_train)
    y_pred_balanced = classifier.predict(X_test)
    assert_greater(f1_score(y_test, y_pred_balanced, average='weighted'),
                   f1_score(y_test, y_pred, average='weighted'))


def check_class_weight_balanced_linear_classifier(name, Classifier):
    """Test class weights with non-contiguous class labels."""
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                  [1.0, 1.0], [1.0, 0.0]])
    y = np.array([1, 1, 1, -1, -1])

    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    if hasattr(classifier, "n_iter"):
        # This is a very small dataset, default n_iter are likely to prevent
        # convergence
        classifier.set_params(n_iter=1000)
    set_random_state(classifier)

    # Let the model compute the class frequencies
    classifier.set_params(class_weight='balanced')
    coef_balanced = classifier.fit(X, y).coef_.copy()

    # Count each label occurrence to reweight manually
    n_samples = len(y)
    n_classes = float(len(np.unique(y)))

    class_weight = {1: n_samples / (np.sum(y == 1) * n_classes),
                    -1: n_samples / (np.sum(y == -1) * n_classes)}
    classifier.set_params(class_weight=class_weight)
    coef_manual = classifier.fit(X, y).coef_.copy()

    assert_array_almost_equal(coef_balanced, coef_manual)


def check_estimators_overwrite_params(name, Estimator):
    X, y = make_blobs(random_state=0, n_samples=9)
    y = multioutput_estimator_convert_y_2d(name, y)
    # some want non-negative input
    X -= X.min()
    with warnings.catch_warnings(record=True):
        # catch deprecation warnings
        estimator = Estimator()

    set_testing_parameters(estimator)
    set_random_state(estimator)

    # Make a physical copy of the original estimator parameters before fitting.
    params = estimator.get_params()
    original_params = deepcopy(params)

    # Fit the model
    estimator.fit(X, y)

    # Compare the state of the model parameters with the original parameters
    new_params = estimator.get_params()
    for param_name, original_value in original_params.items():
        new_value = new_params[param_name]

        # We should never change or mutate the internal state of input
        # parameters by default. To check this we use the joblib.hash function
        # that introspects recursively any subobjects to compute a checksum.
        # The only exception to this rule of immutable constructor parameters
        # is possible RandomState instance but in this check we explicitly
        # fixed the random_state params recursively to be integer seeds.
        assert_equal(hash(new_value), hash(original_value),
                     "Estimator %s should not change or mutate "
                     " the parameter %s from %s to %s during fit."
                     % (name, param_name, original_value, new_value))


def check_sparsify_coefficients(name, Estimator):
    X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1],
                  [-1, -2], [2, 2], [-2, -2]])
    y = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    est = Estimator()

    est.fit(X, y)
    pred_orig = est.predict(X)

    # test sparsify with dense inputs
    est.sparsify()
    assert_true(sparse.issparse(est.coef_))
    pred = est.predict(X)
    assert_array_equal(pred, pred_orig)

    # pickle and unpickle with sparse coef_
    est = pickle.loads(pickle.dumps(est))
    assert_true(sparse.issparse(est.coef_))
    pred = est.predict(X)
    assert_array_equal(pred, pred_orig)


def check_classifier_data_not_an_array(name, Estimator):
    X = np.array([[3, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 1]])
    y = [1, 1, 1, 2, 2, 2]
    y = multioutput_estimator_convert_y_2d(name, y)
    check_estimators_data_not_an_array(name, Estimator, X, y)


def check_regressor_data_not_an_array(name, Estimator):
    X, y = _boston_subset(n_samples=50)
    y = multioutput_estimator_convert_y_2d(name, y)
    check_estimators_data_not_an_array(name, Estimator, X, y)


def check_estimators_data_not_an_array(name, Estimator, X, y):

    if name in CROSS_DECOMPOSITION:
        raise SkipTest
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        # separate estimators to control random seeds
        estimator_1 = Estimator()
        estimator_2 = Estimator()
    set_testing_parameters(estimator_1)
    set_testing_parameters(estimator_2)
    set_random_state(estimator_1)
    set_random_state(estimator_2)

    y_ = NotAnArray(np.asarray(y))
    X_ = NotAnArray(np.asarray(X))

    # fit
    estimator_1.fit(X_, y_)
    pred1 = estimator_1.predict(X_)
    estimator_2.fit(X, y)
    pred2 = estimator_2.predict(X)
    assert_array_almost_equal(pred1, pred2, 2, name)


def check_parameters_default_constructible(name, Estimator):
    classifier = LinearDiscriminantAnalysis()
    # test default-constructibility
    # get rid of deprecation warnings
    with warnings.catch_warnings(record=True):
        if name in META_ESTIMATORS:
            estimator = Estimator(classifier)
        else:
            estimator = Estimator()
        # test cloning
        clone(estimator)
        # test __repr__
        repr(estimator)
        # test that set_params returns self
        assert_true(estimator.set_params() is estimator)

        # test if init does nothing but set parameters
        # this is important for grid_search etc.
        # We get the default parameters from init and then
        # compare these against the actual values of the attributes.

        # this comes from getattr. Gets rid of deprecation decorator.
        init = getattr(estimator.__init__, 'deprecated_original',
                       estimator.__init__)

        try:
            def param_filter(p):
                """Identify hyper parameters of an estimator"""
                return (p.name != 'self'
                        and p.kind != p.VAR_KEYWORD
                        and p.kind != p.VAR_POSITIONAL)

            init_params = [p for p in signature(init).parameters.values()
                           if param_filter(p)]
        except (TypeError, ValueError):
            # init is not a python function.
            # true for mixins
            return
        params = estimator.get_params()
        if name in META_ESTIMATORS:
            # they can need a non-default argument
            init_params = init_params[1:]

        for init_param in init_params:
            assert_not_equal(init_param.default, init_param.empty,
                             "parameter %s for %s has no default value"
                             % (init_param.name, type(estimator).__name__))
            assert_in(type(init_param.default),
                      [str, int, float, bool, tuple, type(None),
                       np.float64, types.FunctionType, Memory])
            if init_param.name not in params.keys():
                # deprecated parameter, not in get_params
                assert_true(init_param.default is None)
                continue

            param_value = params[init_param.name]
            if isinstance(param_value, np.ndarray):
                assert_array_equal(param_value, init_param.default)
            else:
                assert_equal(param_value, init_param.default)


def multioutput_estimator_convert_y_2d(name, y):
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if "MultiTask" in name:
        return np.reshape(y, (-1, 1))
    return y


def check_non_transformer_estimators_n_iter(name, estimator,
                                            multi_output=False):
    # Check if all iterative solvers, run for more than one iteration

    iris = load_iris()
    X, y_ = iris.data, iris.target

    if multi_output:
        y_ = np.reshape(y_, (-1, 1))

    set_random_state(estimator, 0)
    if name == 'AffinityPropagation':
        estimator.fit(X)
    else:
        estimator.fit(X, y_)

    # HuberRegressor depends on scipy.optimize.fmin_l_bfgs_b
    # which doesn't return a n_iter for old versions of SciPy.
    if not (name == 'HuberRegressor' and estimator.n_iter_ is None):
        assert_greater(estimator.n_iter_, 0)


def check_transformer_n_iter(name, estimator):
    if name in CROSS_DECOMPOSITION:
        # Check using default data
        X = [[0., 0., 1.], [1., 0., 0.], [2., 2., 2.], [2., 5., 4.]]
        y_ = [[0.1, -0.2], [0.9, 1.1], [0.1, -0.5], [0.3, -0.2]]

    else:
        X, y_ = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                           random_state=0, n_features=2, cluster_std=0.1)
        X -= X.min() - 0.1
    set_random_state(estimator, 0)
    estimator.fit(X, y_)

    # These return a n_iter per component.
    if name in CROSS_DECOMPOSITION:
        for iter_ in estimator.n_iter_:
            assert_greater(iter_, 1)
    else:
        assert_greater(estimator.n_iter_, 1)


def check_get_params_invariance(name, estimator):
    class T(BaseEstimator):
        """Mock classifier
        """

        def __init__(self):
            pass

        def fit(self, X, y):
            return self

    if name in ('FeatureUnion', 'Pipeline'):
        e = estimator([('clf', T())])

    elif name in ('GridSearchCV', 'RandomizedSearchCV', 'SelectFromModel'):
        return

    else:
        e = estimator()

    shallow_params = e.get_params(deep=False)
    deep_params = e.get_params(deep=True)

    assert_true(all(item in deep_params.items() for item in
                    shallow_params.items()))


def check_classifiers_regression_target(name, Estimator):
    # Check if classifier throws an exception when fed regression targets

    boston = load_boston()
    X, y = boston.data, boston.target
    e = Estimator()
    msg = 'Unknown label type: '
    assert_raises_regex(ValueError, msg, e.fit, X, y)
