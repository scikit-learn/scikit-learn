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
import traceback
import inspect
import pickle
import pkgutil
import struct

import numpy as np
from scipy import sparse

from sklearn.externals.six import PY3
from sklearn.externals.six.moves import zip
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import meta_estimators
from sklearn.utils.testing import set_random_state
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import check_skip_travis
from sklearn.utils.testing import ignore_warnings

import sklearn
from sklearn.base import (clone, ClassifierMixin, RegressorMixin,
                          TransformerMixin, ClusterMixin)
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import (load_iris, load_boston, make_blobs,
                              make_classification)
from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

from sklearn.lda import LDA
from sklearn.svm.base import BaseLibSVM

from sklearn.cross_validation import train_test_split
from sklearn.utils.validation import DataConversionWarning

dont_test = ['SparseCoder', 'EllipticEnvelope', 'DictVectorizer',
             'LabelBinarizer', 'LabelEncoder', 'MultiLabelBinarizer',
             'TfidfTransformer', 'IsotonicRegression', 'OneHotEncoder',
             'RandomTreesEmbedding', 'FeatureHasher', 'DummyClassifier',
             'DummyRegressor', 'TruncatedSVD', 'PolynomialFeatures']

class NotAnArray(object):
    " An object that is convertable to an array"

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype=None):
        return self.data


def multioutput_estimator_convert_y_2d(name, y):
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if name in (['MultiTaskElasticNetCV', 'MultiTaskLassoCV',
                 'MultiTaskLasso', 'MultiTaskElasticNet']):
        return y[:, np.newaxis]
    return y


def test_all_estimators():
    # Test that estimators are default-constructible, clonable
    # and have working repr.
    estimators = all_estimators(include_meta_estimators=True)

    # Meta sanity-check to make sure that the estimator introspection runs
    # properly
    assert_greater(len(estimators), 0)

    for name, Estimator in estimators:
        # some can just not be sensibly default constructed
        if name in dont_test:
            continue
        yield check_parameters_default_constructible, name, Estimator


def check_parameters_default_constructible(name, Estimator):
    classifier = LDA()
    # test default-constructibility
    # get rid of deprecation warnings
    with warnings.catch_warnings(record=True):
        if name in meta_estimators:
            estimator = Estimator(classifier)
        else:
            estimator = Estimator()
        # test cloning
        clone(estimator)
        # test __repr__
        repr(estimator)
        # test that set_params returns self
        assert_true(isinstance(estimator.set_params(), Estimator))

        # test if init does nothing but set parameters
        # this is important for grid_search etc.
        # We get the default parameters from init and then
        # compare these against the actual values of the attributes.

        # this comes from getattr. Gets rid of deprecation decorator.
        init = getattr(estimator.__init__, 'deprecated_original',
                       estimator.__init__)
        try:
            args, varargs, kws, defaults = inspect.getargspec(init)
        except TypeError:
            # init is not a python function.
            # true for mixins
            return
        params = estimator.get_params()
        if name in meta_estimators:
            # they need a non-default argument
            args = args[2:]
        else:
            args = args[1:]
        if args:
            # non-empty list
            assert_equal(len(args), len(defaults))
        else:
            return
        for arg, default in zip(args, defaults):
            if arg not in params.keys():
                # deprecated parameter, not in get_params
                assert_true(default is None)
                continue

            if isinstance(params[arg], np.ndarray):
                assert_array_equal(params[arg], default)
            else:
                assert_equal(params[arg], default)


def test_all_estimator_no_base_class():
    for name, Estimator in all_estimators():
        msg = ("Base estimators such as {0} should not be included"
               " in all_estimators").format(name)
        assert_false(name.lower().startswith('base'), msg=msg)


def test_estimators_sparse_data():
    # All estimators should either deal with sparse data or raise an
    # exception with type TypeError and an intelligible error message
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    estimators = all_estimators()
    estimators = [(name, Estimator) for name, Estimator in estimators
                  if issubclass(Estimator, (ClassifierMixin, RegressorMixin))]
    for name, Estimator in estimators:
        if name in dont_test:
            continue
        yield check_regressors_classifiers_sparse_data, name, Estimator, X, y


def check_regressors_classifiers_sparse_data(name, Estimator, X, y):
    # catch deprecation warnings
    with warnings.catch_warnings():
        estimator = Estimator()
    # fit and predict
    try:
        estimator.fit(X, y)
        estimator.predict(X)
        if hasattr(estimator, 'predict_proba'):
            estimator.predict_proba(X)
    except TypeError as e:
        if not 'sparse' in repr(e):
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


def test_transformers():
    # test if transformers do something sensible on training set
    # also test all shapes / shape errors
    transformers = all_estimators(type_filter='transformer')
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X = StandardScaler().fit_transform(X)
    X -= X.min()

    for name, Transformer in transformers:
        if name in dont_test:
            continue
        # these don't actually fit the data:
        if name in ['AdditiveChi2Sampler', 'Binarizer', 'Normalizer']:
            continue
        yield check_transformer, name, Transformer, X, y


def _is_32bit():
    """Detect if process is 32bit Python."""
    return struct.calcsize('P') * 8 == 32


def check_transformer(name, Transformer, X, y):
    if name in ('CCA', 'LocallyLinearEmbedding', 'KernelPCA') and _is_32bit():
        # Those transformers yield non-deterministic output when executed on
        # a 32bit Python. The same transformers are stable on 64bit Python.
        # FIXME: try to isolate a minimalistic reproduction case only depending
        # on numpy & scipy and/or maybe generate a test dataset that does not
        # cause such unstable behaviors.
        msg = name + ' is non deterministic on 32bit Python'
        raise SkipTest(msg)

    n_samples, n_features = X.shape
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        transformer = Transformer()
    set_random_state(transformer)
    if hasattr(transformer, 'compute_importances'):
        transformer.compute_importances = True

    if name == 'SelectKBest':
        # SelectKBest has a default of k=10
        # which is more feature than we have.
        transformer.k = 1
    elif name in ['GaussianRandomProjection',
                  'SparseRandomProjection']:
        # Due to the jl lemma and very few samples, the number
        # of components of the random matrix projection will be greater
        # than the number of features.
        # So we impose a smaller number (avoid "auto" mode)
        transformer.n_components = 1
    elif name == "MiniBatchDictionaryLearning":
        transformer.set_params(n_iter=5)    # default = 1000

    elif name == "KernelPCA":
        transformer.remove_zero_eig = False

    # fit

    if name in ('PLSCanonical', 'PLSRegression', 'CCA', 'PLSSVD'):
        y_ = np.c_[y, y]
        y_[::2, 1] *= 2
    else:
        y_ = y

    transformer.fit(X, y_)
    X_pred = transformer.fit_transform(X, y=y_)
    if isinstance(X_pred, tuple):
        for x_pred in X_pred:
            assert_equal(x_pred.shape[0], n_samples)
    else:
        assert_equal(X_pred.shape[0], n_samples)

    if hasattr(transformer, 'transform'):
        if name in ('PLSCanonical', 'PLSRegression', 'CCA',
                    'PLSSVD'):
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

        # raises error on malformed input for transform
        assert_raises(ValueError, transformer.transform, X.T)


def test_transformers_sparse_data():
    # All transformers should either deal with sparse data or raise an
    # exception with type TypeError and an intelligible error message
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    estimators = all_estimators(type_filter='transformer')
    for name, Transformer in estimators:
        if name in dont_test:
            continue
        yield check_transformer_sparse_data, name, Transformer, X, y


def check_transformer_sparse_data(name, Transformer, X, y):
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        if name in ['Scaler', 'StandardScaler']:
            transformer = Transformer(with_mean=False)
        elif name in ['GaussianRandomProjection',
                      'SparseRandomProjection']:
            # Due to the jl lemma and very few samples, the number
            # of components of the random matrix projection will be greater
            # than the number of features.
            # So we impose a smaller number (avoid "auto" mode)
            transformer = Transformer(n_components=np.int(X.shape[1] / 4))
        else:
            transformer = Transformer()
    # fit
    try:
        transformer.fit(X, y)
    except TypeError as e:
        if not 'sparse' in repr(e):
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


def test_estimators_nan_inf():
    # Test that all estimators check their input for NaN's and infs
    rnd = np.random.RandomState(0)
    X_train_finite = rnd.uniform(size=(10, 3))
    X_train_nan = rnd.uniform(size=(10, 3))
    X_train_nan[0, 0] = np.nan
    X_train_inf = rnd.uniform(size=(10, 3))
    X_train_inf[0, 0] = np.inf
    y = np.ones(10)
    y[:5] = 0
    estimators = all_estimators()
    estimators = [(name, E) for name, E in estimators
                  if (issubclass(E, ClassifierMixin) or
                      issubclass(E, RegressorMixin) or
                      issubclass(E, TransformerMixin) or
                      issubclass(E, ClusterMixin))]
    for X_train in [X_train_nan, X_train_inf]:
        for name, Estimator in estimators:
            if name in dont_test:
                continue
            if name in ('PLSCanonical', 'PLSRegression', 'CCA',
                        'PLSSVD', 'Imputer'):  # Imputer accepts nan
                continue
            yield (check_estimators_nan_inf, name, Estimator, X_train,
                   X_train_finite,
                   multioutput_estimator_convert_y_2d(name, y))


def check_estimators_nan_inf(name, Estimator, X_train, X_train_finite, y):
    error_string_fit = "Estimator doesn't check for NaN and inf in fit."
    error_string_predict = ("Estimator doesn't check for NaN and inf in"
                            " predict.")
    error_string_transform = ("Estimator doesn't check for NaN and inf in"
                              " transform.")
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        estimator = Estimator()
        if name in ['GaussianRandomProjection',
                    'SparseRandomProjection']:
            # Due to the jl lemma and very few samples, the number
            # of components of the random matrix projection will be
            # greater
            # than the number of features.
            # So we impose a smaller number (avoid "auto" mode)
            estimator = Estimator(n_components=1)
        elif name == "SelectKBest":
            estimator = Estimator(k=1)

        set_random_state(estimator, 1)
        # try to fit
        try:
            if issubclass(Estimator, ClusterMixin):
                estimator.fit(X_train)
            else:
                estimator.fit(X_train, y)
        except ValueError as e:
            if not 'inf' in repr(e) and not 'NaN' in repr(e):
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
        if issubclass(Estimator, ClusterMixin):
            # All estimators except clustering algorithm
            # support fitting with (optional) y
            estimator.fit(X_train_finite)
        else:
            estimator.fit(X_train_finite, y)

        # predict
        if hasattr(estimator, "predict"):
            try:
                estimator.predict(X_train)
            except ValueError as e:
                if not 'inf' in repr(e) and not 'NaN' in repr(e):
                    print(error_string_predict, Estimator, e)
                    traceback.print_exc(file=sys.stdout)
                    raise e
            except Exception as exc:
                print(error_string_predict, Estimator, exc)
                traceback.print_exc(file=sys.stdout)
            else:
                raise AssertionError(error_string_predict, Estimator)

        # transform
        if hasattr(estimator, "transform"):
            try:
                estimator.transform(X_train)
            except ValueError as e:
                if not 'inf' in repr(e) and not 'NaN' in repr(e):
                    print(error_string_transform, Estimator, e)
                    traceback.print_exc(file=sys.stdout)
                    raise e
            except Exception as exc:
                print(error_string_transform, Estimator, exc)
                traceback.print_exc(file=sys.stdout)
            else:
                raise AssertionError(error_string_transform, Estimator)


def test_transformers_pickle():
    # test if transformers do something sensible on training set
    # also test all shapes / shape errors
    transformers = all_estimators(type_filter='transformer')
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    n_samples, n_features = X.shape
    X = StandardScaler().fit_transform(X)
    X -= X.min()

    for name, Transformer in transformers:
        if name in dont_test:
            continue
        yield check_transformer_pickle, name, Transformer, X, y


def check_transformer_pickle(name, Transformer, X, y):
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        transformer = Transformer()
    if not hasattr(transformer, 'transform'):
        return
    set_random_state(transformer)
    if hasattr(transformer, 'compute_importances'):
        transformer.compute_importances = True

    if name == "SelectKBest":
        # SelectKBest has a default of k=10
        # which is more feature than we have.
        transformer.k = 1
    elif name in ['GaussianRandomProjection', 'SparseRandomProjection']:
        # Due to the jl lemma and very few samples, the number
        # of components of the random matrix projection will be greater
        # than the number of features.
        # So we impose a smaller number (avoid "auto" mode)
        transformer.n_components = 1

    # fit
    if name in ('PLSCanonical', 'PLSRegression', 'CCA',
                'PLSSVD'):
        random_state = np.random.RandomState(seed=12345)
        y_ = np.vstack([y, 2 * y + random_state.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y

    transformer.fit(X, y_)
    X_pred = transformer.fit(X, y_).transform(X)
    pickled_transformer = pickle.dumps(transformer)
    unpickled_transformer = pickle.loads(pickled_transformer)
    pickled_X_pred = unpickled_transformer.transform(X)

    assert_array_almost_equal(pickled_X_pred, X_pred)


def test_classifiers_one_label():
    # test classifiers trained on a single label always return this label
    # or raise an sensible error message
    rnd = np.random.RandomState(0)
    X_train = rnd.uniform(size=(10, 3))
    X_test = rnd.uniform(size=(10, 3))
    y = np.ones(10)
    classifiers = all_estimators(type_filter='classifier')
    for name, Classifier in classifiers:
        if name in dont_test:
            continue
        yield check_classifiers_one_label, name, Classifier, X_train, X_test, y


def check_classifiers_one_label(name, Classifier, X_train, X_test, y):
    error_string_fit = "Classifier can't train when only one class is present."
    error_string_predict = ("Classifier can't predict when only one class is "
                            "present.")
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
        # try to fit
        try:
            classifier.fit(X_train, y)
        except ValueError as e:
            if not 'class' in repr(e):
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


def test_clustering():
    # test if clustering algorithms do something sensible
    # also test all shapes / shape errors
    clustering = all_estimators(type_filter='cluster')
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=7)
    X = StandardScaler().fit_transform(X)
    for name, Alg in clustering:
        if name in ('WardAgglomeration', "FeatureAgglomeration"):
            # this is clustering on the features
            # let's not test that here.
            continue
        yield check_clustering, name, Alg, X, y


def check_clustering(name, Alg, X, y):
    n_samples, n_features = X.shape
    # catch deprecation and neighbors warnings
    with warnings.catch_warnings(record=True):
        alg = Alg()
        if hasattr(alg, "n_clusters"):
            alg.set_params(n_clusters=3)
        set_random_state(alg)
        if name == 'AffinityPropagation':
            alg.set_params(preference=-100)
        # fit
        alg.fit(X)

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


def test_classifiers_train():
    # test if classifiers do something sensible on training set
    # also test all shapes / shape errors
    classifiers = all_estimators(type_filter='classifier')
    X_m, y_m = make_blobs(random_state=0)
    X_m, y_m = shuffle(X_m, y_m, random_state=7)
    X_m = StandardScaler().fit_transform(X_m)
    # generate binary problem from multi-class one
    y_b = y_m[y_m != 2]
    X_b = X_m[y_m != 2]
    for (X, y) in [(X_m, y_m), (X_b, y_b)]:
        # do it once with binary, once with multiclass
        for name, Classifier in classifiers:
            if name in dont_test:
                continue
            if name in ['MultinomialNB', 'BernoulliNB']:
                # TODO also test these!
                continue
            yield check_classifiers_train, name, Classifier, X, y


def check_classifiers_train(name, Classifier, X, y):
    # catch deprecation warnings
    classes = np.unique(y)
    n_classes = len(classes)
    n_samples, n_features = X.shape
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    # raises error on malformed input for fit
    assert_raises(ValueError, classifier.fit, X, y[:-1])

    # fit
    classifier.fit(X, y)
    assert_true(hasattr(classifier, "classes_"))
    y_pred = classifier.predict(X)
    assert_equal(y_pred.shape, (n_samples,))
    # training set performance
    assert_greater(accuracy_score(y, y_pred), 0.85)

    # raises error on malformed input for predict
    assert_raises(ValueError, classifier.predict, X.T)
    if hasattr(classifier, "decision_function"):
        try:
            # decision_function agrees with predict:
            decision = classifier.decision_function(X)
            if n_classes is 2:
                assert_equal(decision.ravel().shape, (n_samples,))
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
        # predict_proba agrees with predict:
        y_prob = classifier.predict_proba(X)
        assert_equal(y_prob.shape, (n_samples, n_classes))
        assert_array_equal(np.argmax(y_prob, axis=1), y_pred)
        # check that probas for all classes sum to one
        assert_array_almost_equal(np.sum(y_prob, axis=1), np.ones(n_samples))
        # raises error on malformed input
        assert_raises(ValueError, classifier.predict_proba, X.T)
        # raises error on malformed input for predict_proba
        assert_raises(ValueError, classifier.predict_proba, X.T)


def test_classifiers_classes():
    # test if classifiers can cope with non-consecutive classes
    classifiers = all_estimators(type_filter='classifier')
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=1)
    X = StandardScaler().fit_transform(X)
    y_names = iris.target_names[y]
    for y_names in [y_names, y_names.astype('O')]:
        for name, Classifier in classifiers:
            if name in dont_test:
                continue
            if name in ['MultinomialNB', 'BernoulliNB']:
                # TODO also test these!
                continue
            yield check_classifiers_classes, name, Classifier, X, y, y_names


def check_classifiers_classes(name, Classifier, X, y, y_names):
    if name in ["LabelPropagation", "LabelSpreading"]:
        # TODO some complication with -1 label
        y_ = y
    else:
        y_ = y_names

    classes = np.unique(y_)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    # fit
    try:
        classifier.fit(X, y_)
    except Exception as e:
        print(e)

    y_pred = classifier.predict(X)
    # training set performance
    assert_array_equal(np.unique(y_), np.unique(y_pred))
    accuracy = accuracy_score(y_, y_pred)
    assert_greater(accuracy, 0.78,
                   "accuracy %f of %s not greater than 0.78"
                   % (accuracy, name))
    #assert_array_equal(
        #clf.classes_, classes,
        #"Unexpected classes_ attribute for %r" % clf)
    if np.any(classifier.classes_ != classes):
        print("Unexpected classes_ attribute for %r: "
              "expected %s, got %s" %
              (classifier, classes, classifier.classes_))


def test_classifiers_input_shapes():
    # test if classifiers can cope with y.shape = (n_samples, 1)
    classifiers = all_estimators(type_filter='classifier')
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=1)
    X = StandardScaler().fit_transform(X)
    for name, Classifier in classifiers:
        if name in dont_test:
            continue
        if name in ["MultinomialNB", "LabelPropagation", "LabelSpreading"]:
            # TODO some complication with -1 label
            continue
        if name in ["DecisionTreeClassifier", "ExtraTreeClassifier"]:
            # We don't raise a warning in these classifiers, as
            # the column y interface is used by the forests.
            continue
        yield check_classifiers_input_shapes, name, Classifier, X, y


def check_classifiers_input_shapes(name, Classifier, X, y):
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    set_random_state(classifier)
    # fit
    classifier.fit(X, y)
    y_pred = classifier.predict(X)

    set_random_state(classifier)
    # Check that when a 2D y is given, a DataConversionWarning is
    # raised
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", DataConversionWarning)
        classifier.fit(X, y[:, np.newaxis])
    assert_equal(len(w), 1)
    assert_array_equal(y_pred, classifier.predict(X))


def test_classifiers_pickle():
    # test if classifiers do something sensible on training set
    # also test all shapes / shape errors
    classifiers = all_estimators(type_filter='classifier')
    X_m, y_m = make_blobs(random_state=0)
    X_m, y_m = shuffle(X_m, y_m, random_state=7)
    X_m = StandardScaler().fit_transform(X_m)
    # generate binary problem from multi-class one
    y_b = y_m[y_m != 2]
    X_b = X_m[y_m != 2]
    for (X, y) in [(X_m, y_m), (X_b, y_b)]:
        # do it once with binary, once with multiclass
        n_samples, n_features = X.shape
        for name, Classifier in classifiers:
            if name in dont_test:
                continue
            if name in ['MultinomialNB', 'BernoulliNB']:
                # TODO also test these!
                continue
            yield check_classifiers_pickle, name, Classifier, X, y


def check_classifiers_pickle(name, Classifier, X, y):
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    # raises error on malformed input for fit
    assert_raises(ValueError, classifier.fit, X, y[:-1])

    # fit
    classifier.fit(X, y)
    y_pred = classifier.predict(X)
    pickled_classifier = pickle.dumps(classifier)
    unpickled_classifier = pickle.loads(pickled_classifier)
    pickled_y_pred = unpickled_classifier.predict(X)

    assert_array_almost_equal(pickled_y_pred, y_pred)


BOSTON = None


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


def test_regressors_int():
    # test if regressors can cope with integer labels (by converting them to
    # float)
    regressors = all_estimators(type_filter='regressor')
    X, _ = _boston_subset()
    X = X[:50]
    rnd = np.random.RandomState(0)
    y = rnd.randint(3, size=X.shape[0])
    for name, Regressor in regressors:
        if name in dont_test or name in ('CCA'):
            continue
        elif name in ('OrthogonalMatchingPursuitCV'):
            # FIXME: This test is unstable on Travis, see issue #3190.
            check_skip_travis()
        yield (check_regressors_int, name, Regressor, X,
               multioutput_estimator_convert_y_2d(name, y))


def check_regressors_int(name, Regressor, X, y):
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        # separate estimators to control random seeds
        regressor_1 = Regressor()
        regressor_2 = Regressor()
    set_random_state(regressor_1)
    set_random_state(regressor_2)

    if name in ('_PLS', 'PLSCanonical', 'PLSRegression'):
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


def test_regressors_train():
    regressors = all_estimators(type_filter='regressor')
    # TODO: test with intercept
    # TODO: test with multiple responses
    X, y = _boston_subset()
    y = StandardScaler().fit_transform(y)   # X is already scaled
    for name, Regressor in regressors:
        if name in dont_test:
            continue
        yield (check_regressors_train, name, Regressor, X,
               multioutput_estimator_convert_y_2d(name, y))


def check_regressors_train(name, Regressor, X, y):
    if name == 'OrthogonalMatchingPursuitCV':
        # FIXME: This test is unstable on Travis, see issue #3190.
        check_skip_travis()
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        regressor = Regressor()
    if not hasattr(regressor, 'alphas') and hasattr(regressor, 'alpha'):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01

    # raises error on malformed input for fit
    assert_raises(ValueError, regressor.fit, X, y[:-1])
    # fit
    if name in ('PLSCanonical', 'PLSRegression', 'CCA'):
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y
    set_random_state(regressor)
    regressor.fit(X, y_)
    regressor.predict(X)

      # TODO: find out why PLS and CCA fail. RANSAC is random
      # and furthermore assumes the presence of outliers, hence
      # skipped
    if name not in ('PLSCanonical', 'CCA', 'RANSACRegressor'):
        assert_greater(regressor.score(X, y_), 0.5)


def test_regressor_pickle():
    # Test that estimators can be pickled, and once pickled
    # give the same answer as before.
    regressors = all_estimators(type_filter='regressor')
    X, y = _boston_subset()
    # TODO: test with intercept
    # TODO: test with multiple responses
    y = StandardScaler().fit_transform(y)   # X is already scaled
    for name, Regressor in regressors:
        if name in dont_test:
            continue
        yield (check_regressors_pickle, name, Regressor, X,
               multioutput_estimator_convert_y_2d(name, y))


def check_regressors_pickle(name, Regressor, X, y):
    if name == 'OrthogonalMatchingPursuitCV':
        # FIXME: This test is unstable on Travis, see issue #3190.
        check_skip_travis()
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        regressor = Regressor()
    if not hasattr(regressor, 'alphas') and hasattr(regressor, 'alpha'):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01

    if name in ('PLSCanonical', 'PLSRegression', 'CCA'):
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y
    regressor.fit(X, y_)
    y_pred = regressor.predict(X)
    # store old predictions
    pickled_regressor = pickle.dumps(regressor)
    unpickled_regressor = pickle.loads(pickled_regressor)
    pickled_y_pred = unpickled_regressor.predict(X)
    assert_array_almost_equal(pickled_y_pred, y_pred)


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

    for n_centers in [2, 3]:
        # create a very noisy dataset
        X, y = make_blobs(centers=n_centers, random_state=0, cluster_std=20)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                            random_state=0)
        for name, Classifier in classifiers:
            if name == "NuSVC":
                # the sparse version has a parameter that doesn't do anything
                continue
            if name.endswith("NB"):
                # NaiveBayes classifiers have a somewhat different interface.
                # FIXME SOON!
                continue
            yield (check_class_weight_classifiers, name, Classifier, X_train,
                   y_train, X_test, y_test)


def check_class_weight_classifiers(name, Classifier, X_train, y_train, X_test,
                                   y_test):
    n_centers = len(np.unique(y_train))

    if n_centers == 2:
        class_weight = {0: 1000, 1: 0.0001}
    else:
        class_weight = {0: 1000, 1: 0.0001, 2: 0.0001}

    with warnings.catch_warnings(record=True):
        classifier = Classifier(class_weight=class_weight)
    if hasattr(classifier, "n_iter"):
        classifier.set_params(n_iter=100)

    set_random_state(classifier)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    assert_greater(np.mean(y_pred == 0), 0.9)


def test_class_weight_auto_classifies():
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
            if name == "NuSVC":
                # the sparse version has a parameter that doesn't do anything
                continue

            if name.startswith("RidgeClassifier"):
                # RidgeClassifier behaves unexpected
                # FIXME!
                continue

            if name.endswith("NB"):
                # NaiveBayes classifiers have a somewhat different interface.
                # FIXME SOON!
                continue
            yield (check_class_weight_auto_classifiers, name, Classifier,
                   X_train, y_train, X_test, y_test, weights)


def check_class_weight_auto_classifiers(name, Classifier, X_train, y_train,
                                        X_test, y_test, weights):
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    if hasattr(classifier, "n_iter"):
        classifier.set_params(n_iter=100)

    set_random_state(classifier)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)

    classifier.set_params(class_weight='auto')
    classifier.fit(X_train, y_train)
    y_pred_auto = classifier.predict(X_test)
    assert_greater(f1_score(y_test, y_pred_auto),
                   f1_score(y_test, y_pred))


def test_estimators_overwrite_params():
    # test whether any classifier overwrites his init parameters during fit
    for est_type in ["classifier", "regressor", "transformer"]:
        estimators = all_estimators(type_filter=est_type)
        X, y = make_blobs(random_state=0, n_samples=9)
        # some want non-negative input
        X -= X.min()
        for name, Estimator in estimators:
            if (name in dont_test
                    or name in ['CCA', '_CCA', 'PLSCanonical',
                                'PLSRegression',
                                'PLSSVD', 'GaussianProcess']):
                # FIXME!
                # in particular GaussianProcess!
                continue

            yield (check_estimators_overwrite_params, name, Estimator, X,
                   multioutput_estimator_convert_y_2d(name, y))


def check_estimators_overwrite_params(name, Estimator, X, y):
    with warnings.catch_warnings(record=True):
        # catch deprecation warnings
        estimator = Estimator()

    if hasattr(estimator, 'batch_size'):
        # FIXME
        # for MiniBatchDictLearning
        estimator.batch_size = 1

    if name in ['GaussianRandomProjection',
                'SparseRandomProjection']:
        # Due to the jl lemma and very few samples, the number
        # of components of the random matrix projection will be
        # greater
        # than the number of features.
        # So we impose a smaller number (avoid "auto" mode)
        estimator = Estimator(n_components=1)
    elif name == "SelectKBest":
        estimator = Estimator(k=1)

    set_random_state(estimator)

    params = estimator.get_params()
    estimator.fit(X, y)
    new_params = estimator.get_params()
    for k, v in params.items():
        assert_false(np.any(new_params[k] != v),
                     "Estimator %s changes its parameter %s"
                     " from %s to %s during fit."
                     % (name, k, v, new_params[k]))


def test_cluster_overwrite_params():
    # test whether any classifier overwrites his init parameters during fit
    clusterers = all_estimators(type_filter="cluster")
    X, y = make_blobs(random_state=0, n_samples=9)
    # some want non-negative input
    X
    for name, Clustering in clusterers:
        yield check_cluster_overwrite_params, name, Clustering, X, y


def check_cluster_overwrite_params(name, Clustering, X, y):
    with warnings.catch_warnings(record=True):
        # catch deprecation warnings
        clustering = Clustering()
    params = clustering.get_params()
    clustering.fit(X)
    new_params = clustering.get_params()
    for k, v in params.items():
        assert_false(np.any(new_params[k] != v),
                     "Estimator %s changes its parameter %s"
                     " from %s to %s during fit."
                     % (name, k, v, new_params[k]))


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
    """Test if predict with sparsified estimators works.

    Tests regression, binary classification, and multi-class classification.
    """
    estimators = all_estimators()
    X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
    y = [1, 1, 1, 2, 2, 2]

    # test regression and binary classification
    for name, Estimator in estimators:
        try:
            Estimator.sparsify
        except:
            continue
        yield check_sparsify_binary_classifier, name, Estimator, X, y

    # test multiclass classification
    classifiers = all_estimators(type_filter='classifier')
    y[-1] = 3  # make multi-class
    for name, Classifier in classifiers:
        try:
            Classifier.sparsify
        except:
            continue
        yield check_sparsify_multiclass_classifier, name, Classifier, X, y


def check_sparsify_multiclass_classifier(name, Classifier, X, y):
    est = Classifier()

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


def check_sparsify_binary_classifier(name, Estimator, X, y):
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


def check_estimators_not_an_array(name, Estimator, X, y):
    if name in ('CCA', '_PLS', 'PLSCanonical', 'PLSRegression'):
        raise SkipTest
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        # separate estimators to control random seeds
        regressor_1 = Estimator()
        regressor_2 = Estimator()
    set_random_state(regressor_1)
    set_random_state(regressor_2)

    y_ = NotAnArray(np.asarray(y))

    # fit
    regressor_1.fit(X, y_)
    pred1 = regressor_1.predict(X)
    regressor_2.fit(X, y)
    pred2 = regressor_2.predict(X)
    assert_array_almost_equal(pred1, pred2, 2, name)


def test_regressors_not_an_array():
    regressors = all_estimators(type_filter='regressor')
    X, y = _boston_subset(n_samples=50)
    X = StandardScaler().fit_transform(X)

    for name, Regressor in regressors:
        if name in dont_test:
            continue
        yield (check_estimators_not_an_array, name, Regressor, X,
               multioutput_estimator_convert_y_2d(name, y))


def test_classifiers_not_an_array():
    classifiers = all_estimators(type_filter='classifier')
    X = np.array([[3, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 1]])
    y = [1, 1, 1, 2, 2, 2]

    for name, Classifier in classifiers:
        if name in dont_test:
            continue
        yield (check_estimators_not_an_array, name, Classifier, X,
               multioutput_estimator_convert_y_2d(name, y))


