from __future__ import print_function

import warnings
import sys
import traceback
import inspect
import pickle

import numpy as np
from scipy import sparse
import struct

from sklearn.externals.six.moves import zip
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import META_ESTIMATORS
from sklearn.utils.testing import set_random_state
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import check_skip_travis

from sklearn.base import (clone, ClusterMixin)
from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

from sklearn.lda import LDA
from sklearn.random_projection import BaseRandomProjection
from sklearn.feature_selection import SelectKBest
from sklearn.svm.base import BaseLibSVM

from sklearn.utils.validation import DataConversionWarning
from sklearn.cross_validation import train_test_split

from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris, load_boston, make_blobs


BOSTON = None
CROSS_DECOMPOSITION = ['PLSCanonical', 'PLSRegression', 'CCA', 'PLSSVD']


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


def set_fast_parameters(estimator):
    # speed up some estimators
    params = estimator.get_params()
    if "n_iter" in params:
        estimator.set_params(n_iter=5)
    if "max_iter" in params:
        # NMF
        if estimator.max_iter is not None:
            estimator.set_params(max_iter=min(5, estimator.max_iter))
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


class NotAnArray(object):
    " An object that is convertable to an array"

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype=None):
        return self.data


def _is_32bit():
    """Detect if process is 32bit Python."""
    return struct.calcsize('P') * 8 == 32


def check_regressors_classifiers_sparse_data(name, Estimator):
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    # catch deprecation warnings
    with warnings.catch_warnings():
        estimator = Estimator()
    set_fast_parameters(estimator)
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


def check_transformer(name, Transformer):
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

    if name == "KernelPCA":
        transformer.remove_zero_eig = False

    set_fast_parameters(transformer)

    # fit

    if name in CROSS_DECOMPOSITION:
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

        # raises error on malformed input for transform
        if hasattr(X, 'T'):
            # If it's not an array, it does not have a 'T' property
            assert_raises(ValueError, transformer.transform, X.T)


def check_transformer_sparse_data(name, Transformer):
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        if name in ['Scaler', 'StandardScaler']:
            transformer = Transformer(with_mean=False)
        else:
            transformer = Transformer()

    set_fast_parameters(transformer)

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


def check_estimators_nan_inf(name, Estimator):
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
            set_fast_parameters(estimator)
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


def check_transformer_pickle(name, Transformer):
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    n_samples, n_features = X.shape
    X = StandardScaler().fit_transform(X)
    X -= X.min()
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        transformer = Transformer()
    if not hasattr(transformer, 'transform'):
        return
    set_random_state(transformer)
    set_fast_parameters(transformer)

    # fit
    if name in CROSS_DECOMPOSITION:
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


def check_clustering(name, Alg):
    X, y = make_blobs(n_samples=50, random_state=1)
    X, y = shuffle(X, y, random_state=7)
    X = StandardScaler().fit_transform(X)
    n_samples, n_features = X.shape
    # catch deprecation and neighbors warnings
    with warnings.catch_warnings(record=True):
        alg = Alg()
    set_fast_parameters(alg)
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
        set_fast_parameters(classifier)
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


def check_classifiers_train(name, Classifier):
    X_m, y_m = make_blobs(random_state=0)
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
        set_fast_parameters(classifier)
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
            assert_greater(accuracy_score(y, y_pred), 0.85)

        # raises error on malformed input for predict
        assert_raises(ValueError, classifier.predict, X.T)
        if hasattr(classifier, "decision_function"):
            try:
                # decision_function agrees with predict:
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
            # predict_proba agrees with predict:
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


def check_classifiers_input_shapes(name, Classifier):
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=1)
    X = StandardScaler().fit_transform(X)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    set_fast_parameters(classifier)
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
        set_fast_parameters(classifier)
        # fit
        classifier.fit(X, y_)

        y_pred = classifier.predict(X)
        # training set performance
        assert_array_equal(np.unique(y_), np.unique(y_pred))
        if np.any(classifier.classes_ != classes):
            print("Unexpected classes_ attribute for %r: "
                  "expected %s, got %s" %
                  (classifier, classes, classifier.classes_))


def check_classifiers_pickle(name, Classifier):
    X, y = make_blobs(random_state=0)
    X, y = shuffle(X, y, random_state=7)
    X -= X.min()

    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    set_fast_parameters(classifier)
    # raises error on malformed input for fit
    assert_raises(ValueError, classifier.fit, X, y[:-1])

    # fit
    classifier.fit(X, y)
    y_pred = classifier.predict(X)
    pickled_classifier = pickle.dumps(classifier)
    unpickled_classifier = pickle.loads(pickled_classifier)
    pickled_y_pred = unpickled_classifier.predict(X)

    assert_array_almost_equal(pickled_y_pred, y_pred)


def check_regressors_int(name, Regressor):
    X, _ = _boston_subset()
    X = X[:50]
    rnd = np.random.RandomState(0)
    y = rnd.randint(3, size=X.shape[0])
    y = multioutput_estimator_convert_y_2d(name, y)
    if name == 'OrthogonalMatchingPursuitCV':
        # FIXME: This test is unstable on Travis, see issue #3190.
        check_skip_travis()
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        # separate estimators to control random seeds
        regressor_1 = Regressor()
        regressor_2 = Regressor()
    set_fast_parameters(regressor_1)
    set_fast_parameters(regressor_2)
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
    y = StandardScaler().fit_transform(y)   # X is already scaled
    y = multioutput_estimator_convert_y_2d(name, y)
    if name == 'OrthogonalMatchingPursuitCV':
        # FIXME: This test is unstable on Travis, see issue #3190.
        check_skip_travis()
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        regressor = Regressor()
    set_fast_parameters(regressor)
    if not hasattr(regressor, 'alphas') and hasattr(regressor, 'alpha'):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01

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
    regressor.predict(X)

      # TODO: find out why PLS and CCA fail. RANSAC is random
      # and furthermore assumes the presence of outliers, hence
      # skipped
    if name not in ('PLSCanonical', 'CCA', 'RANSACRegressor'):
        assert_greater(regressor.score(X, y_), 0.5)


def check_regressors_pickle(name, Regressor):
    X, y = _boston_subset()
    y = StandardScaler().fit_transform(y)   # X is already scaled
    y = multioutput_estimator_convert_y_2d(name, y)
    if name == 'OrthogonalMatchingPursuitCV':
        # FIXME: This test is unstable on Travis, see issue #3190.
        check_skip_travis()
    rnd = np.random.RandomState(0)
    # catch deprecation warnings
    with warnings.catch_warnings(record=True):
        regressor = Regressor()
    set_fast_parameters(regressor)
    if not hasattr(regressor, 'alphas') and hasattr(regressor, 'alpha'):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01

    if name in CROSS_DECOMPOSITION:
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


def check_class_weight_classifiers(name, Classifier):
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

        set_random_state(classifier)
        classifier.fit(X_train, y_train)
        y_pred = classifier.predict(X_test)
        assert_greater(np.mean(y_pred == 0), 0.9)


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


def check_class_weight_auto_linear_classifier(name, Classifier):
    """Test class weights with non-contiguous class labels."""
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                  [1.0, 1.0], [1.0, 0.0]])
    y = [1, 1, 1, -1, -1]

    with warnings.catch_warnings(record=True):
        classifier = Classifier()
    if hasattr(classifier, "n_iter"):
        # This is a very small dataset, default n_iter are likely to prevent
        # convergence
        classifier.set_params(n_iter=1000)
    set_random_state(classifier)

    # Let the model compute the class frequencies
    classifier.set_params(class_weight='auto')
    coef_auto = classifier.fit(X, y).coef_.copy()

    # Count each label occurrence to reweight manually
    mean_weight = (1. / 3 + 1. / 2) / 2
    class_weight = {
        1: 1. / 3 / mean_weight,
        -1: 1. / 2 / mean_weight,
    }
    classifier.set_params(class_weight=class_weight)
    coef_manual = classifier.fit(X, y).coef_.copy()

    assert_array_almost_equal(coef_auto, coef_manual)


def check_estimators_overwrite_params(name, Estimator):
    X, y = make_blobs(random_state=0, n_samples=9)
    y = multioutput_estimator_convert_y_2d(name, y)
    # some want non-negative input
    X -= X.min()
    with warnings.catch_warnings(record=True):
        # catch deprecation warnings
        estimator = Estimator()

    if name == 'MiniBatchDictLearning' or name == 'MiniBatchSparsePCA':
        # FIXME
        # for MiniBatchDictLearning and MiniBatchSparsePCA
        estimator.batch_size = 1

    set_fast_parameters(estimator)

    set_random_state(estimator)

    params = estimator.get_params()
    estimator.fit(X, y)
    new_params = estimator.get_params()
    for k, v in params.items():
        assert_false(np.any(new_params[k] != v),
                     "Estimator %s changes its parameter %s"
                     " from %s to %s during fit."
                     % (name, k, v, new_params[k]))


def check_cluster_overwrite_params(name, Clustering):
    X, y = make_blobs(random_state=0, n_samples=9)
    with warnings.catch_warnings(record=True):
        # catch deprecation warnings
        clustering = Clustering()
    set_fast_parameters(clustering)
    params = clustering.get_params()
    clustering.fit(X)
    new_params = clustering.get_params()
    for k, v in params.items():
        assert_false(np.any(new_params[k] != v),
                     "Estimator %s changes its parameter %s"
                     " from %s to %s during fit."
                     % (name, k, v, new_params[k]))


def check_sparsify_multiclass_classifier(name, Classifier):
    X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1],
                  [-1, -2], [2, 2], [-2, -2]])
    y = [1, 1, 1, 2, 2, 2, 3, 3, 3]
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


def check_sparsify_binary_classifier(name, Estimator):
    X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
    y = [1, 1, 1, 2, 2, 2]
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
    set_fast_parameters(estimator_1)
    set_fast_parameters(estimator_2)
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
    classifier = LDA()
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
        if name in META_ESTIMATORS:
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


def multioutput_estimator_convert_y_2d(name, y):
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if name in (['MultiTaskElasticNetCV', 'MultiTaskLassoCV',
                 'MultiTaskLasso', 'MultiTaskElasticNet']):
        return y[:, np.newaxis]
    return y


def check_non_transformer_estimators_n_iter(name, estimator, multi_output=False):
    # Check if all iterative solvers, run for more than one iteratiom

    iris = load_iris()
    X, y_ = iris.data, iris.target

    if multi_output:
        y_ = y_[:, np.newaxis]

    set_random_state(estimator, 0)
    if name == 'AffinityPropagation':
        estimator.fit(X)
    else:
        estimator.fit(X, y_)
    assert_greater(estimator.n_iter_, 0)


def check_transformer_n_iter(name, estimator):
    if name in CROSS_DECOMPOSITION:
        # Check using default data
        X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
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
