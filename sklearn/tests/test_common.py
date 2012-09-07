"""
General tests for all estimators in sklearn.
"""
import os
import warnings
import sys
import traceback

import numpy as np
from scipy import sparse
from nose.tools import assert_raises, assert_equal, assert_true
from numpy.testing import assert_array_equal, \
        assert_array_almost_equal

import sklearn
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_greater
from sklearn.base import clone, ClassifierMixin, RegressorMixin, \
        TransformerMixin, ClusterMixin
from sklearn.utils import shuffle
from sklearn.preprocessing import Scaler
#from sklearn.cross_validation import train_test_split
from sklearn.datasets import load_iris, load_boston, make_blobs
from sklearn.metrics import zero_one_score, adjusted_rand_score
from sklearn.lda import LDA
from sklearn.svm.base import BaseLibSVM

# import "special" estimators
from sklearn.grid_search import GridSearchCV
from sklearn.decomposition import SparseCoder
from sklearn.pipeline import Pipeline
from sklearn.pls import _PLS, PLSCanonical, PLSRegression, CCA, PLSSVD
from sklearn.ensemble import BaseEnsemble
from sklearn.multiclass import OneVsOneClassifier, OneVsRestClassifier,\
        OutputCodeClassifier
from sklearn.feature_selection import RFE, RFECV, SelectKBest
from sklearn.naive_bayes import MultinomialNB, BernoulliNB
from sklearn.covariance import EllipticEnvelope, EllipticEnvelop
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.kernel_approximation import AdditiveChi2Sampler
from sklearn.preprocessing import LabelBinarizer, LabelEncoder, Binarizer, \
        Normalizer
from sklearn.cluster import WardAgglomeration, AffinityPropagation, \
        SpectralClustering

dont_test = [Pipeline, GridSearchCV, SparseCoder, EllipticEnvelope,
        EllipticEnvelop, DictVectorizer, LabelBinarizer, LabelEncoder,
        TfidfTransformer]
meta_estimators = [BaseEnsemble, OneVsOneClassifier, OutputCodeClassifier,
        OneVsRestClassifier, RFE, RFECV]


def test_all_estimators():
    # Test that estimators are default-constructible, clonable
    # and have working repr.
    estimators = all_estimators()
    clf = LDA()

    for name, E in estimators:
        # some can just not be sensibly default constructed
        if E in dont_test:
            continue
        # test default-constructibility
        # get rid of deprecation warnings
        with warnings.catch_warnings(record=True):
            if E in meta_estimators:
                e = E(clf)
            else:
                e = E()
            #test cloning
            clone(e)
            # test __repr__
            repr(e)


def test_estimators_sparse_data():
    # All estimators should either deal with sparse data, or raise an
    # intelligible error message
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    estimators = all_estimators()
    estimators = [(name, E) for name, E in estimators
                        if issubclass(E, (ClassifierMixin, RegressorMixin))]
    for name, Clf in estimators:
        if Clf in dont_test or Clf in meta_estimators:
            continue
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            clf = Clf()
        # fit
        try:
            clf.fit(X, y)
        except TypeError, e:
            if not 'sparse' in repr(e):
                print ("Estimator %s doesn't seem to fail gracefully on "
                    "sparse data" % name)
                traceback.print_exc(file=sys.stdout)
                raise e
        except Exception, exc:
            print ("Estimator %s doesn't seem to fail gracefully on "
                "sparse data" % name)
            traceback.print_exc(file=sys.stdout)
            raise exc


def test_transformers():
    # test if transformers do something sensible on training set
    # also test all shapes / shape errors
    estimators = all_estimators()
    transformers = [(name, E) for name, E in estimators if issubclass(E,
        TransformerMixin)]
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
            random_state=0, n_features=2, cluster_std=0.1)
    n_samples, n_features = X.shape
    X = Scaler().fit_transform(X)
    X -= X.min()

    succeeded = True

    for name, Trans in transformers:
        if Trans in dont_test or Trans in meta_estimators:
            continue
        # these don't actually fit the data:
        if Trans in [AdditiveChi2Sampler, Binarizer, Normalizer]:
            continue
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            trans = Trans()

        if hasattr(trans, 'compute_importances'):
            trans.compute_importances = True

        if Trans is SelectKBest:
            # SelectKBest has a default of k=10
            # which is more feature than we have.
            trans.k = 1

        # fit

        if Trans in (_PLS, PLSCanonical, PLSRegression, CCA, PLSSVD):
            y_ = np.vstack([y, 2 * y + np.random.randint(2, size=len(y))])
            y_ = y_.T
        else:
            y_ = y

        try:
            trans.fit(X, y_)
            X_pred = trans.fit_transform(X, y=y_)
            if isinstance(X_pred, tuple):
                for x_pred in X_pred:
                    assert_equal(x_pred.shape[0], n_samples)
            else:
                assert_equal(X_pred.shape[0], n_samples)
        except Exception as e:
            print trans
            print e
            print
            succeeded = False

        if hasattr(trans, 'transform'):
            if Trans in (_PLS, PLSCanonical, PLSRegression, CCA, PLSSVD):
                X_pred2 = trans.transform(X, y_)
            else:
                X_pred2 = trans.transform(X)
            if isinstance(X_pred, tuple) and isinstance(X_pred2, tuple):
                for x_pred, x_pred2 in zip(X_pred, X_pred2):
                    assert_array_almost_equal(x_pred, x_pred2, 2,
                        "fit_transform not correct in %s" % Trans)
            else:
                assert_array_almost_equal(X_pred, X_pred2, 2,
                    "fit_transform not correct in %s" % Trans)

            # raises error on malformed input for transform
            assert_raises(ValueError, trans.transform, X.T)
    assert_true(succeeded)


def test_transformers_sparse_data():
    # All estimators should either deal with sparse data, or raise an
    # intelligible error message
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    estimators = all_estimators()
    estimators = [(name, E) for name, E in estimators
                        if issubclass(E, TransformerMixin)]
    for name, Trans in estimators:
        if Trans in dont_test or Trans in meta_estimators:
            continue
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            if Trans is Scaler:
                trans = Trans(with_mean=False)
            else:
                trans = Trans()
        # fit
        try:
            trans.fit(X, y)
        except TypeError, e:
            if not 'sparse' in repr(e):
                print ("Estimator %s doesn't seem to fail gracefully on "
                    "sparse data" % name)
                traceback.print_exc(file=sys.stdout)
                raise e
        except Exception, exc:
            print ("Estimator %s doesn't seem to fail gracefully on "
                "sparse data" % name)
            traceback.print_exc(file=sys.stdout)
            raise exc


def test_classifiers_one_label():
    # test classifiers trained on a single label always return this label
    # or raise an sensible error message
    rnd = np.random.RandomState(0)
    X_train = rnd.uniform(size=(10, 3))
    X_test = rnd.uniform(size=(10, 3))
    y = np.ones(10)
    estimators = all_estimators()
    classifiers = [(name, E) for name, E in estimators if issubclass(E,
        ClassifierMixin)]
    error_string_fit = "Classifier can't train when only one class is present."
    error_string_predict = ("Classifier can't predict when only one class is "
        "present.")
    for name, Clf in classifiers:
        if Clf in dont_test or Clf in meta_estimators:
            continue
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            clf = Clf()
            # try to fit
            try:
                clf.fit(X_train, y)
            except ValueError, e:
                if not 'class' in repr(e):
                    print(error_string_fit, Clf, e)
                    traceback.print_exc(file=sys.stdout)
                    raise e
                else:
                    continue
            except Exception, exc:
                    print(error_string_fit, Clf, exc)
                    traceback.print_exc(file=sys.stdout)
                    raise exc
            # predict
            try:
                assert_array_equal(clf.predict(X_test), y)
            except Exception, exc:
                print(error_string_predict, Clf, exc)
                traceback.print_exc(file=sys.stdout)


def test_clustering():
    # test if clustering algorithms do something sensible
    # also test all shapes / shape errors
    estimators = all_estimators()
    clustering = [(name, E) for name, E in estimators if issubclass(E,
        ClusterMixin)]
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=7)
    n_samples, n_features = X.shape
    X = Scaler().fit_transform(X)
    for name, Alg in clustering:
        if Alg is WardAgglomeration:
            # this is clustering on the features
            # let's not test that here.
            continue
        # catch deprecation and neighbors warnings
        with warnings.catch_warnings(record=True):
            alg = Alg()
            if hasattr(alg, "n_clusters"):
                alg.set_params(n_clusters=3)
            if hasattr(alg, "random_state"):
                alg.set_params(random_state=1)
            if Alg is AffinityPropagation:
                alg.set_params(preference=-100)
            # fit
            alg.fit(X)

        assert_equal(alg.labels_.shape, (n_samples,))
        pred = alg.labels_
        assert_greater(adjusted_rand_score(pred, y), 0.4)
        # fit another time with ``fit_predict`` and compare results
        if Alg is SpectralClustering:
            # there is no way to make Spectral clustering deterministic :(
            continue
        if hasattr(alg, "random_state"):
            alg.set_params(random_state=1)
        with warnings.catch_warnings(record=True):
            pred2 = alg.fit_predict(X)
        assert_array_equal(pred, pred2)


def test_classifiers_train():
    # test if classifiers do something sensible on training set
    # also test all shapes / shape errors
    estimators = all_estimators()
    classifiers = [(name, E) for name, E in estimators if issubclass(E,
        ClassifierMixin)]
    iris = load_iris()
    X_m, y_m = iris.data, iris.target
    X_m, y_m = shuffle(X_m, y_m, random_state=7)
    X_m = Scaler().fit_transform(X_m)
    # generate binary problem from multi-class one
    y_b = y_m[y_m != 2]
    X_b = X_m[y_m != 2]
    for (X, y) in [(X_m, y_m), (X_b, y_b)]:
        # do it once with binary, once with multiclass
        n_labels = len(np.unique(y))
        n_samples, n_features = X.shape
        for name, Clf in classifiers:
            if Clf in dont_test or Clf in meta_estimators:
                continue
            if Clf in [MultinomialNB, BernoulliNB]:
                # TODO also test these!
                continue
            # catch deprecation warnings
            with warnings.catch_warnings(record=True):
                clf = Clf()
            # raises error on malformed input for fit
            assert_raises(ValueError, clf.fit, X, y[:-1])

            # fit
            clf.fit(X, y)
            y_pred = clf.predict(X)
            assert_equal(y_pred.shape, (n_samples,))
            # training set performance
            assert_greater(zero_one_score(y, y_pred), 0.78)

            # raises error on malformed input for predict
            assert_raises(ValueError, clf.predict, X.T)
            if hasattr(clf, "decision_function"):
                try:
                    # decision_function agrees with predict:
                    decision = clf.decision_function(X)
                    if n_labels is 2:
                        assert_equal(decision.ravel().shape, (n_samples,))
                        dec_pred = (decision.ravel() > 0).astype(np.int)
                        assert_array_equal(dec_pred, y_pred)
                    if n_labels is 3 and not isinstance(clf, BaseLibSVM):
                        # 1on1 of LibSVM works differently
                        assert_equal(decision.shape, (n_samples, n_labels))
                        assert_array_equal(np.argmax(decision, axis=1), y_pred)

                    # raises error on malformed input
                    assert_raises(ValueError, clf.decision_function, X.T)
                    # raises error on malformed input for decision_function
                    assert_raises(ValueError, clf.decision_function, X.T)
                except NotImplementedError:
                    pass
            if hasattr(clf, "predict_proba"):
                try:
                    # predict_proba agrees with predict:
                    y_prob = clf.predict_proba(X)
                    assert_equal(y_prob.shape, (n_samples, n_labels))
                    # raises error on malformed input
                    assert_raises(ValueError, clf.predict_proba, X.T)
                    assert_array_equal(np.argmax(y_prob, axis=1), y_pred)
                    # raises error on malformed input for predict_proba
                    assert_raises(ValueError, clf.predict_proba, X.T)
                except NotImplementedError:
                    pass


def test_classifiers_classes():
    # test if classifiers can cope with non-consecutive classes
    estimators = all_estimators()
    classifiers = [(name, E) for name, E in estimators if issubclass(E,
        ClassifierMixin)]
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=7)
    X = Scaler().fit_transform(X)
    y = 2 * y + 1
    # TODO: make work with next line :)
    #y = y.astype(np.str)
    for name, Clf in classifiers:
        if Clf in dont_test or Clf in meta_estimators:
            continue
        if Clf in [MultinomialNB, BernoulliNB]:
            # TODO also test these!
            continue

        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            clf = Clf()
        # fit
        clf.fit(X, y)
        y_pred = clf.predict(X)
        # training set performance
        assert_array_equal(np.unique(y), np.unique(y_pred))
        assert_greater(zero_one_score(y, y_pred), 0.78)


def test_regressors_int():
    # test if regressors can cope with integer labels (by converting them to
    # float)
    estimators = all_estimators()
    regressors = [(name, E) for name, E in estimators if issubclass(E,
        RegressorMixin)]
    boston = load_boston()
    X, y = boston.data, boston.target
    X, y = shuffle(X, y, random_state=0)
    X = Scaler().fit_transform(X)
    y = np.random.randint(2, size=X.shape[0])
    for name, Reg in regressors:
        if Reg in dont_test or Reg in meta_estimators or Reg in (CCA,):
            continue
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            # separate estimators to control random seeds
            reg1 = Reg()
            reg2 = Reg()
        if hasattr(reg1, 'alpha'):
            reg1.set_params(alpha=0.01)
            reg2.set_params(alpha=0.01)
        if hasattr(reg1, 'random_state'):
            reg1.set_params(random_state=0)
            reg2.set_params(random_state=0)

        if Reg in (_PLS, PLSCanonical, PLSRegression):
            y_ = np.vstack([y, 2 * y + np.random.randint(2, size=len(y))])
            y_ = y_.T
        else:
            y_ = y

        # fit
        reg1.fit(X, y_)
        pred1 = reg1.predict(X)
        reg2.fit(X, y_.astype(np.float))
        pred2 = reg2.predict(X)
        assert_array_almost_equal(pred1, pred2, 2, name)


def test_regressors_train():
    estimators = all_estimators()
    regressors = [(name, E) for name, E in estimators if issubclass(E,
        RegressorMixin)]
    boston = load_boston()
    X, y = boston.data, boston.target
    X, y = shuffle(X, y, random_state=0)
    # TODO: test with intercept
    # TODO: test with multiple responses
    X = Scaler().fit_transform(X)
    y = Scaler().fit_transform(y)
    succeeded = True
    for name, Reg in regressors:
        if Reg in dont_test or Reg in meta_estimators:
            continue
        # catch deprecation warnings
        with warnings.catch_warnings(record=True):
            reg = Reg()
        if hasattr(reg, 'alpha'):
            reg.set_params(alpha=0.01)

        # raises error on malformed input for fit
        assert_raises(ValueError, reg.fit, X, y[:-1])
        # fit
        try:
            if Reg in (_PLS, PLSCanonical, PLSRegression, CCA):
                y_ = np.vstack([y, 2 * y + np.random.randint(2, size=len(y))])
                y_ = y_.T
            else:
                y_ = y
            reg.fit(X, y_)
            reg.predict(X)

            if Reg not in (PLSCanonical, CCA):  # TODO: find out why
                assert_greater(reg.score(X, y_), 0.5)
        except Exception as e:
            print(reg)
            print e
            print
            succeeded = False

    assert_true(succeeded)


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
            warnings.simplefilter('ignore',  UserWarning)
            execfile('setup.py', dict(__name__='__main__'))
    finally:
        sys.argv = old_argv
        os.chdir(cwd)
