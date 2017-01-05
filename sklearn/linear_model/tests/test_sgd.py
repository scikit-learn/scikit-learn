import pickle
import unittest

import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import raises
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_false, assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import ignore_warnings

from sklearn import linear_model, datasets, metrics
from sklearn.base import clone
from sklearn.linear_model import SGDClassifier, SGDRegressor
from sklearn.preprocessing import LabelEncoder, scale, MinMaxScaler

from sklearn.linear_model import sgd_fast


class SparseSGDClassifier(SGDClassifier):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return super(SparseSGDClassifier, self).fit(X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return super(SparseSGDClassifier, self).partial_fit(X, y, *args, **kw)

    def decision_function(self, X):
        X = sp.csr_matrix(X)
        return super(SparseSGDClassifier, self).decision_function(X)

    def predict_proba(self, X):
        X = sp.csr_matrix(X)
        return super(SparseSGDClassifier, self).predict_proba(X)


class SparseSGDRegressor(SGDRegressor):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDRegressor.fit(self, X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDRegressor.partial_fit(self, X, y, *args, **kw)

    def decision_function(self, X, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDRegressor.decision_function(self, X, *args, **kw)


# Test Data

# test sample 1
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
Y = [1, 1, 1, 2, 2, 2]
T = np.array([[-1, -1], [2, 2], [3, 2]])
true_result = [1, 2, 2]

# test sample 2; string class labels
X2 = np.array([[-1, 1], [-0.75, 0.5], [-1.5, 1.5],
               [1, 1], [0.75, 0.5], [1.5, 1.5],
               [-1, -1], [0, -0.5], [1, -1]])
Y2 = ["one"] * 3 + ["two"] * 3 + ["three"] * 3
T2 = np.array([[-1.5, 0.5], [1, 2], [0, -2]])
true_result2 = ["one", "two", "three"]

# test sample 3
X3 = np.array([[1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0],
               [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 1],
               [0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0]])
Y3 = np.array([1, 1, 1, 1, 2, 2, 2, 2])

# test sample 4 - two more or less redundant feature groups
X4 = np.array([[1, 0.9, 0.8, 0, 0, 0], [1, .84, .98, 0, 0, 0],
               [1, .96, .88, 0, 0, 0], [1, .91, .99, 0, 0, 0],
               [0, 0, 0, .89, .91, 1], [0, 0, 0, .79, .84, 1],
               [0, 0, 0, .91, .95, 1], [0, 0, 0, .93, 1, 1]])
Y4 = np.array([1, 1, 1, 1, 2, 2, 2, 2])

iris = datasets.load_iris()

# test sample 5 - test sample 1 as binary classification problem
X5 = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
Y5 = [1, 1, 1, 2, 2, 2]
true_result5 = [0, 1, 1]


# Classification Test Case

class CommonTest(object):

    def factory(self, **kwargs):
        if "random_state" not in kwargs:
            kwargs["random_state"] = 42
        return self.factory_class(**kwargs)

    # a simple implementation of ASGD to use for testing
    # uses squared loss to find the gradient
    def asgd(self, X, y, eta, alpha, weight_init=None, intercept_init=0.0):
        if weight_init is None:
            weights = np.zeros(X.shape[1])
        else:
            weights = weight_init

        average_weights = np.zeros(X.shape[1])
        intercept = intercept_init
        average_intercept = 0.0
        decay = 1.0

        # sparse data has a fixed decay of .01
        if (isinstance(self, SparseSGDClassifierTestCase) or
                isinstance(self, SparseSGDRegressorTestCase)):
            decay = .01

        for i, entry in enumerate(X):
            p = np.dot(entry, weights)
            p += intercept
            gradient = p - y[i]
            weights *= 1.0 - (eta * alpha)
            weights += -(eta * gradient * entry)
            intercept += -(eta * gradient) * decay

            average_weights *= i
            average_weights += weights
            average_weights /= i + 1.0

            average_intercept *= i
            average_intercept += intercept
            average_intercept /= i + 1.0

        return average_weights, average_intercept

    def _test_warm_start(self, X, Y, lr):
        # Test that explicit warm restart...
        clf = self.factory(alpha=0.01, eta0=0.01, n_iter=5, shuffle=False,
                           learning_rate=lr)
        clf.fit(X, Y)

        clf2 = self.factory(alpha=0.001, eta0=0.01, n_iter=5, shuffle=False,
                            learning_rate=lr)
        clf2.fit(X, Y,
                 coef_init=clf.coef_.copy(),
                 intercept_init=clf.intercept_.copy())

        # ... and implicit warm restart are equivalent.
        clf3 = self.factory(alpha=0.01, eta0=0.01, n_iter=5, shuffle=False,
                            warm_start=True, learning_rate=lr)
        clf3.fit(X, Y)

        assert_equal(clf3.t_, clf.t_)
        assert_array_almost_equal(clf3.coef_, clf.coef_)

        clf3.set_params(alpha=0.001)
        clf3.fit(X, Y)

        assert_equal(clf3.t_, clf2.t_)
        assert_array_almost_equal(clf3.coef_, clf2.coef_)

    def test_warm_start_constant(self):
        self._test_warm_start(X, Y, "constant")

    def test_warm_start_invscaling(self):
        self._test_warm_start(X, Y, "invscaling")

    def test_warm_start_optimal(self):
        self._test_warm_start(X, Y, "optimal")

    def test_input_format(self):
        # Input format tests.
        clf = self.factory(alpha=0.01, n_iter=5,
                           shuffle=False)
        clf.fit(X, Y)
        Y_ = np.array(Y)[:, np.newaxis]

        Y_ = np.c_[Y_, Y_]
        assert_raises(ValueError, clf.fit, X, Y_)

    def test_clone(self):
        # Test whether clone works ok.
        clf = self.factory(alpha=0.01, n_iter=5, penalty='l1')
        clf = clone(clf)
        clf.set_params(penalty='l2')
        clf.fit(X, Y)

        clf2 = self.factory(alpha=0.01, n_iter=5, penalty='l2')
        clf2.fit(X, Y)

        assert_array_equal(clf.coef_, clf2.coef_)

    def test_plain_has_no_average_attr(self):
        clf = self.factory(average=True, eta0=.01)
        clf.fit(X, Y)

        assert_true(hasattr(clf, 'average_coef_'))
        assert_true(hasattr(clf, 'average_intercept_'))
        assert_true(hasattr(clf, 'standard_intercept_'))
        assert_true(hasattr(clf, 'standard_coef_'))

        clf = self.factory()
        clf.fit(X, Y)

        assert_false(hasattr(clf, 'average_coef_'))
        assert_false(hasattr(clf, 'average_intercept_'))
        assert_false(hasattr(clf, 'standard_intercept_'))
        assert_false(hasattr(clf, 'standard_coef_'))

    def test_late_onset_averaging_not_reached(self):
        clf1 = self.factory(average=600)
        clf2 = self.factory()
        for _ in range(100):
            if isinstance(clf1, SGDClassifier):
                clf1.partial_fit(X, Y, classes=np.unique(Y))
                clf2.partial_fit(X, Y, classes=np.unique(Y))
            else:
                clf1.partial_fit(X, Y)
                clf2.partial_fit(X, Y)

        assert_array_almost_equal(clf1.coef_, clf2.coef_, decimal=16)
        assert_almost_equal(clf1.intercept_, clf2.intercept_, decimal=16)

    def test_late_onset_averaging_reached(self):
        eta0 = .001
        alpha = .0001
        Y_encode = np.array(Y)
        Y_encode[Y_encode == 1] = -1.0
        Y_encode[Y_encode == 2] = 1.0

        clf1 = self.factory(average=7, learning_rate="constant",
                            loss='squared_loss', eta0=eta0,
                            alpha=alpha, n_iter=2, shuffle=False)
        clf2 = self.factory(average=0, learning_rate="constant",
                            loss='squared_loss', eta0=eta0,
                            alpha=alpha, n_iter=1, shuffle=False)

        clf1.fit(X, Y_encode)
        clf2.fit(X, Y_encode)

        average_weights, average_intercept = \
            self.asgd(X, Y_encode, eta0, alpha,
                      weight_init=clf2.coef_.ravel(),
                      intercept_init=clf2.intercept_)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  average_weights.ravel(),
                                  decimal=16)
        assert_almost_equal(clf1.intercept_, average_intercept, decimal=16)

    @raises(ValueError)
    def test_sgd_bad_alpha_for_optimal_learning_rate(self):
        # Check whether expected ValueError on bad alpha, i.e. 0
        # since alpha is used to compute the optimal learning rate
        self.factory(alpha=0, learning_rate="optimal")


class DenseSGDClassifierTestCase(unittest.TestCase, CommonTest):
    """Test suite for the dense representation variant of SGD"""
    factory_class = SGDClassifier

    def test_sgd(self):
        # Check that SGD gives any results :-)

        for loss in ("hinge", "squared_hinge", "log", "modified_huber"):
            clf = self.factory(penalty='l2', alpha=0.01, fit_intercept=True,
                               loss=loss, n_iter=10, shuffle=True)
            clf.fit(X, Y)
            # assert_almost_equal(clf.coef_[0], clf.coef_[1], decimal=7)
            assert_array_equal(clf.predict(T), true_result)

    @raises(ValueError)
    def test_sgd_bad_l1_ratio(self):
        # Check whether expected ValueError on bad l1_ratio
        self.factory(l1_ratio=1.1)

    @raises(ValueError)
    def test_sgd_bad_learning_rate_schedule(self):
        # Check whether expected ValueError on bad learning_rate
        self.factory(learning_rate="<unknown>")

    @raises(ValueError)
    def test_sgd_bad_eta0(self):
        # Check whether expected ValueError on bad eta0
        self.factory(eta0=0, learning_rate="constant")

    @raises(ValueError)
    def test_sgd_bad_alpha(self):
        # Check whether expected ValueError on bad alpha
        self.factory(alpha=-.1)

    @raises(ValueError)
    def test_sgd_bad_penalty(self):
        # Check whether expected ValueError on bad penalty
        self.factory(penalty='foobar', l1_ratio=0.85)

    @raises(ValueError)
    def test_sgd_bad_loss(self):
        # Check whether expected ValueError on bad loss
        self.factory(loss="foobar")

    @raises(ValueError)
    def test_sgd_n_iter_param(self):
        # Test parameter validity check
        self.factory(n_iter=-10000)

    @raises(ValueError)
    def test_sgd_shuffle_param(self):
        # Test parameter validity check
        self.factory(shuffle="false")

    @raises(TypeError)
    def test_argument_coef(self):
        # Checks coef_init not allowed as model argument (only fit)
        # Provided coef_ does not match dataset.
        self.factory(coef_init=np.zeros((3,))).fit(X, Y)

    @raises(ValueError)
    def test_provide_coef(self):
        # Checks coef_init shape for the warm starts
        # Provided coef_ does not match dataset.
        self.factory().fit(X, Y, coef_init=np.zeros((3,)))

    @raises(ValueError)
    def test_set_intercept(self):
        # Checks intercept_ shape for the warm starts
        # Provided intercept_ does not match dataset.
        self.factory().fit(X, Y, intercept_init=np.zeros((3,)))

    def test_set_intercept_binary(self):
        # Checks intercept_ shape for the warm starts in binary case
        self.factory().fit(X5, Y5, intercept_init=0)

    def test_average_binary_computed_correctly(self):
        # Checks the SGDClassifier correctly computes the average weights
        eta = .1
        alpha = 2.
        n_samples = 20
        n_features = 10
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)

        clf = self.factory(loss='squared_loss',
                           learning_rate='constant',
                           eta0=eta, alpha=alpha,
                           fit_intercept=True,
                           n_iter=1, average=True, shuffle=False)

        # simple linear function without noise
        y = np.dot(X, w)
        y = np.sign(y)

        clf.fit(X, y)

        average_weights, average_intercept = self.asgd(X, y, eta, alpha)
        average_weights = average_weights.reshape(1, -1)
        assert_array_almost_equal(clf.coef_,
                                  average_weights,
                                  decimal=14)
        assert_almost_equal(clf.intercept_, average_intercept, decimal=14)

    def test_set_intercept_to_intercept(self):
        # Checks intercept_ shape consistency for the warm starts
        # Inconsistent intercept_ shape.
        clf = self.factory().fit(X5, Y5)
        self.factory().fit(X5, Y5, intercept_init=clf.intercept_)
        clf = self.factory().fit(X, Y)
        self.factory().fit(X, Y, intercept_init=clf.intercept_)

    @raises(ValueError)
    def test_sgd_at_least_two_labels(self):
        # Target must have at least two labels
        self.factory(alpha=0.01, n_iter=20).fit(X2, np.ones(9))

    def test_partial_fit_weight_class_balanced(self):
        # partial_fit with class_weight='balanced' not supported"""
        assert_raises_regexp(ValueError,
                             "class_weight 'balanced' is not supported for "
                             "partial_fit. In order to use 'balanced' weights, "
                             "use compute_class_weight\('balanced', classes, y\). "
                             "In place of y you can us a large enough sample "
                             "of the full training set target to properly "
                             "estimate the class frequency distributions. "
                             "Pass the resulting weights as the class_weight "
                             "parameter.",
                             self.factory(class_weight='balanced').partial_fit,
                             X, Y, classes=np.unique(Y))

    def test_sgd_multiclass(self):
        # Multi-class test case
        clf = self.factory(alpha=0.01, n_iter=20).fit(X2, Y2)
        assert_equal(clf.coef_.shape, (3, 2))
        assert_equal(clf.intercept_.shape, (3,))
        assert_equal(clf.decision_function([[0, 0]]).shape, (1, 3))
        pred = clf.predict(T2)
        assert_array_equal(pred, true_result2)

    def test_sgd_multiclass_average(self):
        eta = .001
        alpha = .01
        # Multi-class average test case
        clf = self.factory(loss='squared_loss',
                           learning_rate='constant',
                           eta0=eta, alpha=alpha,
                           fit_intercept=True,
                           n_iter=1, average=True, shuffle=False)

        np_Y2 = np.array(Y2)
        clf.fit(X2, np_Y2)
        classes = np.unique(np_Y2)

        for i, cl in enumerate(classes):
            y_i = np.ones(np_Y2.shape[0])
            y_i[np_Y2 != cl] = -1
            average_coef, average_intercept = self.asgd(X2, y_i, eta, alpha)
            assert_array_almost_equal(average_coef, clf.coef_[i], decimal=16)
            assert_almost_equal(average_intercept,
                                clf.intercept_[i],
                                decimal=16)

    def test_sgd_multiclass_with_init_coef(self):
        # Multi-class test case
        clf = self.factory(alpha=0.01, n_iter=20)
        clf.fit(X2, Y2, coef_init=np.zeros((3, 2)),
                intercept_init=np.zeros(3))
        assert_equal(clf.coef_.shape, (3, 2))
        assert_true(clf.intercept_.shape, (3,))
        pred = clf.predict(T2)
        assert_array_equal(pred, true_result2)

    def test_sgd_multiclass_njobs(self):
        # Multi-class test case with multi-core support
        clf = self.factory(alpha=0.01, n_iter=20, n_jobs=2).fit(X2, Y2)
        assert_equal(clf.coef_.shape, (3, 2))
        assert_equal(clf.intercept_.shape, (3,))
        assert_equal(clf.decision_function([[0, 0]]).shape, (1, 3))
        pred = clf.predict(T2)
        assert_array_equal(pred, true_result2)

    def test_set_coef_multiclass(self):
        # Checks coef_init and intercept_init shape for multi-class
        # problems
        # Provided coef_ does not match dataset
        clf = self.factory()
        assert_raises(ValueError, clf.fit, X2, Y2, coef_init=np.zeros((2, 2)))

        # Provided coef_ does match dataset
        clf = self.factory().fit(X2, Y2, coef_init=np.zeros((3, 2)))

        # Provided intercept_ does not match dataset
        clf = self.factory()
        assert_raises(ValueError, clf.fit, X2, Y2,
                      intercept_init=np.zeros((1,)))

        # Provided intercept_ does match dataset.
        clf = self.factory().fit(X2, Y2, intercept_init=np.zeros((3,)))

    def test_sgd_proba(self):
        # Check SGD.predict_proba

        # Hinge loss does not allow for conditional prob estimate.
        # We cannot use the factory here, because it defines predict_proba
        # anyway.
        clf = SGDClassifier(loss="hinge", alpha=0.01, n_iter=10).fit(X, Y)
        assert_false(hasattr(clf, "predict_proba"))
        assert_false(hasattr(clf, "predict_log_proba"))

        # log and modified_huber losses can output probability estimates
        # binary case
        for loss in ["log", "modified_huber"]:
            clf = self.factory(loss=loss, alpha=0.01, n_iter=10)
            clf.fit(X, Y)
            p = clf.predict_proba([[3, 2]])
            assert_true(p[0, 1] > 0.5)
            p = clf.predict_proba([[-1, -1]])
            assert_true(p[0, 1] < 0.5)

            p = clf.predict_log_proba([[3, 2]])
            assert_true(p[0, 1] > p[0, 0])
            p = clf.predict_log_proba([[-1, -1]])
            assert_true(p[0, 1] < p[0, 0])

        # log loss multiclass probability estimates
        clf = self.factory(loss="log", alpha=0.01, n_iter=10).fit(X2, Y2)

        d = clf.decision_function([[.1, -.1], [.3, .2]])
        p = clf.predict_proba([[.1, -.1], [.3, .2]])
        assert_array_equal(np.argmax(p, axis=1), np.argmax(d, axis=1))
        assert_almost_equal(p[0].sum(), 1)
        assert_true(np.all(p[0] >= 0))

        p = clf.predict_proba([[-1, -1]])
        d = clf.decision_function([[-1, -1]])
        assert_array_equal(np.argsort(p[0]), np.argsort(d[0]))

        l = clf.predict_log_proba([[3, 2]])
        p = clf.predict_proba([[3, 2]])
        assert_array_almost_equal(np.log(p), l)

        l = clf.predict_log_proba([[-1, -1]])
        p = clf.predict_proba([[-1, -1]])
        assert_array_almost_equal(np.log(p), l)

        # Modified Huber multiclass probability estimates; requires a separate
        # test because the hard zero/one probabilities may destroy the
        # ordering present in decision_function output.
        clf = self.factory(loss="modified_huber", alpha=0.01, n_iter=10)
        clf.fit(X2, Y2)
        d = clf.decision_function([[3, 2]])
        p = clf.predict_proba([[3, 2]])
        if not isinstance(self, SparseSGDClassifierTestCase):
            assert_equal(np.argmax(d, axis=1), np.argmax(p, axis=1))
        else:   # XXX the sparse test gets a different X2 (?)
            assert_equal(np.argmin(d, axis=1), np.argmin(p, axis=1))

        # the following sample produces decision_function values < -1,
        # which would cause naive normalization to fail (see comment
        # in SGDClassifier.predict_proba)
        x = X.mean(axis=0)
        d = clf.decision_function([x])
        if np.all(d < -1):  # XXX not true in sparse test case (why?)
            p = clf.predict_proba([x])
            assert_array_almost_equal(p[0], [1 / 3.] * 3)

    def test_sgd_l1(self):
        # Test L1 regularization
        n = len(X4)
        rng = np.random.RandomState(13)
        idx = np.arange(n)
        rng.shuffle(idx)

        X = X4[idx, :]
        Y = Y4[idx]

        clf = self.factory(penalty='l1', alpha=.2, fit_intercept=False,
                           n_iter=2000, shuffle=False)
        clf.fit(X, Y)
        assert_array_equal(clf.coef_[0, 1:-1], np.zeros((4,)))
        pred = clf.predict(X)
        assert_array_equal(pred, Y)

        # test sparsify with dense inputs
        clf.sparsify()
        assert_true(sp.issparse(clf.coef_))
        pred = clf.predict(X)
        assert_array_equal(pred, Y)

        # pickle and unpickle with sparse coef_
        clf = pickle.loads(pickle.dumps(clf))
        assert_true(sp.issparse(clf.coef_))
        pred = clf.predict(X)
        assert_array_equal(pred, Y)

    def test_class_weights(self):
        # Test class weights.
        X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                      [1.0, 1.0], [1.0, 0.0]])
        y = [1, 1, 1, -1, -1]

        clf = self.factory(alpha=0.1, n_iter=1000, fit_intercept=False,
                           class_weight=None)
        clf.fit(X, y)
        assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([1]))

        # we give a small weights to class 1
        clf = self.factory(alpha=0.1, n_iter=1000, fit_intercept=False,
                           class_weight={1: 0.001})
        clf.fit(X, y)

        # now the hyperplane should rotate clock-wise and
        # the prediction on this point should shift
        assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([-1]))

    def test_equal_class_weight(self):
        # Test if equal class weights approx. equals no class weights.
        X = [[1, 0], [1, 0], [0, 1], [0, 1]]
        y = [0, 0, 1, 1]
        clf = self.factory(alpha=0.1, n_iter=1000, class_weight=None)
        clf.fit(X, y)

        X = [[1, 0], [0, 1]]
        y = [0, 1]
        clf_weighted = self.factory(alpha=0.1, n_iter=1000,
                                    class_weight={0: 0.5, 1: 0.5})
        clf_weighted.fit(X, y)

        # should be similar up to some epsilon due to learning rate schedule
        assert_almost_equal(clf.coef_, clf_weighted.coef_, decimal=2)

    @raises(ValueError)
    def test_wrong_class_weight_label(self):
        # ValueError due to not existing class label.
        clf = self.factory(alpha=0.1, n_iter=1000, class_weight={0: 0.5})
        clf.fit(X, Y)

    @raises(ValueError)
    def test_wrong_class_weight_format(self):
        # ValueError due to wrong class_weight argument type.
        clf = self.factory(alpha=0.1, n_iter=1000, class_weight=[0.5])
        clf.fit(X, Y)

    def test_weights_multiplied(self):
        # Tests that class_weight and sample_weight are multiplicative
        class_weights = {1: .6, 2: .3}
        rng = np.random.RandomState(0)
        sample_weights = rng.random_sample(Y4.shape[0])
        multiplied_together = np.copy(sample_weights)
        multiplied_together[Y4 == 1] *= class_weights[1]
        multiplied_together[Y4 == 2] *= class_weights[2]

        clf1 = self.factory(alpha=0.1, n_iter=20, class_weight=class_weights)
        clf2 = self.factory(alpha=0.1, n_iter=20)

        clf1.fit(X4, Y4, sample_weight=sample_weights)
        clf2.fit(X4, Y4, sample_weight=multiplied_together)

        assert_almost_equal(clf1.coef_, clf2.coef_)

    def test_balanced_weight(self):
        # Test class weights for imbalanced data"""
        # compute reference metrics on iris dataset that is quite balanced by
        # default
        X, y = iris.data, iris.target
        X = scale(X)
        idx = np.arange(X.shape[0])
        rng = np.random.RandomState(6)
        rng.shuffle(idx)
        X = X[idx]
        y = y[idx]
        clf = self.factory(alpha=0.0001, n_iter=1000,
                           class_weight=None, shuffle=False).fit(X, y)
        assert_almost_equal(metrics.f1_score(y, clf.predict(X), average='weighted'), 0.96,
                            decimal=1)

        # make the same prediction using balanced class_weight
        clf_balanced = self.factory(alpha=0.0001, n_iter=1000,
                                    class_weight="balanced",
                                    shuffle=False).fit(X, y)
        assert_almost_equal(metrics.f1_score(y, clf_balanced.predict(X), average='weighted'), 0.96,
                            decimal=1)

        # Make sure that in the balanced case it does not change anything
        # to use "balanced"
        assert_array_almost_equal(clf.coef_, clf_balanced.coef_, 6)

        # build an very very imbalanced dataset out of iris data
        X_0 = X[y == 0, :]
        y_0 = y[y == 0]

        X_imbalanced = np.vstack([X] + [X_0] * 10)
        y_imbalanced = np.concatenate([y] + [y_0] * 10)

        # fit a model on the imbalanced data without class weight info
        clf = self.factory(n_iter=1000, class_weight=None, shuffle=False)
        clf.fit(X_imbalanced, y_imbalanced)
        y_pred = clf.predict(X)
        assert_less(metrics.f1_score(y, y_pred, average='weighted'), 0.96)

        # fit a model with balanced class_weight enabled
        clf = self.factory(n_iter=1000, class_weight="balanced", shuffle=False)
        clf.fit(X_imbalanced, y_imbalanced)
        y_pred = clf.predict(X)
        assert_greater(metrics.f1_score(y, y_pred, average='weighted'), 0.96)

        # fit another using a fit parameter override
        clf = self.factory(n_iter=1000, class_weight="balanced", shuffle=False)
        clf.fit(X_imbalanced, y_imbalanced)
        y_pred = clf.predict(X)
        assert_greater(metrics.f1_score(y, y_pred, average='weighted'), 0.96)

    def test_sample_weights(self):
        # Test weights on individual samples
        X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                      [1.0, 1.0], [1.0, 0.0]])
        y = [1, 1, 1, -1, -1]

        clf = self.factory(alpha=0.1, n_iter=1000, fit_intercept=False)
        clf.fit(X, y)
        assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([1]))

        # we give a small weights to class 1
        clf.fit(X, y, sample_weight=[0.001] * 3 + [1] * 2)

        # now the hyperplane should rotate clock-wise and
        # the prediction on this point should shift
        assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([-1]))

    @raises(ValueError)
    def test_wrong_sample_weights(self):
        # Test if ValueError is raised if sample_weight has wrong shape
        clf = self.factory(alpha=0.1, n_iter=1000, fit_intercept=False)
        # provided sample_weight too long
        clf.fit(X, Y, sample_weight=np.arange(7))

    @raises(ValueError)
    def test_partial_fit_exception(self):
        clf = self.factory(alpha=0.01)
        # classes was not specified
        clf.partial_fit(X3, Y3)

    def test_partial_fit_binary(self):
        third = X.shape[0] // 3
        clf = self.factory(alpha=0.01)
        classes = np.unique(Y)

        clf.partial_fit(X[:third], Y[:third], classes=classes)
        assert_equal(clf.coef_.shape, (1, X.shape[1]))
        assert_equal(clf.intercept_.shape, (1,))
        assert_equal(clf.decision_function([[0, 0]]).shape, (1, ))
        id1 = id(clf.coef_.data)

        clf.partial_fit(X[third:], Y[third:])
        id2 = id(clf.coef_.data)
        # check that coef_ haven't been re-allocated
        assert_true(id1, id2)

        y_pred = clf.predict(T)
        assert_array_equal(y_pred, true_result)

    def test_partial_fit_multiclass(self):
        third = X2.shape[0] // 3
        clf = self.factory(alpha=0.01)
        classes = np.unique(Y2)

        clf.partial_fit(X2[:third], Y2[:third], classes=classes)
        assert_equal(clf.coef_.shape, (3, X2.shape[1]))
        assert_equal(clf.intercept_.shape, (3,))
        assert_equal(clf.decision_function([[0, 0]]).shape, (1, 3))
        id1 = id(clf.coef_.data)

        clf.partial_fit(X2[third:], Y2[third:])
        id2 = id(clf.coef_.data)
        # check that coef_ haven't been re-allocated
        assert_true(id1, id2)

    def test_partial_fit_multiclass_average(self):
        third = X2.shape[0] // 3
        clf = self.factory(alpha=0.01, average=X2.shape[0])
        classes = np.unique(Y2)

        clf.partial_fit(X2[:third], Y2[:third], classes=classes)
        assert_equal(clf.coef_.shape, (3, X2.shape[1]))
        assert_equal(clf.intercept_.shape, (3,))

        clf.partial_fit(X2[third:], Y2[third:])
        assert_equal(clf.coef_.shape, (3, X2.shape[1]))
        assert_equal(clf.intercept_.shape, (3,))

    def test_fit_then_partial_fit(self):
        # Partial_fit should work after initial fit in the multiclass case.
        # Non-regression test for #2496; fit would previously produce a
        # Fortran-ordered coef_ that subsequent partial_fit couldn't handle.
        clf = self.factory()
        clf.fit(X2, Y2)
        clf.partial_fit(X2, Y2)     # no exception here

    def _test_partial_fit_equal_fit(self, lr):
        for X_, Y_, T_ in ((X, Y, T), (X2, Y2, T2)):
            clf = self.factory(alpha=0.01, eta0=0.01, n_iter=2,
                               learning_rate=lr, shuffle=False)
            clf.fit(X_, Y_)
            y_pred = clf.decision_function(T_)
            t = clf.t_

            classes = np.unique(Y_)
            clf = self.factory(alpha=0.01, eta0=0.01, learning_rate=lr,
                               shuffle=False)
            for i in range(2):
                clf.partial_fit(X_, Y_, classes=classes)
            y_pred2 = clf.decision_function(T_)

            assert_equal(clf.t_, t)
            assert_array_almost_equal(y_pred, y_pred2, decimal=2)

    def test_partial_fit_equal_fit_constant(self):
        self._test_partial_fit_equal_fit("constant")

    def test_partial_fit_equal_fit_optimal(self):
        self._test_partial_fit_equal_fit("optimal")

    def test_partial_fit_equal_fit_invscaling(self):
        self._test_partial_fit_equal_fit("invscaling")

    def test_regression_losses(self):
        clf = self.factory(alpha=0.01, learning_rate="constant",
                           eta0=0.1, loss="epsilon_insensitive")
        clf.fit(X, Y)
        assert_equal(1.0, np.mean(clf.predict(X) == Y))

        clf = self.factory(alpha=0.01, learning_rate="constant",
                           eta0=0.1, loss="squared_epsilon_insensitive")
        clf.fit(X, Y)
        assert_equal(1.0, np.mean(clf.predict(X) == Y))

        clf = self.factory(alpha=0.01, loss="huber")
        clf.fit(X, Y)
        assert_equal(1.0, np.mean(clf.predict(X) == Y))

        clf = self.factory(alpha=0.01, learning_rate="constant", eta0=0.01,
                           loss="squared_loss")
        clf.fit(X, Y)
        assert_equal(1.0, np.mean(clf.predict(X) == Y))

    def test_warm_start_multiclass(self):
        self._test_warm_start(X2, Y2, "optimal")

    def test_multiple_fit(self):
        # Test multiple calls of fit w/ different shaped inputs.
        clf = self.factory(alpha=0.01, n_iter=5,
                           shuffle=False)
        clf.fit(X, Y)
        assert_true(hasattr(clf, "coef_"))

        # Non-regression test: try fitting with a different label set.
        y = [["ham", "spam"][i] for i in LabelEncoder().fit_transform(Y)]
        clf.fit(X[:, :-1], y)


class SparseSGDClassifierTestCase(DenseSGDClassifierTestCase):
    """Run exactly the same tests using the sparse representation variant"""

    factory_class = SparseSGDClassifier


###############################################################################
# Regression Test Case

class DenseSGDRegressorTestCase(unittest.TestCase, CommonTest):
    """Test suite for the dense representation variant of SGD"""

    factory_class = SGDRegressor

    def test_sgd(self):
        # Check that SGD gives any results.
        clf = self.factory(alpha=0.1, n_iter=2,
                           fit_intercept=False)
        clf.fit([[0, 0], [1, 1], [2, 2]], [0, 1, 2])
        assert_equal(clf.coef_[0], clf.coef_[1])

    @raises(ValueError)
    def test_sgd_bad_penalty(self):
        # Check whether expected ValueError on bad penalty
        self.factory(penalty='foobar', l1_ratio=0.85)

    @raises(ValueError)
    def test_sgd_bad_loss(self):
        # Check whether expected ValueError on bad loss
        self.factory(loss="foobar")

    def test_sgd_averaged_computed_correctly(self):
        # Tests the average regressor matches the naive implementation

        eta = .001
        alpha = .01
        n_samples = 20
        n_features = 10
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)

        # simple linear function without noise
        y = np.dot(X, w)

        clf = self.factory(loss='squared_loss',
                           learning_rate='constant',
                           eta0=eta, alpha=alpha,
                           fit_intercept=True,
                           n_iter=1, average=True, shuffle=False)

        clf.fit(X, y)
        average_weights, average_intercept = self.asgd(X, y, eta, alpha)

        assert_array_almost_equal(clf.coef_,
                                  average_weights,
                                  decimal=16)
        assert_almost_equal(clf.intercept_, average_intercept, decimal=16)

    def test_sgd_averaged_partial_fit(self):
        # Tests whether the partial fit yields the same average as the fit
        eta = .001
        alpha = .01
        n_samples = 20
        n_features = 10
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)

        # simple linear function without noise
        y = np.dot(X, w)

        clf = self.factory(loss='squared_loss',
                           learning_rate='constant',
                           eta0=eta, alpha=alpha,
                           fit_intercept=True,
                           n_iter=1, average=True, shuffle=False)

        clf.partial_fit(X[:int(n_samples / 2)][:], y[:int(n_samples / 2)])
        clf.partial_fit(X[int(n_samples / 2):][:], y[int(n_samples / 2):])
        average_weights, average_intercept = self.asgd(X, y, eta, alpha)

        assert_array_almost_equal(clf.coef_,
                                  average_weights,
                                  decimal=16)
        assert_almost_equal(clf.intercept_[0], average_intercept, decimal=16)

    def test_average_sparse(self):
        # Checks the average weights on data with 0s

        eta = .001
        alpha = .01
        clf = self.factory(loss='squared_loss',
                           learning_rate='constant',
                           eta0=eta, alpha=alpha,
                           fit_intercept=True,
                           n_iter=1, average=True, shuffle=False)

        n_samples = Y3.shape[0]

        clf.partial_fit(X3[:int(n_samples / 2)][:], Y3[:int(n_samples / 2)])
        clf.partial_fit(X3[int(n_samples / 2):][:], Y3[int(n_samples / 2):])
        average_weights, average_intercept = self.asgd(X3, Y3, eta, alpha)

        assert_array_almost_equal(clf.coef_,
                                  average_weights,
                                  decimal=16)
        assert_almost_equal(clf.intercept_, average_intercept, decimal=16)

    def test_sgd_least_squares_fit(self):
        xmin, xmax = -5, 5
        n_samples = 100
        rng = np.random.RandomState(0)
        X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

        # simple linear function without noise
        y = 0.5 * X.ravel()

        clf = self.factory(loss='squared_loss', alpha=0.1, n_iter=20,
                           fit_intercept=False)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_greater(score, 0.99)

        # simple linear function with noise
        y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

        clf = self.factory(loss='squared_loss', alpha=0.1, n_iter=20,
                           fit_intercept=False)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_greater(score, 0.5)

    def test_sgd_epsilon_insensitive(self):
        xmin, xmax = -5, 5
        n_samples = 100
        rng = np.random.RandomState(0)
        X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

        # simple linear function without noise
        y = 0.5 * X.ravel()

        clf = self.factory(loss='epsilon_insensitive', epsilon=0.01,
                           alpha=0.1, n_iter=20,
                           fit_intercept=False)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_true(score > 0.99)

        # simple linear function with noise
        y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

        clf = self.factory(loss='epsilon_insensitive', epsilon=0.01,
                           alpha=0.1, n_iter=20,
                           fit_intercept=False)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_true(score > 0.5)

    def test_sgd_huber_fit(self):
        xmin, xmax = -5, 5
        n_samples = 100
        rng = np.random.RandomState(0)
        X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

        # simple linear function without noise
        y = 0.5 * X.ravel()

        clf = self.factory(loss="huber", epsilon=0.1, alpha=0.1, n_iter=20,
                           fit_intercept=False)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_greater(score, 0.99)

        # simple linear function with noise
        y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

        clf = self.factory(loss="huber", epsilon=0.1, alpha=0.1, n_iter=20,
                           fit_intercept=False)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_greater(score, 0.5)

    def test_elasticnet_convergence(self):
        # Check that the SGD output is consistent with coordinate descent

        n_samples, n_features = 1000, 5
        rng = np.random.RandomState(0)
        X = rng.randn(n_samples, n_features)
        # ground_truth linear model that generate y from X and to which the
        # models should converge if the regularizer would be set to 0.0
        ground_truth_coef = rng.randn(n_features)
        y = np.dot(X, ground_truth_coef)

        # XXX: alpha = 0.1 seems to cause convergence problems
        for alpha in [0.01, 0.001]:
            for l1_ratio in [0.5, 0.8, 1.0]:
                cd = linear_model.ElasticNet(alpha=alpha, l1_ratio=l1_ratio,
                                             fit_intercept=False)
                cd.fit(X, y)
                sgd = self.factory(penalty='elasticnet', n_iter=50,
                                   alpha=alpha, l1_ratio=l1_ratio,
                                   fit_intercept=False)
                sgd.fit(X, y)
                err_msg = ("cd and sgd did not converge to comparable "
                           "results for alpha=%f and l1_ratio=%f"
                           % (alpha, l1_ratio))
                assert_almost_equal(cd.coef_, sgd.coef_, decimal=2,
                                    err_msg=err_msg)

    @ignore_warnings
    def test_partial_fit(self):
        third = X.shape[0] // 3
        clf = self.factory(alpha=0.01)

        clf.partial_fit(X[:third], Y[:third])
        assert_equal(clf.coef_.shape, (X.shape[1], ))
        assert_equal(clf.intercept_.shape, (1,))
        assert_equal(clf.predict([[0, 0]]).shape, (1, ))
        id1 = id(clf.coef_.data)

        clf.partial_fit(X[third:], Y[third:])
        id2 = id(clf.coef_.data)
        # check that coef_ haven't been re-allocated
        assert_true(id1, id2)

    def _test_partial_fit_equal_fit(self, lr):
        clf = self.factory(alpha=0.01, n_iter=2, eta0=0.01,
                           learning_rate=lr, shuffle=False)
        clf.fit(X, Y)
        y_pred = clf.predict(T)
        t = clf.t_

        clf = self.factory(alpha=0.01, eta0=0.01,
                           learning_rate=lr, shuffle=False)
        for i in range(2):
            clf.partial_fit(X, Y)
        y_pred2 = clf.predict(T)

        assert_equal(clf.t_, t)
        assert_array_almost_equal(y_pred, y_pred2, decimal=2)

    def test_partial_fit_equal_fit_constant(self):
        self._test_partial_fit_equal_fit("constant")

    def test_partial_fit_equal_fit_optimal(self):
        self._test_partial_fit_equal_fit("optimal")

    def test_partial_fit_equal_fit_invscaling(self):
        self._test_partial_fit_equal_fit("invscaling")

    def test_loss_function_epsilon(self):
        clf = self.factory(epsilon=0.9)
        clf.set_params(epsilon=0.1)
        assert clf.loss_functions['huber'][1] == 0.1


class SparseSGDRegressorTestCase(DenseSGDRegressorTestCase):
    # Run exactly the same tests using the sparse representation variant

    factory_class = SparseSGDRegressor


def test_l1_ratio():
    # Test if l1 ratio extremes match L1 and L2 penalty settings.
    X, y = datasets.make_classification(n_samples=1000,
                                        n_features=100, n_informative=20,
                                        random_state=1234)

    # test if elasticnet with l1_ratio near 1 gives same result as pure l1
    est_en = SGDClassifier(alpha=0.001, penalty='elasticnet',
                           l1_ratio=0.9999999999, random_state=42).fit(X, y)
    est_l1 = SGDClassifier(alpha=0.001, penalty='l1', random_state=42).fit(X, y)
    assert_array_almost_equal(est_en.coef_, est_l1.coef_)

    # test if elasticnet with l1_ratio near 0 gives same result as pure l2
    est_en = SGDClassifier(alpha=0.001, penalty='elasticnet',
                           l1_ratio=0.0000000001, random_state=42).fit(X, y)
    est_l2 = SGDClassifier(alpha=0.001, penalty='l2', random_state=42).fit(X, y)
    assert_array_almost_equal(est_en.coef_, est_l2.coef_)


def test_underflow_or_overlow():
    with np.errstate(all='raise'):
        # Generate some weird data with hugely unscaled features
        rng = np.random.RandomState(0)
        n_samples = 100
        n_features = 10

        X = rng.normal(size=(n_samples, n_features))
        X[:, :2] *= 1e300
        assert_true(np.isfinite(X).all())

        # Use MinMaxScaler to scale the data without introducing a numerical
        # instability (computing the standard deviation naively is not possible
        # on this data)
        X_scaled = MinMaxScaler().fit_transform(X)
        assert_true(np.isfinite(X_scaled).all())

        # Define a ground truth on the scaled data
        ground_truth = rng.normal(size=n_features)
        y = (np.dot(X_scaled, ground_truth) > 0.).astype(np.int32)
        assert_array_equal(np.unique(y), [0, 1])

        model = SGDClassifier(alpha=0.1, loss='squared_hinge', n_iter=500)

        # smoke test: model is stable on scaled data
        model.fit(X_scaled, y)
        assert_true(np.isfinite(model.coef_).all())

        # model is numerically unstable on unscaled data
        msg_regxp = (r"Floating-point under-/overflow occurred at epoch #.*"
                     " Scaling input data with StandardScaler or MinMaxScaler"
                     " might help.")
        assert_raises_regexp(ValueError, msg_regxp, model.fit, X, y)


def test_numerical_stability_large_gradient():
    # Non regression test case for numerical stability on scaled problems
    # where the gradient can still explode with some losses
    model = SGDClassifier(loss='squared_hinge', n_iter=10, shuffle=True,
                          penalty='elasticnet', l1_ratio=0.3, alpha=0.01,
                          eta0=0.001, random_state=0)
    with np.errstate(all='raise'):
        model.fit(iris.data, iris.target)
    assert_true(np.isfinite(model.coef_).all())


def test_large_regularization():
    # Non regression tests for numerical stability issues caused by large
    # regularization parameters
    for penalty in ['l2', 'l1', 'elasticnet']:
        model = SGDClassifier(alpha=1e5, learning_rate='constant', eta0=0.1,
                              n_iter=5, penalty=penalty, shuffle=False)
        with np.errstate(all='raise'):
            model.fit(iris.data, iris.target)
        assert_array_almost_equal(model.coef_, np.zeros_like(model.coef_))


def _test_gradient_common(loss_function, cases):
    # Test gradient of different loss functions
    # cases is a list of (p, y, expected)
    for p, y, expected in cases:
        assert_almost_equal(loss_function.dloss(p, y), expected)


def test_gradient_hinge():
    # Test Hinge (hinge / perceptron)
    # hinge
    loss = sgd_fast.Hinge(1.0)
    cases = [
        # (p, y, expected)
        (1.1, 1.0, 0.0), (-2.0, -1.0, 0.0),
        (1.0, 1.0, -1.0), (-1.0, -1.0, 1.0), (0.5, 1.0, -1.0),
        (2.0, -1.0, 1.0), (-0.5, -1.0, 1.0), (0.0, 1.0, -1.0)
    ]
    _test_gradient_common(loss, cases)

    # perceptron
    loss = sgd_fast.Hinge(0.0)
    cases = [
        # (p, y, expected)
        (1.0, 1.0, 0.0), (-0.1, -1.0, 0.0),
        (0.0, 1.0, -1.0), (0.0, -1.0, 1.0), (0.5, -1.0, 1.0),
        (2.0, -1.0, 1.0), (-0.5, 1.0, -1.0), (-1.0, 1.0, -1.0),
    ]
    _test_gradient_common(loss, cases)


def test_gradient_squared_hinge():
    # Test SquaredHinge
    loss = sgd_fast.SquaredHinge(1.0)
    cases = [
        # (p, y, expected)
        (1.0, 1.0, 0.0), (-2.0, -1.0, 0.0), (1.0, -1.0, 4.0),
        (-1.0, 1.0, -4.0), (0.5, 1.0, -1.0), (0.5, -1.0, 3.0)
    ]
    _test_gradient_common(loss, cases)


def test_gradient_log():
    # Test Log (logistic loss)
    loss = sgd_fast.Log()
    cases = [
        # (p, y, expected)
        (1.0, 1.0, -1.0 / (np.exp(1.0) + 1.0)),
        (1.0, -1.0, 1.0 / (np.exp(-1.0) + 1.0)),
        (-1.0, -1.0, 1.0 / (np.exp(1.0) + 1.0)),
        (-1.0, 1.0, -1.0 / (np.exp(-1.0) + 1.0)),
        (0.0, 1.0, -0.5), (0.0, -1.0, 0.5),
        (17.9, -1.0, 1.0), (-17.9, 1.0, -1.0),
    ]
    _test_gradient_common(loss, cases)
    assert_almost_equal(loss.dloss(18.1, 1.0), np.exp(-18.1) * -1.0, 16)
    assert_almost_equal(loss.dloss(-18.1, -1.0), np.exp(-18.1) * 1.0, 16)


def test_gradient_squared_loss():
    # Test SquaredLoss
    loss = sgd_fast.SquaredLoss()
    cases = [
        # (p, y, expected)
        (0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 1.0),
        (0.5, -1.0, 1.5), (-2.5, 2.0, -4.5)
    ]
    _test_gradient_common(loss, cases)


def test_gradient_huber():
    # Test Huber
    loss = sgd_fast.Huber(0.1)
    cases = [
        # (p, y, expected)
        (0.0, 0.0, 0.0), (0.1, 0.0, 0.1), (0.0, 0.1, -0.1),
        (3.95, 4.0, -0.05), (5.0, 2.0, 0.1), (-1.0, 5.0, -0.1)
    ]
    _test_gradient_common(loss, cases)


def test_gradient_modified_huber():
    # Test ModifiedHuber
    loss = sgd_fast.ModifiedHuber()
    cases = [
        # (p, y, expected)
        (1.0, 1.0, 0.0), (-1.0, -1.0, 0.0), (2.0, 1.0, 0.0),
        (0.0, 1.0, -2.0), (-1.0, 1.0, -4.0), (0.5, -1.0, 3.0),
        (0.5, -1.0, 3.0), (-2.0, 1.0, -4.0), (-3.0, 1.0, -4.0)
    ]
    _test_gradient_common(loss, cases)


def test_gradient_epsilon_insensitive():
    # Test EpsilonInsensitive
    loss = sgd_fast.EpsilonInsensitive(0.1)
    cases = [
        (0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (-2.05, -2.0, 0.0),
        (3.05, 3.0, 0.0), (2.2, 2.0, 1.0), (2.0, -1.0, 1.0),
        (2.0, 2.2, -1.0), (-2.0, 1.0, -1.0)
    ]
    _test_gradient_common(loss, cases)


def test_gradient_squared_epsilon_insensitive():
    # Test SquaredEpsilonInsensitive
    loss = sgd_fast.SquaredEpsilonInsensitive(0.1)
    cases = [
        (0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (-2.05, -2.0, 0.0),
        (3.05, 3.0, 0.0), (2.2, 2.0, 0.2), (2.0, -1.0, 5.8),
        (2.0, 2.2, -0.2), (-2.0, 1.0, -5.8)
    ]
    _test_gradient_common(loss, cases)
