import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_equal, assert_approx_equal
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from sklearn import linear_model, datasets, metrics
from sklearn import preprocessing
from sklearn.linear_model import SGDClassifier, SGDRegressor
from sklearn.utils.testing import assert_greater, assert_less
from sklearn.base import clone

import unittest
from nose.tools import raises
from nose.tools import assert_raises, assert_true, assert_equal


class SparseSGDClassifier(SGDClassifier):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDClassifier.fit(self, X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDClassifier.partial_fit(self, X, y, *args, **kw)

    def decision_function(self, X, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDClassifier.decision_function(self, X, *args, **kw)

    def predict_proba(self, X, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDClassifier.predict_proba(self, X, *args, **kw)

    def predict_log_proba(self, X, *args, **kw):
        X = sp.csr_matrix(X)
        return SGDClassifier.predict_log_proba(self, X, *args, **kw)


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

##
## Test Data
##


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

# test sample 4 - two more or less redundent feature groups
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

##
## Classification Test Case
##


class CommonTest(object):

    def _test_warm_start(self, lr):
        # Test that explicit warm restart...
        clf = self.factory(alpha=0.01, eta0=0.01, n_iter=5, shuffle=False,
                           learning_rate=lr)
        clf.fit(X, Y)

        clf2 = self.factory(alpha=0.001, eta0=0.01, n_iter=5, shuffle=False,
                            learning_rate=lr)
        clf2.fit(X, Y,
                 coef_init=clf.coef_.copy(),
                 intercept_init=clf.intercept_.copy())

        #... and implicit warm restart are equivalent.
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
        self._test_warm_start("constant")

    def test_warm_start_invscaling(self):
        self._test_warm_start("invscaling")

    def test_warm_start_optimal(self):
        self._test_warm_start("optimal")

    def test_multiple_fit(self):
        """Test multiple calls of fit w/ different shaped inputs."""
        clf = self.factory(alpha=0.01, n_iter=5,
                           shuffle=False)
        clf.fit(X, Y)
        assert_true(hasattr(clf, "coef_"))

        clf.fit(X[:, :-1], Y)

    def test_input_format(self):
        """Input format tests. """
        clf = self.factory(alpha=0.01, n_iter=5,
                           shuffle=False)
        Y_ = np.array(Y)[:, np.newaxis]
        clf.fit(X, Y_)

        Y_ = np.c_[Y_, Y_]
        assert_raises(ValueError, clf.fit, X, Y_)

    def test_clone(self):
        """Test whether clone works ok. """
        clf = self.factory(alpha=0.01, n_iter=5, penalty='l1')
        clf = clone(clf)
        clf.set_params(penalty='l2')
        clf.fit(X, Y)

        clf2 = self.factory(alpha=0.01, n_iter=5, penalty='l2')
        clf2.fit(X, Y)

        assert_array_equal(clf.coef_, clf2.coef_)


class DenseSGDClassifierTestCase(unittest.TestCase, CommonTest):
    """Test suite for the dense representation variant of SGD"""

    factory = SGDClassifier

    def test_sgd(self):
        """Check that SGD gives any results :-)"""

        clf = self.factory(penalty='l2', alpha=0.01, fit_intercept=True,
                      n_iter=10, shuffle=True)
        clf.fit(X, Y)
        #assert_almost_equal(clf.coef_[0], clf.coef_[1], decimal=7)
        assert_array_equal(clf.predict(T), true_result)

    @raises(ValueError)
    def test_sgd_bad_penalty(self):
        """Check whether expected ValueError on bad penalty"""
        self.factory(penalty='foobar', rho=0.85)

    @raises(ValueError)
    def test_sgd_bad_loss(self):
        """Check whether expected ValueError on bad loss"""
        self.factory(loss="foobar")

    @raises(ValueError)
    def test_sgd_n_iter_param(self):
        """Test parameter validity check"""
        self.factory(n_iter=-10000)

    @raises(ValueError)
    def test_sgd_shuffle_param(self):
        """Test parameter validity check"""
        self.factory(shuffle="false")

    @raises(TypeError)
    def test_arument_coef(self):
        """Checks coef_init not allowed as model argument (only fit)"""
        # Provided coef_ does not match dataset.
        self.factory(coef_init=np.zeros((3,))).fit(X, Y)

    @raises(ValueError)
    def test_provide_coef(self):
        """Checks coef_init shape for the warm starts"""
        # Provided coef_ does not match dataset.
        self.factory().fit(X, Y, coef_init=np.zeros((3,)))

    @raises(ValueError)
    def test_set_intercept(self):
        """Checks intercept_ shape for the warm starts"""
        # Provided intercept_ does not match dataset.
        self.factory().fit(X, Y, intercept_init=np.zeros((3,)))

    def test_set_intercept_binary(self):
        """Checks intercept_ shape for the warm starts in binary case"""
        self.factory().fit(X5, Y5, intercept_init=0)

    def test_set_intercept_to_intercept(self):
        """Checks intercept_ shape consistency for the warm starts"""
        # Inconsistent intercept_ shape.
        clf = self.factory().fit(X5, Y5)
        self.factory().fit(X5, Y5, intercept_init=clf.intercept_)
        clf = self.factory().fit(X, Y)
        self.factory().fit(X, Y, intercept_init=clf.intercept_)

    @raises(ValueError)
    def test_sgd_at_least_two_labels(self):
        """Target must have at least two labels"""
        self.factory(alpha=0.01, n_iter=20).fit(X2, np.ones(9))

    def test_sgd_multiclass(self):
        """Multi-class test case"""
        clf = self.factory(alpha=0.01, n_iter=20).fit(X2, Y2)
        assert_equal(clf.coef_.shape, (3, 2))
        assert_equal(clf.intercept_.shape, (3,))
        assert_equal(clf.decision_function([0, 0]).shape, (1, 3))
        pred = clf.predict(T2)
        assert_array_equal(pred, true_result2)

    def test_sgd_multiclass_with_init_coef(self):
        """Multi-class test case"""
        clf = self.factory(alpha=0.01, n_iter=20)
        clf.fit(X2, Y2, coef_init=np.zeros((3, 2)),
                intercept_init=np.zeros(3))
        assert_equal(clf.coef_.shape, (3, 2))
        assert_true(clf.intercept_.shape, (3,))
        pred = clf.predict(T2)
        assert_array_equal(pred, true_result2)

    def test_sgd_multiclass_njobs(self):
        """Multi-class test case with multi-core support"""
        clf = self.factory(alpha=0.01, n_iter=20, n_jobs=2).fit(X2, Y2)
        assert_equal(clf.coef_.shape, (3, 2))
        assert_equal(clf.intercept_.shape, (3,))
        assert_equal(clf.decision_function([0, 0]).shape, (1, 3))
        pred = clf.predict(T2)
        assert_array_equal(pred, true_result2)

    def test_set_coef_multiclass(self):
        """Checks coef_init and intercept_init shape for for multi-class
        problems"""
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
        """Check SGD.predict_proba"""

        # hinge loss does not allow for conditional prob estimate
        clf = self.factory(loss="hinge", alpha=0.01, n_iter=10).fit(X, Y)
        assert_raises(NotImplementedError, clf.predict_proba, [3, 2])

        # the log and modified_huber losses can output "probability" estimates
        for loss in ("log", "modified_huber"):
            clf = self.factory(loss=loss, alpha=0.01, n_iter=10).fit(X, Y)
            p = clf.predict_proba([3, 2])
            assert_true(p[0, 1] > 0.5)
            p = clf.predict_proba([-1, -1])
            assert_true(p[0, 1] < 0.5)

            p = clf.predict_log_proba([3, 2])
            assert_true(p[0, 1] > p[0, 0])
            p = clf.predict_log_proba([-1, -1])
            assert_true(p[0, 1] < p[0, 0])

    def test_sgd_l1(self):
        """Test L1 regularization"""
        n = len(X4)
        rng = np.random.RandomState(13)
        idx = np.arange(n)
        rng.shuffle(idx)

        X = X4[idx, :]
        Y = Y4[idx, :]

        clf = self.factory(penalty='l1', alpha=.2, fit_intercept=False,
                           n_iter=2000)
        clf.fit(X, Y)
        assert_array_equal(clf.coef_[0, 1:-1], np.zeros((4,)))
        pred = clf.predict(X)
        assert_array_equal(pred, Y)

    def test_class_weights(self):
        """
        Test class weights.
        """
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
        """Test if equal class weights approx. equals no class weights. """
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
        """ValueError due to not existing class label."""
        clf = self.factory(alpha=0.1, n_iter=1000, class_weight={0: 0.5})
        clf.fit(X, Y)

    @raises(ValueError)
    def test_wrong_class_weight_format(self):
        """ValueError due to wrong class_weight argument type."""
        clf = self.factory(alpha=0.1, n_iter=1000, class_weight=[0.5])
        clf.fit(X, Y)

    def test_auto_weight(self):
        """Test class weights for imbalanced data"""
        # compute reference metrics on iris dataset that is quite balanced by
        # default
        X, y = iris.data, iris.target
        X = preprocessing.scale(X)
        idx = np.arange(X.shape[0])
        rng = np.random.RandomState(0)
        rng.shuffle(idx)
        X = X[idx]
        y = y[idx]
        clf = self.factory(alpha=0.0001, n_iter=1000,
                           class_weight=None).fit(X, y)
        assert_approx_equal(metrics.f1_score(y, clf.predict(X)), 0.96, 2)

        # make the same prediction using automated class_weight
        clf_auto = self.factory(alpha=0.0001, n_iter=1000,
                                class_weight="auto").fit(X, y)
        assert_approx_equal(metrics.f1_score(y, clf_auto.predict(X)), 0.96, 2)

        # Make sure that in the balanced case it does not change anything
        # to use "auto"
        assert_array_almost_equal(clf.coef_, clf_auto.coef_, 6)

        # build an very very imbalanced dataset out of iris data
        X_0 = X[y == 0, :]
        y_0 = y[y == 0]

        X_imbalanced = np.vstack([X] + [X_0] * 10)
        y_imbalanced = np.concatenate([y] + [y_0] * 10)

        # fit a model on the imbalanced data without class weight info
        clf = self.factory(n_iter=1000, class_weight=None)
        clf.fit(X_imbalanced, y_imbalanced)
        y_pred = clf.predict(X)
        assert_less(metrics.f1_score(y, y_pred), 0.96)

        # fit a model with auto class_weight enabled
        clf = self.factory(n_iter=1000, class_weight="auto")
        clf.fit(X_imbalanced, y_imbalanced)
        y_pred = clf.predict(X)
        assert_greater(metrics.f1_score(y, y_pred), 0.96)

        # fit another using a fit parameter override
        clf = self.factory(n_iter=1000, class_weight="auto")
        clf.fit(X_imbalanced, y_imbalanced)
        y_pred = clf.predict(X)
        assert_greater(metrics.f1_score(y, y_pred), 0.96)

    def test_sample_weights(self):
        """Test weights on individual samples"""
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
        """Test if ValueError is raised if sample_weight has wrong shape"""
        clf = self.factory(alpha=0.1, n_iter=1000, fit_intercept=False)
        # provided sample_weight too long
        clf.fit(X, Y, sample_weight=range(7))

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
        assert_equal(clf.decision_function([0, 0]).shape, (1, ))
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
        assert_equal(clf.decision_function([0, 0]).shape, (1, 3))
        id1 = id(clf.coef_.data)

        clf.partial_fit(X2[third:], Y2[third:])
        id2 = id(clf.coef_.data)
        # check that coef_ haven't been re-allocated
        assert_true(id1, id2)

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

        clf = self.factory(alpha=0.01, loss="huber")
        clf.fit(X, Y)
        assert_equal(1.0, np.mean(clf.predict(X) == Y))

        clf = self.factory(alpha=0.01, learning_rate="constant", eta0=0.01,
                           loss="squared_loss")
        clf.fit(X, Y)
        assert_equal(1.0, np.mean(clf.predict(X) == Y))


class SparseSGDClassifierTestCase(DenseSGDClassifierTestCase):
    """Run exactly the same tests using the sparse representation variant"""

    factory = SparseSGDClassifier


###############################################################################
# Regression Test Case

class DenseSGDRegressorTestCase(unittest.TestCase):
    """Test suite for the dense representation variant of SGD"""

    factory = SGDRegressor

    def test_sgd(self):
        """Check that SGD gives any results."""
        clf = self.factory(alpha=0.1, n_iter=2,
                           fit_intercept=False)
        clf.fit([[0, 0], [1, 1], [2, 2]], [0, 1, 2])
        assert_equal(clf.coef_[0], clf.coef_[1])

    @raises(ValueError)
    def test_sgd_bad_penalty(self):
        """Check whether expected ValueError on bad penalty"""
        self.factory(penalty='foobar', rho=0.85)

    @raises(ValueError)
    def test_sgd_bad_loss(self):
        """Check whether expected ValueError on bad loss"""
        self.factory(loss="foobar")

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
        y = 0.5 * X.ravel() \
            + np.random.randn(n_samples, 1).ravel()

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
        """Check that the SGD ouput is consistent with coordinate descent"""

        n_samples, n_features = 1000, 5
        rng = np.random.RandomState(0)
        X = np.random.randn(n_samples, n_features)
        # ground_truth linear model that generate y from X and to which the
        # models should converge if the regularizer would be set to 0.0
        ground_truth_coef = rng.randn(n_features)
        y = np.dot(X, ground_truth_coef)

        # XXX: alpha = 0.1 seems to cause convergence problems
        for alpha in [0.01, 0.001]:
            for rho in [0.5, 0.8, 1.0]:
                cd = linear_model.ElasticNet(alpha=alpha, rho=rho,
                                             fit_intercept=False)
                cd.fit(X, y)
                sgd = self.factory(penalty='elasticnet', n_iter=50,
                                   alpha=alpha, rho=rho, fit_intercept=False)
                sgd.fit(X, y)
                err_msg = ("cd and sgd did not converge to comparable "
                           "results for alpha=%f and rho=%f" % (alpha, rho))
                assert_almost_equal(cd.coef_, sgd.coef_, decimal=2,
                                    err_msg=err_msg)

    def test_partial_fit(self):
        third = X.shape[0] // 3
        clf = self.factory(alpha=0.01)

        clf.partial_fit(X[:third], Y[:third])
        assert_equal(clf.coef_.shape, (X.shape[1], ))
        assert_equal(clf.intercept_.shape, (1,))
        assert_equal(clf.decision_function([0, 0]).shape, (1, ))
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
    """Run exactly the same tests using the sparse representation variant"""

    factory = SparseSGDRegressor
