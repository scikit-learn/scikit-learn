import unittest
import math

import numpy as np
import scipy.sparse as sp
from sklearn.linear_model import SAGRegressor, SAGClassifier
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises_regexp
from sklearn.datasets import make_blobs

# this is used for sag classification
def log_dloss(p, y):
    z = p * y
    # approximately equal and saves the computation of the log
    if z > 18.0:
        return math.exp(-z) * -y
    if z < -18.0:
        return -y
    return -y / (math.exp(z) + 1.0)


# this is used for sag regression
def squared_dloss(p, y):
    return p - y


class SparseSAGClassifier(SAGClassifier):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGClassifier.fit(self, X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGClassifier.partial_fit(self, X, y, *args, **kw)

    def decision_function(self, X):
        X = sp.csr_matrix(X)
        return SAGClassifier.decision_function(self, X)

    def predict_proba(self, X):
        X = sp.csr_matrix(X)
        return SAGClassifier.predict_proba(self, X)


class SparseSAGRegressor(SAGRegressor):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGRegressor.fit(self, X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGRegressor.partial_fit(self, X, y, *args, **kw)

    def decision_function(self, X, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGRegressor.decision_function(self, X, *args, **kw)


class CommonTest(object):

    # def sag(self, X, y, eta, alpha, intercept_init=0.0,
    #         indexes=None, dloss=None):
    #     n_samples, n_features = X.shape[0], X.shape[1]

    #     weights = np.zeros(X.shape[1])
    #     sum_gradient = np.zeros(X.shape[1])
    #     gradient_memory = np.zeros((n_samples, n_features))
    #     intercept = intercept_init
    #     rng = np.random.RandomState(77)
    #     decay = 1.0
    #     seen = set()
    #     if indexes is None:
    #         indexes = range(n_samples)

    #     # sparse data has a fixed decay of .01
    #     # if (isinstance(self, SparseSGDClassifierTestCase) or
    #     #         isinstance(self, SparseSGDRegressorTestCase)):
    #     #     decay = .01

    #     for i in indexes:
    #         # idx = int(rng.rand(1) * n_samples)
    #         idx = i
    #         entry = X[idx]
    #         seen.add(idx)
    #         p = np.dot(entry, weights) + intercept
    #         gradient = dloss(p, y[idx])
    #         update = entry * gradient + alpha * weights
    #         sum_gradient += update - gradient_memory[idx]
    #         gradient_memory[idx] = update

    #         intercept -= eta * gradient * decay
    #         weights -= eta * sum_gradient / len(seen)

    #     return weights, intercept

    def sag_sparse(self, X, y, eta, alpha, intercept_init=0.0, n_iter=1,
                   indexes=None, dloss=None, class_weight=None):
        n_samples, n_features = X.shape[0], X.shape[1]

        weights = np.zeros(n_features)
        sum_gradient = np.zeros(n_features)
        last_updated = np.zeros(n_features, dtype=np.int)
        gradient_memory = np.zeros((n_samples, n_features))
        intercept = intercept_init
        wscale = 1.0
        decay = 1.0
        seen = set()

        if indexes is None:
            indexes = np.tile(np.arange(n_samples), n_iter)
        c_sum = np.zeros(len(indexes) + 1)

        # sparse data has a fixed decay of .01
        if (isinstance(self, SparseSAGClassifierTestCase) or
                isinstance(self, SparseSAGRegressorTestCase)):
            decay = .01

        counter = 0
        for k in indexes:
            # idx = int(rng.rand(1) * n_samples)
            idx = k
            entry = X[idx]
            seen.add(idx)

            for j in range(n_features):
                weights[j] -= ((c_sum[counter] - c_sum[last_updated[j]]) *
                               sum_gradient[j])
                last_updated[j] = counter

            p = (wscale * np.dot(entry, weights)) + intercept
            gradient = dloss(p, y[idx])

            if class_weight:
                gradient *= class_weight[y[idx]]

            update = entry * gradient
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            wscale *= (1.0 - alpha * eta)
            c_sum[counter + 1] = c_sum[counter] + eta / (wscale * len(seen))

            intercept -= eta * gradient * decay
            counter += 1

        for k in range(n_features):
            weights[k] -= (c_sum[counter] -
                           c_sum[last_updated[k]]) * sum_gradient[k]
        weights *= wscale

        return weights, intercept


class DenseSAGRegressorTestCase(unittest.TestCase, CommonTest):
    factory = SAGRegressor

    def test_sag_computed_correctly(self):
        """tests the if the sag regressor computed correctly"""
        eta = .001
        alpha = .1
        n_features = 20
        n_samples = 100
        n_iter = 60
        clf1 = self.factory(eta0=eta, alpha=alpha,
                            n_iter=n_iter, random_state=77)
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)
        y = np.dot(X, w)

        clf1.fit(X, y)

        indexes = np.random.random_integers(0, n_samples - 1,
                                            n_iter * n_samples)

        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes,
                                                 dloss=squared_dloss)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    def test_regressor_warm_start(self):
        """tests the regressor warmstart"""
        eta = .001
        alpha = .1
        n_features = 20
        n_samples = 100
        n_iter = 100
        clf1 = self.factory(eta0=eta, alpha=alpha,
                            n_iter=n_iter, random_state=77)
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)
        y = np.dot(X, w)

        clf1.fit(X, y)
        clf1.fit(X, y)

        indexes = np.random.random_integers(0, n_samples - 1,
                                            n_iter * n_samples)

        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes,
                                                 dloss=squared_dloss)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    def test_auto_eta(self):
        """tests the auto eta computed correctly"""
        indexes = [2, 0, 0, 0, 2, 2]
        X = np.array([[1, 2, 3], [2, 3, 4], [2, 3, 2]])
        y = [.5, .6, .7]
        alpha = 1.0
        n_iter = 600
        indexes = np.random.random_integers(0, 2, n_iter * 3)
        # sum the squares of the second sample because that's the largest
        eta = 4 + 9 + 16
        eta = 1.0 / (eta + alpha)

        clf1 = self.factory(eta0='auto', alpha=alpha,
                            n_iter=n_iter, random_state=77)
        clf1.fit(X, y)
        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes,
                                                 dloss=squared_dloss)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    def test_sag_regressor(self):
        """tests the if the sag regressor performs well"""
        xmin, xmax = -5, 5
        n_samples = 100
        rng = np.random.RandomState(0)
        X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

        # simple linear function without noise
        y = 0.5 * X.ravel()

        clf = self.factory(eta0=.001, alpha=0.1, n_iter=20)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_greater(score, 0.99)

        # simple linear function with noise
        y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

        clf = self.factory(eta0=.001, alpha=0.1, n_iter=20)
        clf.fit(X, y)
        score = clf.score(X, y)
        assert_greater(score, 0.5)


class SparseSAGRegressorTestCase(DenseSAGRegressorTestCase):
    """Run exactly the same tests using the sparse representation variant"""

    factory = SparseSAGRegressor


class DenseSAGClassifierTestCase(unittest.TestCase, CommonTest):
    factory = SAGClassifier

    def test_sag_computed_correctly(self):
        """tests the binary classifier computed correctly"""
        eta = .001
        alpha = .1
        n_samples = 50
        n_iter = 20
        clf1 = self.factory(eta0=eta, alpha=alpha,
                            n_iter=n_iter, random_state=77)
        X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                          cluster_std=0.1)

        classes = np.unique(y)
        y_tmp = np.ones(n_samples)
        y_tmp[y != classes[1]] = -1
        y = y_tmp

        clf1.fit(X, y)

        indexes = np.random.random_integers(0, n_samples - 1,
                                            n_iter * n_samples)

        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes,
                                                 dloss=log_dloss)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    def test_sag_multiclass_computed_correctly(self):
        """tests the multiclass classifier is computed correctly"""
        eta = .001
        alpha = .1
        n_samples = 50
        n_iter = 40
        clf1 = self.factory(eta0=eta, alpha=alpha,
                            n_iter=n_iter, random_state=77)
        rng = np.random.RandomState(0)
        X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                          cluster_std=0.1)
        classes = np.unique(y)

        clf1.fit(X, y)

        indexes = rng.random_integers(0, n_samples - 1,
                                      n_iter * n_samples)

        coef = []
        intercept = []
        for cl in classes:
            y_encoded = np.ones(n_samples)
            y_encoded[y != cl] = -1

            spweights, spintercept = self.sag_sparse(X, y_encoded, eta, alpha,
                                                     indexes=indexes,
                                                     dloss=log_dloss)
            coef.append(spweights)
            intercept.append(spintercept)

        coef = np.vstack(coef)
        intercept = np.array(intercept)

        for i, cl in enumerate(classes):
            assert_array_almost_equal(clf1.coef_[i].ravel(),
                                      coef[i].ravel(),
                                      decimal=2)
            assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)

    def test_classifier_results(self):
        """tests classifier results match target"""
        eta = .2
        alpha = .1
        n_features = 20
        n_samples = 10
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)
        y = np.dot(X, w)
        y = np.sign(y)
        clf = self.factory(eta0=eta, alpha=alpha, n_iter=20)
        clf.fit(X, y)
        pred = clf.predict(X)
        assert_almost_equal(pred, y, decimal=12)

    def test_binary_classifier_warm_start(self):
        """tests binary classifier with a warm start"""
        eta = .001
        alpha = .1
        n_samples = 50
        n_iter = 30
        clf1 = self.factory(eta0=eta, alpha=alpha,
                            n_iter=n_iter, random_state=77, warm_start=True)
        rng = np.random.RandomState(0)
        X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                          cluster_std=0.1)

        classes = np.unique(y)
        y_tmp = np.ones(n_samples)
        y_tmp[y != classes[1]] = -1
        y = y_tmp

        clf1.fit(X, y)
        clf1.fit(X, y)

        indexes = rng.random_integers(0, n_samples - 1,
                                      2 * n_iter * n_samples)

        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes,
                                                 dloss=log_dloss)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    def test_multiclass_classifier_warm_start(self):
        """tests multiclass classifier with a warm start"""
        eta = .001
        alpha = .1
        n_samples = 20
        n_iter = 150
        clf1 = self.factory(eta0=eta, alpha=alpha,
                            n_iter=n_iter, random_state=77, warm_start=True)
        rng = np.random.RandomState(0)
        X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                          cluster_std=0.1)
        classes = np.unique(y)

        clf1.fit(X, y)
        clf1.fit(X, y)

        indexes = rng.random_integers(0, n_samples - 1,
                                      2 * n_iter * n_samples)

        coef = []
        intercept = []
        for cl in classes:
            y_encoded = np.ones(n_samples)
            y_encoded[y != cl] = -1

            spweights, spintercept = self.sag_sparse(X, y_encoded, eta, alpha,
                                                     indexes=indexes,
                                                     dloss=log_dloss)
            coef.append(spweights)
            intercept.append(spintercept)

        coef = np.vstack(coef)
        intercept = np.array(intercept)

        for i, cl in enumerate(classes):
            assert_array_almost_equal(clf1.coef_[i].ravel(),
                                      coef[i].ravel(),
                                      decimal=2)
            assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)

    def test_binary_classifier_class_weight(self):
        """tests binary classifier with classweights for each class"""
        eta = .001
        alpha = .1
        n_samples = 50
        n_iter = 20
        rng = np.random.RandomState(0)
        X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                          cluster_std=0.1)
        classes = np.unique(y)
        y_tmp = np.ones(n_samples)
        y_tmp[y != classes[1]] = -1
        y = y_tmp

        class_weight = {1: .45, -1: .55}

        clf1 = self.factory(eta0=eta, alpha=alpha, n_iter=n_iter,
                            random_state=77, warm_start=True,
                            class_weight=class_weight)

        clf1.fit(X, y)

        indexes = rng.random_integers(0, n_samples - 1,
                                      n_iter * n_samples)

        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes,
                                                 dloss=log_dloss,
                                                 class_weight=class_weight)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    def test_multiclass_classifier_class_weight(self):
        """tests multiclass with classweights for each class"""

        eta = .001
        alpha = .1
        n_samples = 20
        n_iter = 100
        class_weight = {0: .45, 1: .55, 2: .75}

        rng = np.random.RandomState(0)
        X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                          cluster_std=0.1)
        classes = np.unique(y)

        clf1 = self.factory(eta0=eta, alpha=alpha, n_iter=n_iter,
                            random_state=77, warm_start=True,
                            class_weight=class_weight)

        clf1.fit(X, y)

        indexes = rng.random_integers(0, n_samples - 1, n_iter * n_samples)

        coef = []
        intercept = []
        for cl in classes:
            y_encoded = np.ones(n_samples)
            y_encoded[y != cl] = -1
            cl_weight = {-1: 1.0, 1: 1.0}
            cl_weight[1] = class_weight[cl]

            spweights, spintercept = self.sag_sparse(X, y_encoded, eta, alpha,
                                                     indexes=indexes,
                                                     dloss=log_dloss,
                                                     class_weight=cl_weight)
            coef.append(spweights)
            intercept.append(spintercept)

        coef = np.vstack(coef)
        intercept = np.array(intercept)

        for i, cl in enumerate(classes):
            assert_array_almost_equal(clf1.coef_[i].ravel(),
                                      coef[i].ravel(),
                                      decimal=2)
            assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)

    def test_classifier_single_class(self):
        """tests value error thrown with only one class"""
        X = [[1, 2], [3, 4]]
        Y = [1, 1]

        assert_raises_regexp(ValueError,
                             "The number of class labels must be "
                             "greater than one.",
                             self.factory().fit,
                             X, Y)


class SparseSAGClassifierTestCase(DenseSAGClassifierTestCase):
    """Run exactly the same tests using the sparse representation variant"""

    factory = SparseSAGClassifier
