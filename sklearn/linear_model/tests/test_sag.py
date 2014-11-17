import pickle
import unittest

import numpy as np
import scipy.sparse as sp
from sklearn.linear_model import SAGRegressor, SAGClassifier
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal


class SparseSAGClassifier(SAGClassifier):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGClassifier.fit(X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGClassifier.partial_fit(X, y, *args, **kw)

    def decision_function(self, X):
        X = sp.csr_matrix(X)
        return SAGClassifier.decision_function(X)

    def predict_proba(self, X):
        X = sp.csr_matrix(X)
        return SAGClassifier.predict_proba(X)


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

    def sag(self, X, y, eta, alpha, intercept_init=0.0, indexes=None):
        n_samples, n_features = X.shape[0], X.shape[1]

        weights = np.zeros(X.shape[1])
        sum_gradient = np.zeros(X.shape[1])
        gradient_memory = np.zeros((n_samples, n_features))
        intercept = intercept_init
        rng = np.random.RandomState(77)
        decay = 1.0
        seen = set()
        if indexes is None:
            indexes = range(n_samples)

        # sparse data has a fixed decay of .01
        # if (isinstance(self, SparseSGDClassifierTestCase) or
        #         isinstance(self, SparseSGDRegressorTestCase)):
        #     decay = .01

        for i in indexes:
            # idx = int(rng.rand(1) * n_samples)
            idx = i
            entry = X[idx]
            seen.add(idx)
            p = np.dot(entry, weights) + intercept
            gradient = p - y[idx]
            update = entry * gradient + alpha * weights
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            intercept -= eta * gradient * decay
            weights -= eta * sum_gradient / len(seen)

        return weights, intercept

    def sag_sparse(self, X, y, eta, alpha, intercept_init=0.0, indexes=None):
        n_samples, n_features = X.shape[0], X.shape[1]

        weights = np.zeros(n_features)
        sum_gradient = np.zeros(n_features)
        last_updated = np.zeros(n_features, dtype=np.int)
        c_sum = np.zeros(n_samples + 1)
        gradient_memory = np.zeros((n_samples, n_features))
        intercept = intercept_init
        wscale = 1.0
        decay = 1.0
        seen = set()
        if indexes is None:
            indexes = range(n_samples)

        # sparse data has a fixed decay of .01
        # if (isinstance(self, SparseSGDClassifierTestCase) or
        #         isinstance(self, SparseSGDRegressorTestCase)):
        #     decay = .01

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
            gradient = p - y[idx]

            update = entry * gradient
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            wscale *= (1.0 - alpha * eta)
            c_sum[counter + 1] = c_sum[counter] + eta / (wscale * len(seen))

            intercept -= eta * gradient * decay
            counter += 1

        for k in range(n_features):
            weights[k] -= (c_sum[counter] - c_sum[last_updated[k]]) * sum_gradient[k]
        weights *= wscale

        return weights, intercept


class DenseSAGRegressorTestCase(unittest.TestCase, CommonTest):
    factory = SAGRegressor

    def test_sag_computed_correctly(self):
        eta = .2
        alpha = .1
        n_features = 20
        n_samples = 10
        clf1 = self.factory(learning_rate="constant", eta0=eta,
                            alpha=alpha, n_iter=1, random_state=77)
        rng = np.random.RandomState(0)
        X = rng.normal(size=(n_samples, n_features))
        w = rng.normal(size=n_features)
        y = np.dot(X, w)

        clf1.fit(X, y)
        # weights, intercept = self.sag(X, y, eta, alpha,
        #                               indexes=[0, 1, 1, 3, 0, 8, 6, 2, 3, 9])

        indexes = [0, 1, 1, 3, 0, 8, 6, 2, 3, 9]

        spweights, spintercept = self.sag_sparse(X, y, eta, alpha,
                                                 indexes=indexes)

        assert_array_almost_equal(clf1.coef_.ravel(),
                                  spweights.ravel(),
                                  decimal=14)
        assert_almost_equal(clf1.intercept_, spintercept, decimal=14)

    # def test_sag_binary_computed_correctly(self):
    #     eta = .1
    #     alpha = 0.1
    #     n_samples = 20
    #     n_features = 10
    #     rng = np.random.RandomState(0)
    #     X = rng.normal(size=(n_samples, n_features))
    #     w = rng.normal(size=n_features)

    #     clf = self.factory(loss='squared_loss',
    #                        learning_rate='constant',
    #                        eta0=eta, alpha=alpha,
    #                        fit_intercept=True,
    #                        n_iter=1, sag=True)

    #     # simple linear function without noise
    #     y = np.dot(X, w)
    #     y = np.sign(y)

    #     clf.fit(X, y)

    #     sag_weights, sag_intercept = self.sag(X, y, eta, alpha)
    #     sag_weights = sag_weights.reshape(1, -1)
    #     assert_array_almost_equal(clf.coef_,
    #                               sag_weights,
    #                               decimal=14)
    #     assert_almost_equal(clf.intercept_, sag_intercept, decimal=14)

class SparseSAGRegressorTestCase(DenseSAGRegressorTestCase):
    """Run exactly the same tests using the sparse representation variant"""

    factory = SparseSAGRegressor
