from ..base import BaseEstimator, TransformerMixin
from ..utils.validation import check_is_fitted

import numpy as np
import scipy as sp
import scipy.misc
import scipy.optimize
import time

import sys



def nca_vectorized_oracle(X, y, n_components, loss, threshold=0):
    dx = X[:, np.newaxis, :] - X[np.newaxis, :, :]  # n_samples x n_samples x n_features
    outer = dx[:, :, np.newaxis, :] * dx[:, :, :, np.newaxis]  # n_samples x n_samples x n_features x n_features

    n_samples, n_features = X.shape

    def oracle(L):
        L = L.reshape((n_components, n_features))

        Lx = np.dot(X, L.T)  # n_samples x n_components
        assert Lx.shape == (n_samples, n_components)

        A = Lx[np.newaxis, :, :] - Lx[:, np.newaxis, :]  # n_samples x n_samples x n_components
        assert A.shape == (n_samples, n_samples, n_components)

        logp = -(A * A).sum(axis=2)  # n_samples x n_samples
        np.fill_diagonal(logp, -np.inf)
        logp -= sp.misc.logsumexp(logp, axis=1)[:, np.newaxis]
        assert logp.shape == (n_samples, n_samples)

        p = np.exp(logp)  # n_samples x n_samples
        assert p.shape == (n_samples, n_samples)

        class_neighbours = y[:, np.newaxis] == y[np.newaxis, :]

        p_i = (p * class_neighbours).sum(axis=1)


        proba = np.copy(p)
        if loss == 'l1':
            proba *= p_i[:, np.newaxis]
        mask = proba > threshold
        proba *= mask

        grad = np.sum(proba[:, :, np.newaxis, np.newaxis] * outer, axis=1).sum(axis=0)
        assert grad.shape == (n_features, n_features)

        neighbours_proba = class_neighbours * p
        if loss == 'kl':
            neighbours_proba /= p_i[:, np.newaxis]

        neighbours_mask = neighbours_proba > threshold
        neighbours_proba *= neighbours_mask

        neighbours_term = np.sum(neighbours_proba[:, :, np.newaxis, np.newaxis] * outer, axis=1).sum(axis=0)

        assert neighbours_term.shape == (n_features, n_features)
        grad -= neighbours_term

        if loss == 'l1':
            function_value = p_i.sum()
        elif loss == 'kl':
            function_value = np.log(p_i).sum()

        grad = 2 * np.dot(L, grad)

        return [-function_value, -grad.flatten()]

    return oracle


def nca_semivectorized_oracle(X, y, n_components, loss, threshold=0):
    n_samples, n_features = X.shape

    def oracle(L):
        L = L.reshape((n_components, n_features))

        Lx = np.dot(X, L.T)  # n_samples x n_components
        assert Lx.shape == (n_samples, n_components)

        grad = np.zeros((n_features, n_features))
        function_value = 0

        for i in range(n_samples):
            A = Lx[i, :] - Lx  # n_samples x n_components
            assert A.shape == (n_samples, n_components)

            logp = -(A * A).sum(axis=1)  # n_samples
            logp[i] = -np.inf
            logp -= sp.misc.logsumexp(logp)
            assert logp.shape == (n_samples, )

            p = np.exp(logp)  # n_samples
            assert p.shape == (n_samples, )

            class_neighbours = y == y[i]
            p_i = (p * class_neighbours).sum()

            if loss == 'l1':
                samples_proba = p * p_i
            elif loss == 'kl':
                samples_proba = p
            mask = samples_proba > threshold  # n_samples

            Xij = X[i] - X[mask, :]  # n_relevant x n_features
            Xij_outer = (Xij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
            grad += (samples_proba[mask, np.newaxis, np.newaxis] * Xij_outer).sum(axis=0)

            if loss == 'l1':
                neighbours_proba = p * (y == y[i])
            elif loss == 'kl':
                neighbours_proba = p / p_i * (y == y[i])
            mask = neighbours_proba > threshold  # n_samples

            Xij = X[i] - X[mask, :]  # n_relevant x n_features
            Xij_outer = (Xij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
            grad -= (neighbours_proba[mask, np.newaxis, np.newaxis] * Xij_outer).sum(axis=0)

            if loss == 'l1':
                function_value += p_i
            elif loss == 'kl':
                function_value += np.log(p_i)

        grad = 2 * np.dot(L, grad)

        return [-function_value, -grad.flatten()]

    return oracle


def optimize_nca(X, y, learning_rate, n_components, loss, n_init, max_iter,
                 method, tol, random_state, cache_size, verbose):
    n_samples, n_features = X.shape

    rng = np.random.RandomState(random_state)
    best_value = np.inf
    best_L = None


    threshold = 1e-3 / np.abs(X).max() ** 2

    outer_products_size = (n_samples * n_features) ** 2
    memory_permits = outer_products_size < cache_size
    if memory_permits:
        oracle = nca_vectorized_oracle(X, y, n_components, loss, threshold=threshold)
    else:
        oracle = nca_semivectorized_oracle(X, y, n_components, loss, threshold=threshold)

    for _ in range(n_init):
        start = time.clock()
        L = rng.uniform(0, 1, (n_components, n_features))

        if method in ['gd', 'adagrad']:
            is_adagrad = method == 'adagrad'
            grad_history = 1e-3 * np.ones_like(L) if is_adagrad else None

            for it in range(max_iter):
                value, grad = oracle(L)
                grad = grad.reshape(L.shape)
                if is_adagrad:
                    grad_history += grad * grad
                    grad /= np.sqrt(grad_history)

                L -= learning_rate * grad
                if verbose > 1:
                    percent = 100. * (it + 1) / max_iter
                    sys.stdout.write("\rProgress {:.2f}% :: Target = {}".format(percent, value))
                    sys.stdout.flush()

                # FIXME?
                if np.abs(grad).max() < tol:
                    break

        elif method == 'scipy':
            state = {'it': 0}

            def callback(L):
                state['it'] += 1
                if verbose > 1:
                    value, _ = oracle(L)
                    percent = 100. * (state['it'] + 1) / max_iter
                    sys.stdout.write("\rProgress {:.2f}% :: Target = {}".format(percent, value))
                    sys.stdout.flush()

            options = {'maxiter': max_iter, "disp": verbose > 1}
            res = sp.optimize.minimize(fun=oracle, x0=L, jac=True, tol=tol,
                                       options=options, callback=callback)
            L = res.x.reshape(L.shape)
            value, grad = oracle(L)

        exec_time = time.clock() - start
        if verbose > 0:
            sys.stdout.write("\rInit {} / {} completed"
                             " :: Time: {:.2f}s :: Loss = {}\n".format(_ + 1, n_init, exec_time, value))
            sys.stdout.flush()

        if value < best_value:
            best_value = value
            best_L = L

    if verbose > 0:
        print("\rCompleted" + " " * 80)

    return best_L


class BaseNCA(BaseEstimator):
    def __init__(self, n_components, method, learning_rate, tol, loss, max_iter,
                 n_init, random_state, verbose, cache_size):
        self.n_components = n_components
        self.random_state = random_state
        self.max_iter = max_iter
        self.learning_rate = learning_rate  # TODO: learning rate scheduling from sgd_fast?
        self.n_init = n_init
        self.loss = loss
        self.verbose = verbose
        self.method = method
        self.tol = tol
        self.cache_size = cache_size

    def fit(self, X, y):
        n_features = X.shape[1]

        if self.method not in ['adagrad', 'gd', 'scipy']:
            raise ValueError("Unsupported optimization method: %s" % self.method)

        # TODO: more validation

        n_components = self.n_components or n_features

        self.matrix_ = optimize_nca(X, y, self.learning_rate, n_components, self.loss,
                                    self.n_init, self.max_iter, self.method, self.tol,
                                    self.random_state, self.cache_size, self.verbose)
        return self


class NCATransformer(BaseNCA, TransformerMixin):
    # TODO
    """Neighbourhood Components Analysis-based transformer

    Transforms the data space using the linear transformation learned
    by NCA

    Parameters
    ----------
    n_components: int, or None (default)
        Specifies dimensionality of the target space. If None,
        n_features will be used.

    method: string

    learning_rate: float

    tol: float

    loss: string

    max_iter: int

    n_init: int

    random_state: int seed, RandomState instance, or None (default)

    cache_size:

    verbose: int

    """
    def __init__(self, n_components=None, method='adagrad', learning_rate=1.0, tol=1e-5, loss='kl',
                 max_iter=100, n_init=1, random_state=None, verbose=0, cache_size=100 * 1024 * 1024):
        super(NCATransformer, self).__init__(n_components, method, learning_rate, tol, loss,
                                             max_iter, n_init, random_state, verbose, cache_size)

    def transform(self, X):
        check_is_fitted(self, 'matrix_')
        return np.dot(X, self.matrix_.T)

# TODO: add NCASimilarity