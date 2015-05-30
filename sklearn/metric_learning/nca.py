from ..base import BaseEstimator, TransformerMixin

import numpy as np
import scipy as sp
import scipy.misc
import scipy.optimize
import time

import sys

def nca_vectorized_oracle(X, y, params):
    dX = X[:, None, :] - X[None, :, :]  # n_samples x n_samples x n_comp
    outer = dX[:, :, None, :] * dX[:, :, :, None]  # n_samples x n_samples x n_comp x n_comp

    n_samples, n_features = X.shape
    n_components = params.n_components

    def oracle(L):
        L = L.reshape((n_components, n_features))

        Lx = np.dot(X, L.T)  # n_samples x n_comp
        assert Lx.shape == (n_samples, n_components)

        A = Lx[None, :, :] - Lx[:, None, :]  # n_samples x n_samples x n_comp
        assert A.shape == (n_samples, n_samples, n_components)

        logp = -(A*A).sum(axis=2)  # n_samples x n_samples
        np.fill_diagonal(logp, -np.inf)
        logp -= sp.misc.logsumexp(logp, axis=1)[:, None]
        assert logp.shape == (n_samples, n_samples)

        p = np.exp(logp)  # n_samples x n_samples
        assert p.shape == (n_samples, n_samples)

        class_neighbours = y[:, None] == y[None, :]

        p_i = (p * class_neighbours).sum(axis=1)

        if params.loss == 'l1':
            grad = np.sum(np.sum(p[:, :, None, None] * outer, axis=0) * p_i[:, None, None], axis=0)
        else:
            grad = np.sum(np.sum(p[:, :, None, None] * outer, axis=0), axis=0)

        assert grad.shape == (n_features, n_features)

        if params.loss == 'l1':
            snd = np.sum(class_neighbours[:, :, None, None] * p[:, :, None, None] * outer, axis=(0, 1))
        else:
            snd = np.sum(class_neighbours[:, :, None, None] * p[:, :, None, None] * outer, axis=0) / p_i[:, None, None]
            snd = np.sum(snd, axis=0)

        assert snd.shape == (n_features, n_features)
        grad -= snd

        fnc = p_i.sum() if params.loss == 'l1' else np.log(p_i).sum()

        grad = 2 * np.dot(L, grad)

        return [-fnc, -grad.flatten()]

    return oracle


def nca_semivectorized_oracle(X, y, params, threshold=None):
    n_samples, n_features = X.shape
    n_components = params.n_components

    if threshold is None:
        threshold = 1e-3 / np.abs(X).max() ** 2

    def oracle(L):
        L = L.reshape((n_components, n_features))

        Lx = np.dot(X, L.T)  # n_samples x n_comp
        assert Lx.shape == (n_samples, n_components)

        grad = np.zeros((n_features, n_features))
        fnc = 0

        for i in range(n_samples):
            A = Lx[i, :] - Lx  # n_samples x n_comp
            assert A.shape == (n_samples, n_components)

            logp = -(A*A).sum(axis=1)  # n_samples
            logp[i] = -np.inf
            logp -= sp.misc.logsumexp(logp)
            assert logp.shape == (n_samples, )

            p = np.exp(logp)  # n_samples
            assert p.shape == (n_samples, )

            class_neighbours = y == y[i]
            p_i = (p * class_neighbours).sum()

            a = p * (p_i if params.loss == 'l1' else 1)  # n_samples
            a_mask = a > threshold

            Xij = X[i] - X[a_mask, :]  # n_relevant x n_features
            grad += (a[a_mask, None, None] * (Xij[:, :, None] * Xij[:, None, :])).sum(axis=0)

            b = p / (1 if params.loss == 'l1' else p_i) * (y == y[i])  # n_samples
            b_mask = b > threshold

            Xij = X[i] - X[b_mask, :]  # n_relevant x n_features
            grad -= (b[b_mask, None, None] * (Xij[:, :, None] * Xij[:, None, :])).sum(axis=0)

            fnc += p_i if params.loss == 'l1' else np.log(p_i)

        grad = 2 * np.dot(L, grad)

        return [-fnc, -grad.flatten()]

    return oracle

def optimize_nca(X, y, params):
    n_samples, n_features = X.shape

    rng = np.random.RandomState(params.random_state)
    best_fnc = np.inf
    best_L = None

    outer_products_size = (n_samples * n_features) ** 2
    memory_permits = outer_products_size < params.cache_size
    if memory_permits:
        oracle = nca_vectorized_oracle(X, y, params)
    else:
        oracle = nca_semivectorized_oracle(X, y, params)

    for _ in range(params.n_init):
        start = time.clock()
        L = rng.uniform(0, 1, (params.n_components, n_features))

        if params.method in ['gd', 'adagrad']:
            is_adagrad = params.method == 'adagrad'
            grad_history = 1e-3 * np.ones_like(L) if is_adagrad else None

            for it in range(params.max_iter):
                fnc, grad = oracle(L)
                grad = grad.reshape(L.shape)
                if is_adagrad:
                    grad_history += grad * grad
                    grad /= np.sqrt(grad_history)

                L -= params.learning_rate * grad
                if params.verbose > 1:
                    percent = 100. * (it + 1) / params.max_iter
                    sys.stdout.write("\rProgress {:.2f}% :: Target = {}".format(percent, fnc))
                    sys.stdout.flush()

                # FIXME?
                if np.abs(grad).max() < params.tol:
                    break

        elif params.method == 'scipy':
            state = {'it' : 0}
            def callback(L):
                fnc, _ = oracle(L)
                state['it'] += 1
                if params.verbose > 1:
                    percent = 100. * (state['it'] + 1) / params.max_iter
                    sys.stdout.write("\rProgress {:.2f}% :: Target = {}".format(percent, fnc))
                    sys.stdout.flush()

            options = {'maxiter' : params.max_iter, "disp" : params.verbose > 1}
            res = sp.optimize.minimize(fun=oracle, x0=L, jac=True, tol=params.tol,
                                       options=options, callback=callback)
            L = res.x.reshape(L.shape)
            fnc, grad = oracle(L)

        exec_time = time.clock() - start
        if params.verbose > 0:
            sys.stdout.write("\rInit {} / {} completed"
                             " :: Time: {:.2f}s :: Loss = {}\n".format(_+1, params.n_init, exec_time, fnc))
            sys.stdout.flush()

        if fnc < best_fnc:
            best_fnc = fnc
            best_L = L

    if params.verbose > 0:
        print("\rCompleted" + " " * 80)

    return best_L

class BaseNCA(BaseEstimator):
    def __init__(self, n_components=None, method='adagrad', learning_rate=1.0, tol=1e-5, loss='kl',
                 max_iter=100, n_init=1, random_state=None, verbose=0, cache_size=100*1024*1024):
        self.n_components = n_components
        self.random_state = random_state
        self.max_iter = max_iter
        self.learning_rate = learning_rate
        self.n_init = n_init
        self.loss = loss
        self.verbose = verbose
        self.method = method
        self.tol = tol
        self.cache_size = cache_size

    def fit(self, X, y):
        n_features = X.shape[1]
        params = self.get_params()

        if params['method'] not in ['adagrad', 'gd', 'scipy']:
            raise ValueError("Unsupported optimization method: %s" % params['method'])

        if self.n_components is None:
            params['n_components'] = n_features

        from collections import namedtuple
        params = namedtuple('ParamSet', params.keys())(*params.values())

        self.L = optimize_nca(X, y, params)
        return self

    def transform(self, X):
        return np.dot(X, self.L.T)