from ..base import BaseEstimator, TransformerMixin

import numpy as np
import scipy as sp
import scipy.misc
import scipy.optimize

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

        grad = np.sum(np.sum(p[:, :, None, None] * outer, axis=0) * p_i[:, None, None], axis=0)
        assert grad.shape == (n_features, n_features)

        snd = np.sum(class_neighbours[:, :, None, None] * p[:, :, None, None] * outer, axis=(0, 1))
        assert snd.shape == (n_features, n_features)
        grad -= snd

        fnc = p_i.sum() if params.loss == 'l1' else np.log(p_i).sum()

        grad = 2 * np.dot(L, grad)

        return [-fnc, -grad.flatten()]

    return oracle


def nca_loopy_oracle(X, y, params, threshold=None):
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

            logp = -(A*A).sum(axis=2)  # n_samples
            logp[i] = -np.inf
            logp -= sp.misc.logsumexp(logp)
            assert logp.shape == (n_samples, )

            p = np.exp(logp)  # n_samples
            assert p.shape == (n_samples, )

            class_neighbours = y == y[i]
            p_i = (p * class_neighbours).sum()

            for j in range(n_samples):
                a = p[j] * (p_i if params.loss == 'l1' else 1)
                if a > threshold:
                    Xij = X[i] - X[j]
                    grad += a * np.outer(Xij, Xij)

                b = p[j] / (1 if params.loss == 'l1' else p_i)
                if b > threshold and y[j] == y[i]:
                    Xij = X[i] - X[j]
                    grad -= b * np.outer(Xij, Xij)

            fnc += p_i if params.loss == 'l1' else np.log(p_i)

        grad = 2 * np.dot(L, grad)

        return [-fnc, -grad.flatten()]

    return oracle

def optimize_nca(X, y, params, cache_size=200*1024*1024):
    n_samples, n_features = X.shape

    rng = np.random.RandomState(params.random_state)
    best_fnc = None
    best_L = None

    outer_products_size = (n_samples * n_features) ** 2
    memory_permits = outer_products_size < cache_size
    if memory_permits:
        oracle = nca_vectorized_oracle(X, y, params)
    else:
        oracle = nca_loopy_oracle(X, y, params)

    for _ in range(params.n_init):
        L = rng.uniform(0, 1, (params.n_components, n_features))

        if params.method == 'gd':
            for it in range(params.max_iter):
                fnc, grad = oracle(L)
                grad = grad.reshape(L.shape)
                L -= params.learning_rate * grad
                if params.verbose:
                    print("Iteration {} :: Target = {}".format(it + 1, fnc))

        elif params.method == 'adagrad':
            grad_history = np.zeros_like(L)
            for it in range(params.max_iter):
                fnc, grad = oracle(L)
                grad = grad.reshape(L.shape)
                grad_history += grad * grad
                L -= params.learning_rate / np.sqrt(grad_history) * grad
                if params.verbose > 0:
                    print("Iteration {} :: Target = {}".format(it + 1, fnc))

        elif params.method == 'scipy':
            state = {'it' : 0}
            def callback(L):
                fnc, _ = oracle(L)
                state['it'] += 1
                if params.verbose:
                    print("Iteration {} :: Target = {}".format(state['it'], fnc))

            options = {'maxiter' : params.max_iter, "disp" : params.verbose > 0}
            res = sp.optimize.minimize(fun=oracle, x0=L, jac=True, tol=params.tol,
                                       options=options, callback=callback)
            L = res.x.reshape(L.shape)
            fnc, grad = oracle(L)

        if fnc > best_fnc:
            best_fnc = fnc
            best_L = L

    return best_L


class BaseNCA(BaseEstimator):
    def __init__(self, n_components, method, tol, loss, learning_rate, max_iter, n_init, random_state, verbose):
        self.n_components = n_components
        self.random_state = random_state
        self.max_iter = max_iter
        self.learning_rate = learning_rate
        self.n_init = n_init
        self.loss = loss
        self.verbose = verbose
        self.method = method
        self.tol = tol

    def fit(self, X, y):
        self.L = optimize_nca(X, y, self)
        return self

    def transform(self, X):
        return np.dot(X, self.L.T)