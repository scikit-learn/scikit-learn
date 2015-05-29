from ..base import BaseEstimator, TransformerMixin

import numpy as np
import scipy as sp
import scipy.misc
import scipy.optimize

def vectorized_nca(X, y, params):
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

    n_samples, n_features = X.shape

    rng = np.random.RandomState(params.random_state)
    best_fnc = None
    best_L = None

    for _ in range(params.n_init):
        L = rng.uniform(0, 1, (params.n_components, n_features))

        if params.method == 'gd':
            for it in range(params.max_iter):
                fnc, grad = oracle(L)
                L -= params.learning_rate * grad.reshape(L.shape)
                if params.verbose:
                    print("Iteration {} :: Target = {}".format(it + 1, fnc))

        elif params.method == 'scipy':
            state = {'it' : 0}
            def callback(L):
                fnc, _ = oracle(L)
                state['it'] += 1
                if params.verbose:
                    print("Iteration {} :: Target = {}".format(state['it'], fnc))

            options = {'maxiter' : params.max_iter, "disp" : True}
            res = sp.optimize.minimize(fun=oracle, method='CG', x0=L, jac=True, tol=params.tol,
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
        self.L = vectorized_nca(X, y, self)
        return self

    def transform(self, X):
        return np.dot(X, self.L.T)