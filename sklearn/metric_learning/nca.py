from ..base import BaseEstimator, TransformerMixin
from ..metrics import euclidean_distances
from ..utils.validation import check_is_fitted

import numpy as np
import scipy as sp
import scipy.misc
import scipy.optimize
import time

import sys


def nca_vectorized_oracle(L, X, y, n_components, loss, outer, threshold=0.0):
    n_samples, n_features = X.shape

    L = L.reshape((n_components, n_features))

    Lx = np.dot(X, L.T)  # n_samples x n_components
    assert Lx.shape == (n_samples, n_components)

    p = -euclidean_distances(Lx, squared=True)
    np.fill_diagonal(p, -np.inf)
    p -= sp.misc.logsumexp(p, axis=1)[:, np.newaxis]
    assert p.shape == (n_samples, n_samples)

    np.exp(p, out=p)  # n_samples x n_samples

    class_neighbours = y[:, np.newaxis] == y[np.newaxis, :]

    p_i = np.einsum("ij,ij->i", p, class_neighbours)

    proba = np.copy(p)
    if loss == 'l1':
        proba *= p_i[:, np.newaxis]
    mask = proba > threshold
    proba *= mask

    grad = np.einsum("ij,ijkl", proba, outer)
    assert grad.shape == (n_features, n_features)

    neighbours_proba = class_neighbours * p
    if loss == 'kl':
        nonzero = p_i > 1e-10
        neighbours_proba[nonzero, :] /= p_i[nonzero, np.newaxis]

    neighbours_mask = neighbours_proba > threshold
    neighbours_proba *= neighbours_mask

    neighbours_term = np.einsum("ij,ijkl", neighbours_proba, outer)

    assert neighbours_term.shape == (n_features, n_features)
    grad -= neighbours_term

    function_value = 0
    if loss == 'l1':
        function_value = p_i.sum()
    elif loss == 'kl':
        function_value = np.log(p_i).sum()

    grad = 2 * np.dot(L, grad)

    return [-function_value, -grad.flatten()]


def nca_semivectorized_oracle(L, X, y, n_components, loss, threshold=0.0):
    n_samples, n_features = X.shape
    L = L.reshape((n_components, n_features))

    Lx = np.dot(X, L.T)  # n_samples x n_components
    assert Lx.shape == (n_samples, n_components)

    grad = np.zeros((n_features, n_features))
    function_value = 0

    p = -euclidean_distances(Lx, squared=True)
    np.fill_diagonal(p, -np.inf)
    p -= sp.misc.logsumexp(p, axis=1)[:, np.newaxis]
    np.exp(p, out=p)

    for i in range(n_samples):
        class_neighbours = y == y[i]
        p_i = (p[i] * class_neighbours).sum()

        samples_proba = np.copy(p[i])
        if loss == 'l1':
            samples_proba *= p_i
        mask = samples_proba > threshold  # n_samples

        Xij = X[i] - X[mask, :]  # n_relevant x n_features
        Xij_outer = (Xij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
        grad += np.einsum("i,ijk", samples_proba[mask], Xij_outer)

        neighbours_proba = p[i] * (y == y[i])
        if loss == 'kl' and p_i > 1e-10:
            neighbours_proba /= p_i
        mask = neighbours_proba > threshold  # n_samples

        Xij = X[i] - X[mask, :]  # n_relevant x n_features
        Xij_outer = (Xij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
        grad -= np.einsum("i,ijk", neighbours_proba[mask], Xij_outer)

        if loss == 'l1':
            function_value += p_i
        elif loss == 'kl':
            if p_i > 1e-10:
                function_value += np.log(p_i)
            else:
                # FIXME
                function_value += 0

    grad = 2 * np.dot(L, grad)

    return [-function_value, -grad.flatten()]


def nca_stochastic_oracle(L, X, y, n_components, loss, minibatch_size, threshold=0.0):
    n_samples, n_features = X.shape
    L = L.reshape((n_components, n_features))

    Lx = np.dot(X, L.T)  # n_samples x n_components
    assert Lx.shape == (n_samples, n_components)

    grad = np.zeros((n_features, n_features))
    function_value = 0


    samples = np.arange(n_samples)
    np.random.shuffle(samples)  # FIXME use seeded rng
    for i in samples[:minibatch_size]:
        A = Lx[i] - Lx
        p = -np.einsum("ij,ij->i", A, A)
        p[i] = -np.inf
        p -= sp.misc.logsumexp(p)
        np.exp(p, out=p)
        assert p.shape == (n_samples,)

        class_neighbours = y == y[i]
        p_i = (p * class_neighbours).sum()

        samples_proba = np.copy(p)
        if loss == 'l1':
            samples_proba *= p_i
        mask = samples_proba > threshold  # n_samples

        Xij = X[i] - X[mask, :]  # n_relevant x n_features
        Xij_outer = (Xij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
        grad += np.einsum("i,ijk", samples_proba[mask], Xij_outer)

        neighbours_proba = p * (y == y[i])
        if loss == 'kl' and p_i > 1e-10:
            neighbours_proba /= p_i
        mask = neighbours_proba > threshold  # n_samples

        Xij = X[i] - X[mask, :]  # n_relevant x n_features
        Xij_outer = (Xij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
        grad -= np.einsum("i,ijk", neighbours_proba[mask], Xij_outer)

        if loss == 'l1':
            function_value += p_i
        elif loss == 'kl':
            if p_i > 1e-10:
                function_value += np.log(p_i)
            else:
                # FIXME
                function_value += 0

    grad = 2 * np.dot(L, grad)

    return [-function_value, -grad.flatten()]


def optimize_nca(X, y, learning_rate, n_components, loss, n_init, max_iter,
                 solver, tol, random_state, method, minibatch_size, verbose):
    n_samples, n_features = X.shape

    rng = np.random.RandomState(random_state)
    best_value = np.inf
    best_L = None

    threshold = 1e-3 / np.abs(X).max() ** 2  # FIXME

    if method == 'vectorized':
        dx = X[:, np.newaxis, :] - X[np.newaxis, :, :]  # n_samples x n_samples x n_features
        outer = dx[:, :, np.newaxis, :] * dx[:, :, :, np.newaxis]  # n_samples x n_samples x n_features x n_features
        extra_args = (X, y, n_components, loss, outer, threshold)
        oracle = nca_vectorized_oracle
    elif method == 'semivectorized':
        extra_args = (X, y, n_components, loss, threshold)
        oracle = nca_semivectorized_oracle
    elif method == 'stochastic':
        extra_args = (X, y, n_components, loss, minibatch_size, threshold)
        oracle = nca_stochastic_oracle

    for _ in range(n_init):
        start = time.clock()
        L = rng.uniform(0, 1, (n_components, n_features))

        if solver in ['gd', 'adagrad']:
            is_adagrad = solver == 'adagrad'
            grad_history = 1e-3 * np.ones_like(L) if is_adagrad else None

            for it in range(max_iter):
                value, grad = oracle(L, *extra_args)
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

        elif solver == 'scipy':
            state = {'it': 0}

            def callback(L):
                state['it'] += 1
                if verbose > 1:
                    value, _ = oracle(L, *extra_args)
                    percent = 100. * (state['it'] + 1) / max_iter
                    sys.stdout.write("\rProgress {:.2f}% :: Target = {}".format(percent, value))
                    sys.stdout.flush()

            options = {'maxiter': max_iter, "disp": verbose > 1, "gtol": tol}
            res = sp.optimize.minimize(fun=oracle, x0=L, args=extra_args, jac=True, options=options, callback=callback)
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
    def __init__(self, n_components, solver, learning_rate, tol, loss, max_iter,
                 n_init, random_state, verbose, method):
        self.n_components = n_components
        self.random_state = random_state
        self.max_iter = max_iter
        self.learning_rate = learning_rate  # TODO: learning rate scheduling from sgd_fast?
        self.n_init = n_init
        self.loss = loss
        self.verbose = verbose
        self.solver = solver
        self.tol = tol
        self.method = method

    def fit(self, X, y):
        n_features = X.shape[1]

        if self.solver not in ['adagrad', 'gd', 'scipy']:
            raise ValueError("Unsupported solver: %s" % self.solver)

        # TODO: more validation

        n_components = self.n_components or n_features

        self.matrix_ = optimize_nca(X, y, self.learning_rate, n_components, self.loss,
                                    self.n_init, self.max_iter, self.solver, self.tol,
                                    self.random_state, self.method, self.verbose)
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

    solver: string

    learning_rate: float

    tol: float

    loss: string

    max_iter: int

    n_init: int

    random_state: int seed, RandomState instance, or None (default)

    method: 'vectorized' | 'semivectorized' | 'stochastic'

    verbose: int

    """
    def __init__(self, n_components=None, solver='adagrad', learning_rate=1.0, tol=1e-5, loss='kl',
                 max_iter=100, n_init=1, random_state=None, verbose=0, method='stochastic', minibatch_size=100):
        super(NCATransformer, self).__init__(n_components, solver, learning_rate, tol, loss,
                                             max_iter, n_init, random_state, verbose, method, minibatch_size)

    def transform(self, X):
        check_is_fitted(self, 'matrix_')
        return np.dot(X, self.matrix_.T)

# TODO: add NCASimilarity
