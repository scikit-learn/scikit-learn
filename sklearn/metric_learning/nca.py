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
        mask = p_i > 0
        function_value = np.log(p_i[mask]).sum()

    grad = 2 * np.dot(L, grad)

    return [-function_value, -grad.flatten()]


def nca_semivectorized_oracle(L, X, y, n_components, loss, threshold=0.0):
    n_samples, n_features = X.shape
    L = L.reshape((n_components, n_features))

    Lx = np.dot(X, L.T)  # n_samples x n_components
    assert Lx.shape == (n_samples, n_components)

    grad = np.zeros((n_components, n_features))
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

        Lxij = Lx[i] - Lx[mask, :]  # n_relevant x n_components
        Xij = X[i] - X[mask, :]  # n_relevant x n_features

        # n_relevant x n_components x n_features
        Xij_outer = Lxij[:, :, np.newaxis] * Xij[:, np.newaxis, :]

        grad += np.einsum("i,ijk", samples_proba[mask], Xij_outer)

        neighbours_proba = p[i] * (y == y[i])
        if loss == 'kl' and p_i > 1e-10:
            neighbours_proba /= p_i
        mask = neighbours_proba > threshold  # n_samples

        Lxij = Lx[i] - Lx[mask, :]
        Xij = X[i] - X[mask, :]  # n_relevant x n_features
        Xij_outer = (Lxij[:, :, np.newaxis] * Xij[:, np.newaxis, :])
        grad -= np.einsum("i,ijk", neighbours_proba[mask], Xij_outer)

        if loss == 'l1':
            function_value += p_i
        elif loss == 'kl':
            if p_i > 1e-10:
                function_value += np.log(p_i)
            else:
                # FIXME
                function_value += 0

    grad *= 2  # np.dot(L, grad)

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
                 solver, tol, random_state, method, minibatch_size, verbose, threshold):
    n_samples, n_features = X.shape

    rng = np.random.RandomState(random_state)
    best_value = np.inf
    best_L = None

    if threshold is None:
        threshold = min(1e-6, 1e-3 / np.abs(X).max() ** 2)  # FIXME

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

    exec_times = []
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
        exec_times.append(exec_time)
        if verbose > 0:
            sys.stdout.write("\rInit {} / {} completed"
                             " :: Time: {:.2f}s :: Loss = {}\n".format(_ + 1, n_init, exec_time, value))
            sys.stdout.flush()

        if value < best_value:
            best_value = value
            best_L = L

    if verbose > 0:
        mean = np.mean(exec_times)
        std = np.std(exec_times)
        print("\rCompleted. Avg. time: {:.2f} +/- {:.2f}s".format(mean, std) + " " * 80)

    return best_L


class BaseNCA(BaseEstimator):
    def __init__(self, n_components, solver, learning_rate, tol, loss, max_iter,
                 n_init, random_state, verbose, method, minitbatch_size, threshold):
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
        self.minitbatch_size = minitbatch_size
        self.threshold = threshold

    def fit(self, X, y):
        n_features = X.shape[1]

        if self.solver not in ['adagrad', 'gd', 'scipy']:
            raise ValueError("Unsupported solver: %s" % self.solver)

        # TODO: more validation

        n_components = self.n_components or n_features

        self.matrix_ = optimize_nca(X, y, self.learning_rate, n_components, self.loss,
                                    self.n_init, self.max_iter, self.solver, self.tol,
                                    self.random_state, self.method, self.minitbatch_size,
                                    self.verbose, self.threshold)
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

    solver: 'adagrad', 'gd' or 'scipy'
        Solver to use for optimization. 'gd' stands for gradient descent,
        'adagrad' stands for gradient descent with AdaGrad heuristic,
        'scipy' means forwarding optimization to scipy.optimization,
        which usually selects LBFGS.

    learning_rate: float
        The learning rate. Used only by gradient descent-based solvers

    tol: float
        Stopping criteria for gradient

    max_iter: int
        Maximal number of iterations for solver

    n_init: int
        Number of random restarts

    random_state: int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        initializing the matrix.

    method: 'vectorized' | 'semivectorized' | 'stochastic'
        Level of vectorization of cost and gradient functions.
        'vectorized' is the fastest, but requires O(N^2 M^2) memory and
        time which doesn't work for medium-size datasets.
        'semivectorized' uses O(N M^2) memory and has the same time
        complexity as the fuly vectorized one.
        'stochastic' is designed to work on really big datasets,
        FIXME!

    minibatch_size: int
        Size of a minibatch for gradient descent-based method.

    verbose: int
        Level of verbosity of debug output. 0 means no output, 1 - basic
        information about initialization runs, 2 - detailed output of each iteration

    """
    def __init__(self, n_components=None, solver='adagrad', learning_rate=1.0, tol=1e-5, loss='kl',
                 max_iter=100, n_init=1, random_state=None, verbose=0, method='stochastic', minibatch_size=100,
                 threshold=None):
        super(NCATransformer, self).__init__(n_components, solver, learning_rate, tol, loss,
                                             max_iter, n_init, random_state, verbose, method,
                                             minibatch_size, threshold)

    def transform(self, X):
        check_is_fitted(self, 'matrix_')
        return np.dot(X, self.matrix_.T)

# TODO: add NCASimilarity
