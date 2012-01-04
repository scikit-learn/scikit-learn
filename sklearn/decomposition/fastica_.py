"""
Python implementation of the fast ICA algorithms.

Reference: Tables 8.3 and 8.4 page 196 in the book:
Independent Component Analysis, by  Hyvarinen et al.
"""

# Author: Pierre Lafaye de Micheaux, Stefan van der Walt, Gael Varoquaux,
#         Bertrand Thirion, Alexandre Gramfort
# License: BSD 3 clause
import numpy as np
from scipy import linalg

from ..base import BaseEstimator
from ..utils import array2d, as_float_array

__all__ = ['fastica', 'FastICA']


def _gs_decorrelation(w, W, j):
    """
    Orthonormalize w wrt the first j rows of W

    Parameters
    ----------
    w: array of shape(n), to be orthogonalized
    W: array of shape(p, n), null space definition
    j: int < p

    caveats
    -------
    assumes that W is orthogonal
    w changed in place
    """
    w -= np.dot(np.dot(w, W[:j].T), W[:j])
    return w


def _sym_decorrelation(W):
    """ Symmetric decorrelation
    i.e. W <- (W * W.T) ^{-1/2} * W
    """
    K = np.dot(W, W.T)
    s, u = linalg.eigh(K)
    # u (resp. s) contains the eigenvectors (resp. square roots of
    # the eigenvalues) of W * W.T
    W = np.dot(np.dot(np.dot(u, np.diag(1.0 / np.sqrt(s))), u.T), W)
    return W


def _ica_def(X, tol, g, gprime, fun_args, max_iter, w_init):
    """Deflationary FastICA using fun approx to neg-entropy function

    Used internally by FastICA.
    """

    n_components = w_init.shape[0]
    W = np.zeros((n_components, n_components), dtype=float)

    # j is the index of the extracted component
    for j in range(n_components):
        w = w_init[j, :].copy()
        w /= np.sqrt((w ** 2).sum())

        n_iterations = 0
        # we set lim to tol+1 to be sure to enter at least once in next while
        lim = tol + 1
        while ((lim > tol) & (n_iterations < (max_iter - 1))):
            wtx = np.dot(w.T, X)
            gwtx = g(wtx, fun_args)
            g_wtx = gprime(wtx, fun_args)
            w1 = (X * gwtx).mean(axis=1) - g_wtx.mean() * w

            _gs_decorrelation(w1, W, j)

            w1 /= np.sqrt((w1 ** 2).sum())

            lim = np.abs(np.abs((w1 * w).sum()) - 1)
            w = w1
            n_iterations = n_iterations + 1

        W[j, :] = w

    return W


def _ica_par(X, tol, g, gprime, fun_args, max_iter, w_init):
    """Parallel FastICA.

    Used internally by FastICA --main loop

    """
    n, p = X.shape

    W = _sym_decorrelation(w_init)

    # we set lim to tol+1 to be sure to enter at least once in next while
    lim = tol + 1
    it = 0
    while ((lim > tol) and (it < (max_iter - 1))):
        wtx = np.dot(W, X)
        gwtx = g(wtx, fun_args)
        g_wtx = gprime(wtx, fun_args)
        W1 = np.dot(gwtx, X.T) / float(p) \
             - np.dot(np.diag(g_wtx.mean(axis=1)), W)

        W1 = _sym_decorrelation(W1)

        lim = max(abs(abs(np.diag(np.dot(W1, W.T))) - 1))
        W = W1
        it += 1

    return W


def fastica(X, n_components=None, algorithm="parallel", whiten=True,
            fun="logcosh", fun_prime='', fun_args={}, max_iter=200,
            tol=1e-04, w_init=None):
    """Perform Fast Independent Component Analysis.

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.
    n_components : int, optional
        Number of components to extract. If None no dimension reduction
        is performed.
    algorithm : {'parallel', 'deflation'}, optional
        Apply a parallel or deflational FASTICA algorithm.
    whiten: boolean, optional
        If true perform an initial whitening of the data. Do not set to
        false unless the data is already white, as you will get incorrect
        results.
        If whiten is true, the data is assumed to have already been
        preprocessed: it should be centered, normed and white.
    fun : string or function, optional
        The functional form of the G function used in the
        approximation to neg-entropy. Could be either 'logcosh', 'exp',
        or 'cube'.
        You can also provide your own function but in this case, its
        derivative should be provided via argument fun_prime
    fun_prime : empty string ('') or function, optional
        See fun.
    fun_args: dictionary, optional
        If empty and if fun='logcosh', fun_args will take value
        {'alpha' : 1.0}
    max_iter: int, optional
        Maximum number of iterations to perform
    tol: float, optional
        A positive scalar giving the tolerance at which the
        un-mixing matrix is considered to have converged
    w_init: (n_components, n_components) array, optional
        Initial un-mixing array of dimension (n.comp,n.comp).
        If None (default) then an array of normal r.v.'s is used
    source_only: boolean, optional
        if True, only the sources matrix is returned

    Returns
    -------
    K: (n_components, p) array or None.
        If whiten is 'True', K is the pre-whitening matrix that projects data
        onto the first n.comp principal components. If whiten is 'False', K is
        'None'.

    W: (n_components, n_components) array
        estimated un-mixing matrix
        The mixing matrix can be obtained by::

            w = np.dot(W, K.T)
            A = w.T * (w * w.T).I

    S: (n_components, n) array
        estimated source matrix


    Notes
    -----

    The data matrix X is considered to be a linear combination of
    non-Gaussian (independent) components i.e. X = AS where columns of S
    contain the independent components and A is a linear mixing
    matrix. In short ICA attempts to `un-mix' the data by estimating an
    un-mixing matrix W where ``S = W K X.``

    This implementation was originally made for data of shape
    [n_features, n_samples]. Now the input is transposed
    before the algorithm is applied. This makes it slightly
    faster for Fortran-ordered input.

    Implemented using FastICA:
    `A. Hyvarinen and E. Oja, Independent Component Analysis:
    Algorithms and Applications, Neural Networks, 13(4-5), 2000,
    pp. 411-430`

    """
    # make interface compatible with other decompositions
    X = array2d(X).T

    algorithm_funcs = {'parallel': _ica_par,
                       'deflation': _ica_def}

    alpha = fun_args.get('alpha', 1.0)
    if (alpha < 1) or (alpha > 2):
        raise ValueError("alpha must be in [1,2]")

    if isinstance(fun, str):
        # Some standard nonlinear functions
        # XXX: these should be optimized, as they can be a bottleneck.
        if fun == 'logcosh':
            def g(x, fun_args):
                alpha = fun_args.get('alpha', 1.0)
                return np.tanh(alpha * x)

            def gprime(x, fun_args):
                alpha = fun_args.get('alpha', 1.0)
                return alpha * (1 - (np.tanh(alpha * x)) ** 2)

        elif fun == 'exp':
            def g(x, fun_args):
                return x * np.exp(-(x ** 2) / 2)

            def gprime(x, fun_args):
                return (1 - x ** 2) * np.exp(-(x ** 2) / 2)

        elif fun == 'cube':
            def g(x, fun_args):
                return x ** 3

            def gprime(x, fun_args):
                return 3 * x ** 2
        else:
            raise ValueError(
                        'fun argument should be one of logcosh, exp or cube')
    elif callable(fun):
        def g(x, fun_args):
            return fun(x, **fun_args)

        def gprime(x, fun_args):
            return fun_prime(x, **fun_args)
    else:
        raise ValueError('fun argument should be either a string '
                         '(one of logcosh, exp or cube) or a function')

    n, p = X.shape

    if n_components is None:
        n_components = min(n, p)
    if (n_components > min(n, p)):
        n_components = min(n, p)
        print("n_components is too large: it will be set to %s" % n_components)

    if whiten:
        # Centering the columns (ie the variables)
        X = X - X.mean(axis=-1)[:, np.newaxis]

        # Whitening and preprocessing by PCA
        u, d, _ = linalg.svd(X, full_matrices=False)

        del _
        K = (u / d).T[:n_components]  # see (6.33) p.140
        del u, d
        X1 = np.dot(K, X)
        # see (13.6) p.267 Here X1 is white and data
        # in X has been projected onto a subspace by PCA
    else:
        # X must be casted to floats to avoid typing issues with numpy
        # 2.0 and the line below
        X1 = as_float_array(X, copy=True)
    X1 *= np.sqrt(p)

    if w_init is None:
        w_init = np.random.normal(size=(n_components, n_components))
    else:
        w_init = np.asarray(w_init)
        if w_init.shape != (n_components, n_components):
            raise ValueError("w_init has invalid shape -- should be %(shape)s"
                             % {'shape': (n_components, n_components)})

    kwargs = {'tol': tol,
              'g': g,
              'gprime': gprime,
              'fun_args': fun_args,
              'max_iter': max_iter,
              'w_init': w_init}

    func = algorithm_funcs.get(algorithm, 'parallel')

    W = func(X1, **kwargs)
    del X1

    if whiten:
        S = np.dot(np.dot(W, K), X)
        return K, W, S.T
    else:
        S = np.dot(W, X)
        return None, W, S.T


class FastICA(BaseEstimator):
    """FastICA; a fast algorithm for Independent Component Analysis

    Parameters
    ----------
    n_components : int, optional
        Number of components to use. If none is passed, all are used.
    algorithm : {'parallel', 'deflation'}
        Apply parallel or deflational algorithm for FastICA
    whiten : boolean, optional
        If whiten is false, the data is already considered to be
        whitened, and no whitening is performed.
    fun : {'logcosh', 'exp', or 'cube'}, or a callable
        The non-linear function used in the FastICA loop to approximate
        negentropy. If a function is passed, it derivative should be
        passed as the 'fun_prime' argument.
    fun_prime : None or a callable
        The derivative of the non-linearity used.
    max_iter : int, optional
        Maximum number of iterations during fit
    tol : float, optional
        Tolerance on update at each iteration
    w_init : None of an (n_components, n_components) ndarray
        The mixing matrix to be used to initialize the algorithm.

    Attributes
    ----------
    `unmixing_matrix_` : 2D array, [n_components, n_samples]
        The unmixing matrix

    Notes
    -----

    Implementation based on
    `A. Hyvarinen and E. Oja, Independent Component Analysis:
    Algorithms and Applications, Neural Networks, 13(4-5), 2000,
    pp. 411-430`

    """

    def __init__(self, n_components=None, algorithm='parallel', whiten=True,
                 fun='logcosh', fun_prime='', fun_args=None, max_iter=200,
                 tol=1e-4, w_init=None):
        super(FastICA, self).__init__()
        self.n_components = n_components
        self.algorithm = algorithm
        self.whiten = whiten
        self.fun = fun
        self.fun_prime = fun_prime
        self.fun_args = {} if fun_args is None else fun_args
        self.max_iter = max_iter
        self.tol = tol
        self.w_init = w_init

    def fit(self, X):
        whitening_, unmixing_, sources_ = fastica(X, self.n_components,
                        self.algorithm, self.whiten,
                        self.fun, self.fun_prime, self.fun_args, self.max_iter,
                        self.tol, self.w_init)
        if self.whiten == True:
            self.unmixing_matrix_ = np.dot(unmixing_, whitening_)
        else:
            self.unmixing_matrix_ = unmixing_
        self.components_ = sources_
        return self

    def transform(self, X):
        """Apply un-mixing matrix "W" to X to recover the sources

        S = X * W.T
        """
        return np.dot(X, self.unmixing_matrix_.T)

    def get_mixing_matrix(self):
        """Compute the mixing matrix
        """
        return linalg.pinv(self.unmixing_matrix_)
