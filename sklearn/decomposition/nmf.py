""" Non-negative matrix factorization
"""
# Author: Vlad Niculae
#         Lars Buitinck <L.J.Buitinck@uva.nl>
# Author: Chih-Jen Lin, National Taiwan University (original projected gradient
#     NMF implementation)
# Author: Anthony Di Franco (original Python and NumPy port)
# License: BSD 3 clause


from __future__ import division

from math import sqrt
import warnings
import numbers

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin
from ..utils import atleast2d_or_csr, check_random_state
from ..utils.extmath import randomized_svd, safe_sparse_dot


def safe_vstack(Xs):
    if any(sp.issparse(X) for X in Xs):
        return sp.vstack(Xs)
    else:
        return np.vstack(Xs)


def norm(x):
    """Dot product-based Euclidean norm implementation

    See: http://fseoane.net/blog/2011/computing-the-vector-norm/
    """
    x = x.ravel()
    return np.sqrt(np.dot(x.T, x))


def _sparseness(x):
    """Hoyer's measure of sparsity for a vector"""
    sqrt_n = np.sqrt(len(x))
    return (sqrt_n - np.linalg.norm(x, 1) / norm(x)) / (sqrt_n - 1)


def check_non_negative(X, whom):
    X = X.data if sp.issparse(X) else X
    if (X < 0).any():
        raise ValueError("Negative values in data passed to %s" % whom)


def _initialize_nmf(X, n_components, variant=None, eps=1e-6,
                    random_state=None):
    """NNDSVD algorithm for NMF initialization.

    Computes a good initial guess for the non-negative
    rank k matrix approximation for X: X = WH

    Parameters
    ----------

    X : array, [n_samples, n_features]
        The data matrix to be decomposed.

    n_components : array, [n_components, n_features]
        The number of components desired in the approximation.

    variant : None | 'a' | 'ar'
        The variant of the NNDSVD algorithm.
        Accepts None, 'a', 'ar'
        None: leaves the zero entries as zero
        'a': Fills the zero entries with the average of X
        'ar': Fills the zero entries with standard normal random variates.
        Default: None

    eps: float
        Truncate all values less then this in output to zero.

    random_state : numpy.RandomState | int, optional
        The generator used to fill in the zeros, when using variant='ar'
        Default: numpy.random

    Returns
    -------

    (W, H) :
        Initial guesses for solving X ~= WH such that
        the number of columns in W is n_components.

    Remarks
    -------

    This implements the algorithm described in
    C. Boutsidis, E. Gallopoulos: SVD based
    initialization: A head start for nonnegative
    matrix factorization - Pattern Recognition, 2008

    http://tinyurl.com/nndsvd
    """
    check_non_negative(X, "NMF initialization")
    if variant not in (None, 'a', 'ar'):
        raise ValueError("Invalid variant name")

    U, S, V = randomized_svd(X, n_components)
    W, H = np.zeros(U.shape), np.zeros(V.shape)

    # The leading singular triplet is non-negative
    # so it can be used as is for initialization.
    W[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
    H[0, :] = np.sqrt(S[0]) * np.abs(V[0, :])

    for j in range(1, n_components):
        x, y = U[:, j], V[j, :]

        # extract positive and negative parts of column vectors
        x_p, y_p = np.maximum(x, 0), np.maximum(y, 0)
        x_n, y_n = np.abs(np.minimum(x, 0)), np.abs(np.minimum(y, 0))

        # and their norms
        x_p_nrm, y_p_nrm = norm(x_p), norm(y_p)
        x_n_nrm, y_n_nrm = norm(x_n), norm(y_n)

        m_p, m_n = x_p_nrm * y_p_nrm, x_n_nrm * y_n_nrm

        # choose update
        if m_p > m_n:
            u = x_p / x_p_nrm
            v = y_p / y_p_nrm
            sigma = m_p
        else:
            u = x_n / x_n_nrm
            v = y_n / y_n_nrm
            sigma = m_n

        lbd = np.sqrt(S[j] * sigma)
        W[:, j] = lbd * u
        H[j, :] = lbd * v

    W[W < eps] = 0
    H[H < eps] = 0

    if variant == "a":
        avg = X.mean()
        W[W == 0] = avg
        H[H == 0] = avg
    elif variant == "ar":
        random_state = check_random_state(random_state)
        avg = X.mean()
        W[W == 0] = abs(avg * random_state.randn(len(W[W == 0])) / 100)
        H[H == 0] = abs(avg * random_state.randn(len(H[H == 0])) / 100)

    return W, H


def _nls_lbfgs(X, Y, init=None, tol=1e-3, max_iter=100):
    """Non-negative least squares solver using L-BFGS.

    Solves for w in Xw = Y
    """
    X, Y = check_arrays(sparse_format='csr')
    n_samples, n_features = X.shape
    n_targets = Y.shape[1]
    last_grad = np.empty((n_features))
    G = safe_sparse_dot(X.T, X)
    Xy = safe_sparse_dot(X.T, Y)

    def f(w, *args):
        return 0.5 * np.sum((safe_sparse_dot(X, w) - Y[:, k]) ** 2)

    def fprime(w, last_grad, *args):
        last_grad[:] = np.dot(G, w) - Xy[:, k]
        # last_grad[:] = safe_sparse_dot(X.T, safe_sparse_dot(X, w) - y)
        return last_grad

    if init is None:
        init = np.zeros(n_features, dtype=np.float64)
    w, _, d = fmin_l_bfgs_b(f, x0=init, fprime=fprime, pgtol=tol,
                            bounds=[(0, None)] * n_features,
                            maxiter=max_iter,
                            args=[last_grad],
                            disp=0)

    return w, last_grad, d['nit']


class ProjectedGradientNMF(BaseEstimator, TransformerMixin):
    """
    Non-negative Matrix Factorization by Alternating NNLS.

    Parameters
    ----------
    n_components : int or None
        Number of components, if n_components is not set all components
        are kept

    init :  'nndsvd' |  'nndsvda' | 'nndsvdar' | 'random'
        Method used to initialize the procedure.
        Default: 'nndsvdar' if n_components < n_features, otherwise random.
        Valid options::

            'nndsvd': Nonnegative Double Singular Value Decomposition (NNDSVD)
                initialization (better for sparseness)
            'nndsvda': NNDSVD with zeros filled with the average of X
                (better when sparsity is not desired)
            'nndsvdar': NNDSVD with zeros filled with small random values
                (generally faster, less accurate alternative to NNDSVDa
                for when sparsity is not desired)
            'random': non-negative random matrices

    nnls_solver : {'pg', 'lbfgs'}
        Non-negative least squares solver to use for learning.  Note that
        L-bfgs is always used for transforming.

    sparseness : 'data' | 'components' | None, default: None
        Where to enforce sparsity in the model.

    beta : double, default: 1
        Degree of sparseness, if sparseness is not None. Larger values mean
        more sparseness.

    eta : double, default: 0.1
        Degree of correctness to maintain, if sparsity is not None. Smaller
        values mean larger error.

    tol : double, default: 1e-4
        Tolerance value used in stopping conditions.

    max_iter : int, default: 200
        Number of iterations to compute.

    nls_max_iter : int, default: 2000
        Number of iterations in NLS subproblem.

    random_state : int or RandomState
        Random number generator seed control.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Non-negative components of the data.

    `reconstruction_err_` : number
        Frobenius norm of the matrix difference between
        the training data and the reconstructed data from
        the fit produced by the model. ``|| X - WH ||_2``

    Examples
    --------

    >>> import numpy as np
    >>> X = np.array([[1,1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from sklearn.decomposition import ProjectedGradientNMF
    >>> model = ProjectedGradientNMF(n_components=2, init='random',
    ...                              random_state=0)
    >>> model.fit(X) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ProjectedGradientNMF(beta=1, eta=0.1, init='random', max_iter=200,
            n_components=2, nls_max_iter=2000, random_state=0, sparseness=None,
            tol=0.0001)
    >>> model.components_
    array([[ 0.77032744,  0.11118662],
           [ 0.38526873,  0.38228063]])
    >>> model.reconstruction_err_ #doctest: +ELLIPSIS
    0.00746...
    >>> model = ProjectedGradientNMF(n_components=2,
    ...              sparseness='components', init='random', random_state=0)
    >>> model.fit(X) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ProjectedGradientNMF(beta=1, eta=0.1, init='random', max_iter=200,
                n_components=2, nls_max_iter=2000, random_state=0,
                sparseness='components', tol=0.0001)
    >>> model.components_
    array([[ 1.67481991,  0.29614922],
           [ 0.        ,  0.4681982 ]])
    >>> model.reconstruction_err_ #doctest: +ELLIPSIS
    0.513...

    References
    ----------
    This implements

    C.-J. Lin. Projected gradient methods
    for non-negative matrix factorization. Neural
    Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    P. Hoyer. Non-negative Matrix Factorization with
    Sparseness Constraints. Journal of Machine Learning
    Research 2004.

    NNDSVD is introduced in

    C. Boutsidis, E. Gallopoulos: SVD based
    initialization: A head start for nonnegative
    matrix factorization - Pattern Recognition, 2008
    http://tinyurl.com/nndsvd
    """

    def __init__(self, n_components=None, init=None, nnls_solver='pg',
                 sparseness=None, beta=1, eta=0.1, tol=1e-4, max_iter=200,
                 nls_max_iter=2000, adaptive=True, random_state=None):
        self.n_components = n_components
        self.init = init
        self.tol = tol
        self.nnls_solver = nnls_solver
        self.sparseness = sparseness
        self.beta = beta
        self.eta = eta
        self.max_iter = max_iter
        self.nls_max_iter = nls_max_iter
        self.adaptive = adaptive
        self.random_state = random_state

    def _init(self, X):
        n_samples, n_features = X.shape
        init = self.init
        if init is None:
            if self.n_components_ < n_features:
                init = 'nndsvd'
            else:
                init = 'random'

        if isinstance(init, (numbers.Integral, np.random.RandomState)):
            random_state = check_random_state(init)
            init = "random"
            warnings.warn("Passing a random seed or generator as init "
                          "is deprecated and will be removed in 0.15. Use "
                          "init='random' and random_state instead.",
                          DeprecationWarning)
        else:
            random_state = self.random_state

        if init == 'nndsvd':
            W, H = _initialize_nmf(X, self.n_components_)
        elif init == 'nndsvda':
            W, H = _initialize_nmf(X, self.n_components_, variant='a')
        elif init == 'nndsvdar':
            W, H = _initialize_nmf(X, self.n_components_, variant='ar')
        elif init == "random":
            rng = check_random_state(random_state)
            W = rng.randn(n_samples, self.n_components_)
            # we do not write np.abs(W, out=W) to stay compatible with
            # numpy 1.5 and earlier where the 'out' keyword is not
            # supported as a kwarg on ufuncs
            np.abs(W, W)
            H = rng.randn(self.n_components_, n_features)
            np.abs(H, H)
        else:
            raise ValueError(
                'Invalid init parameter: got %r instead of one of %r' %
                (init, (None, 'nndsvd', 'nndsvda', 'nndsvdar', 'random')))
        return W, H

    def _update_W(self, X, H, W, tolW):
        n_samples, n_features = X.shape

        if self.sparseness is None:
            W, gradW, iterW = self._nnls(H.T, X.T, W.T, tolW,
                                         self.nls_max_iter)
        elif self.sparseness == 'data':
            W, gradW, iterW = self._nnls(
                safe_vstack([H.T, np.sqrt(self.beta) * np.ones((1,
                             self.n_components_))]),
                safe_vstack([X.T, np.zeros((1, n_samples))]),
                W.T, tolW, self.nls_max_iter)
        elif self.sparseness == 'components':
            W, gradW, iterW = self._nnls(
                safe_vstack([H.T,
                             np.sqrt(self.eta) * np.eye(self.n_components_)]),
                safe_vstack([X.T,
                             np.zeros((self.n_components_, n_samples))]),
                W.T, tolW, self.nls_max_iter)

        return W, gradW, iterW

    def _update_H(self, X, H, W, tolH):
        n_samples, n_features = X.shape
        if self.sparseness is None:
            H, gradH, iterH = self._nnls(W, X, H, tolH,
                                         self.nls_max_iter)
        elif self.sparseness == 'data':
            H, gradH, iterH = self._nnls(
                safe_vstack([W,
                             np.sqrt(self.eta) * np.eye(self.n_components_)]),
                safe_vstack([X, np.zeros((self.n_components_, n_features))]),
                H, tolH, self.nls_max_iter)

        elif self.sparseness == 'components':
            H, gradH, iterH = self._nnls(
                safe_vstack([W,
                             np.sqrt(self.beta)
                             * np.ones((1, self.n_components_))]),
                safe_vstack([X, np.zeros((1, n_features))]),
                H, tolH, self.nls_max_iter)

        return H, gradH, iterH

    def fit_transform(self, X, y=None):
        """Learn a NMF model for the data X and returns the transformed data.

        This is more efficient than calling fit followed by transform.

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be decomposed

        Returns
        -------
        data: array, [n_samples, n_components]
            Transformed data
        """
        X = atleast2d_or_csr(X)
        check_non_negative(X, "NMF.fit")

        if self.sparseness not in (None, 'data', 'components'):
            raise ValueError(
                'Invalid sparseness parameter: got %r instead of one of %r' %
                (self.sparseness, (None, 'data', 'components')))

        n_samples, n_features = X.shape

        if not self.n_components:
            self.n_components_ = n_features
        else:
            self.n_components_ = self.n_components

        W, H = self._init(X)

        gradW = (np.dot(W, np.dot(H, H.T))
                 - safe_sparse_dot(X, H.T, dense_output=True))
        gradH = (np.dot(np.dot(W.T, W), H)
                 - safe_sparse_dot(W.T, X, dense_output=True))
        init_grad = sqrt(norm(gradW) ** 2 + norm(gradH.T) ** 2)

        if self.adaptive:
            tolW = tolH = max(0.001, self.tol) * init_grad  # why max?
        else:
            tolW = tolH = self.tol

        for n_iter in range(1, self.max_iter + 1):
            # stopping condition
            # as discussed in paper
            proj_norm = sqrt(
                norm(gradW[np.logical_or(gradW < 0, W > 0)]) ** 2 +
                norm(gradH[np.logical_or(gradH < 0, H > 0)]) ** 2)
            if proj_norm < self.tol * init_grad:
                break

            # update W
            W, gradW, iterW = self._update_W(X, H, W, tolW)

            W = W.T
            gradW = gradW.T

            if self.adaptive and iterW == 1:
                tolW = 0.1 * tolW

            # update H
            H, gradH, iterH = self._update_H(X, H, W, tolH)

            if self.adaptive and iterH == 1:
                tolH = 0.1 * tolH

            self.comp_sparseness_ = _sparseness(H.ravel())
            self.data_sparseness_ = _sparseness(W.ravel())

            if not sp.issparse(X):
                self.reconstruction_err_ = norm(X - np.dot(W, H))
            else:
                norm2X = np.sum(X.data ** 2)  # Ok because X is CSR
                normWHT = np.trace(np.dot(np.dot(H.T, np.dot(W.T, W)), H))
                cross_prod = np.trace(np.dot((X * H.T).T, W))
                self.reconstruction_err_ = sqrt(norm2X + normWHT
                                                - 2. * cross_prod)

            self.components_ = H

        if n_iter == self.max_iter:
            warnings.warn("Iteration limit reached during fit")

        return W

    def fit(self, X, y=None, **params):
        """Learn a NMF model for the data X.

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be decomposed

        Returns
        -------
        self
        """
        self.fit_transform(X, **params)
        return self

    def transform(self, X):
        """Transform the data X according to the fitted NMF model

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be transformed by the model

        Returns
        -------
        data: array, [n_samples, n_components]
            Transformed data
        """
        X = atleast2d_or_csr(X)
        W = np.zeros((X.shape[0], self.n_components_))
        W, _, _ = self._nnls_lbfgs(self.components_.T, X.T, W.T, tol=1e-7)
        return W.T

    def _nnls(self, X, Y, W, tol=1e-3, max_iter=100, *args):
        if self.nnls_solver == 'pg':
            return self._nnls_pg(X, Y, W, tol, max_iter, *args)
        elif self.nnls_solver == 'lbfgs':
            return self._nnls_lbfgs(X, Y, W, tol, max_iter, *args)
        else:
            raise ValueError('Unsupported solver: %s' % self.nnls_solver)

    def _nnls_pg(self, X, Y, W, tol=1e-3, max_iter=100, sigma=0.01, beta=0.1):
        """Non-negative least square solver

        Solves a non-negative least squares subproblem using the
        projected gradient descent algorithm.
        min || XW - Y ||_2

        Parameters
        ----------
        V, W : array-like
            Constant matrices.

        H_init : array-like
            Initial guess for the solution.

        tol : float
            Tolerance of the stopping condition.

        max_iter : int
            Maximum number of iterations before timing out.

        sigma : float
            Constant used in the sufficient decrease condition checked by the
            line search.  Smaller values lead to a looser sufficient decrease
            condition, thus reducing the time taken by the line search, but
            potentially increasing the number of iterations of the projected
            gradient procedure. 0.01 is a commonly used value in the
            optimization literature.

        beta : float
            Factor by which the step size is decreased (resp. increased) until
            (resp. as long as) the sufficient decrease condition is satisfied.
            Larger values allow to find a better step size but lead to longer
            line search. 0.1 is a commonly used value in the optimization
            literature.

        Returns
        -------
        W : array-like
            Solution to the non-negative least squares problem.

        grad : array-like
            The gradient.

        n_iter : int
            The number of iterations done by the algorithm.

        Reference
        ---------

        C.-J. Lin. Projected gradient methods
        for non-negative matrix factorization. Neural
        Computation, 19(2007), 2756-2779.
        http://www.csie.ntu.edu.tw/~cjlin/nmf/

        """
        if (W < 0).any():
            raise ValueError("Negative values in H_init passed to NLS solver.")

        Xy = safe_sparse_dot(X.T, Y, dense_output=True)
        G = safe_sparse_dot(X.T, X, dense_output=True)

        # values justified in the paper
        alpha = 1
        for n_iter in range(1, max_iter + 1):
            grad = np.dot(G, W) - Xy
            proj_gradient = norm(grad[np.logical_or(grad < 0, W > 0)])
            if proj_gradient < tol:
                break

            for inner_iter in range(1, 20):
                # Gradient step.
                Wn = W - alpha * grad
                # Projection step.
                Wn = np.maximum(Wn, 0)
                d = Wn - W
                gradd = np.sum(grad * d)
                dGd = np.sum(np.dot(G, d) * d)
                suff_decr = (1 - sigma) * gradd + 0.5 * dGd < 0
                if inner_iter == 1:
                    decr_alpha = not suff_decr
                    Wp = W

                if decr_alpha:
                    if suff_decr:
                        W = Wn
                        break
                    else:
                        alpha *= beta
                elif not suff_decr or (Wp == Wn).all():
                    W = Wp
                    break
                else:
                    alpha /= beta
                    Wp = Wn

        if n_iter == max_iter:
            warnings.warn("Iteration limit reached in nls subproblem.")

        return W, grad, n_iter

    def _nnls_lbfgs(self, X, Y, W=None, tol=1e-3, max_iter=100):
        n_samples, n_features = X.shape
        n_targets = Y.shape[1]
        n_iter = np.empty(n_targets, dtype=np.int)
        if W is None:
            W = np.zeros((n_features, n_targets))
        grad = np.empty((n_features, n_targets))
        for k in xrange(n_targets):
            W[:, k], grad[:, k], n_iter[k] = _nls_lbfgs(X, Y[:, k],
                                                        init=W[:, k],
                                                        tol=tol)
        return W, grad, n_iter.max() + 1


class NMF(ProjectedGradientNMF):
    __doc__ = ProjectedGradientNMF.__doc__
    pass
