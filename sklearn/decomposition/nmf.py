""" Non-negative matrix factorization
"""
# Author: Vlad Niculae
#         Lars Buitinck
#         Mathieu Blondel <mathieu@mblondel.org>
#         Tom Dupre la Tour
# Author: Chih-Jen Lin, National Taiwan University (original projected gradient
#                                                   NMF implementation)
# Author: Anthony Di Franco (Projected gradient, Python and NumPy port)
# License: BSD 3 clause


from __future__ import division, print_function

from math import sqrt
import warnings
import numbers
import time

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_random_state, check_array
from ..utils.extmath import randomized_svd, safe_sparse_dot, squared_norm
from ..utils.extmath import fast_dot
from ..utils.validation import check_is_fitted, check_non_negative
from ..utils import deprecated
from ..exceptions import ConvergenceWarning
from .cdnmf_fast import _update_cdnmf_fast

EPSILSON = np.finfo(np.float32).eps


INTEGER_TYPES = (numbers.Integral, np.integer)


def safe_vstack(Xs):
    if any(sp.issparse(X) for X in Xs):
        return sp.vstack(Xs)
    else:
        return np.vstack(Xs)


def norm(x):
    """Dot product-based Euclidean norm implementation

    See: http://fseoane.net/blog/2011/computing-the-vector-norm/
    """
    return sqrt(squared_norm(x))


def trace_dot(X, Y):
    """Trace of np.dot(X, Y.T)."""
    return np.dot(X.ravel(), Y.ravel())


def _sparseness(x):
    """Hoyer's measure of sparsity for a vector"""
    sqrt_n = np.sqrt(len(x))
    return (sqrt_n - np.linalg.norm(x, 1) / norm(x)) / (sqrt_n - 1)


def _check_init(A, shape, whom):
    A = check_array(A)
    if np.shape(A) != shape:
        raise ValueError('Array with wrong shape passed to %s. Expected %s, '
                         'but got %s ' % (whom, shape, np.shape(A)))
    check_non_negative(A, whom)
    if np.max(A) == 0:
        raise ValueError('Array passed to %s is full of zeros.' % whom)


def _safe_compute_error(X, W, H):
    """Frobenius norm between X and WH, safe for sparse array"""
    if not sp.issparse(X):
        error = norm(X - np.dot(W, H))
    else:
        norm_X = np.dot(X.data, X.data)
        norm_WH = trace_dot(np.dot(np.dot(W.T, W), H), H)
        cross_prod = trace_dot((X * H.T), W)
        error = sqrt(norm_X + norm_WH - 2. * cross_prod)
    return error


def beta_divergence(X, W, H, beta):
    """Compute the beta-divergence of X and dot(W, H).

    Parameters
    ----------
    X : float or array-like, shape (n_samples, n_features)

    W : float or dense array-like, shape (n_samples, n_components)

    H : float or dense array-like, shape (n_components, n_features)

    beta : float, string in {'frobenius', 'kullback-leibler', 'itakura-saito'}
        Parameter of the beta-divergence.
        If beta == 2, this is the Frobenius squared norm.
        If beta == 1, this is the generalized Kullback-Leibler divergence.
        If beta == 0, this is the Itakura-Saito divergence.
        Else, this is the general beta-divergence.
    """
    beta = _beta_loss_to_float(beta)

    # The method can be called with scalars
    if not sp.issparse(X):
        X = np.atleast_2d(X)
    W = np.atleast_2d(W)
    H = np.atleast_2d(H)

    # Frobenius norm
    if beta == 2:
        # Avoid the creation of the dense np.dot(W, H) if X is sparse.
        if sp.issparse(X):
            norm_X = np.dot(X.data, X.data)
            norm_WH = trace_dot(np.dot(np.dot(W.T, W), H), H)
            cross_prod = trace_dot((X * H.T), W)
            res = (norm_X + norm_WH - 2. * cross_prod) / 2.
        else:
            res = squared_norm(X - np.dot(W, H)) / 2.
        return res

    if sp.issparse(X):
        # compute np.dot(W, H) only where X is nonzero
        WH_data = _special_sparse_dot(W, H, X).data
        X_data = X.data
    else:
        WH = fast_dot(W, H)
        WH_data = WH.ravel()
        X_data = X.ravel()

    # do not affect the zeros: here 0 ** (-1) = 0 and not infinity
    WH_data = WH_data[X_data != 0]
    X_data = X_data[X_data != 0]

    # used to avoid division by zero
    WH_data[WH_data == 0] = EPSILSON

    # generalized Kullback-Leibler divergence
    if beta == 1:
        # fast and memory efficient computation of np.sum(np.dot(W, H))
        sum_WH = np.dot(np.sum(W, axis=0), np.sum(H, axis=1))
        # computes np.sum(X * log(X / WH)) only where X is nonzero
        div = X_data / WH_data
        res = np.dot(X_data, np.log(div))
        # add full np.sum(np.dot(W, H)) - np.sum(X)
        res += sum_WH - X_data.sum()

    # Itakura-Saito divergence
    elif beta == 0:
        div = X_data / WH_data
        res = np.sum(div) - np.product(X.shape) - np.sum(np.log(div))

    # beta-divergence, beta not in (0, 1, 2)
    else:
        if sp.issparse(X):
            # slow loop, but memory efficient computation of :
            # np.sum(np.dot(W, H) ** beta)
            sum_WH_beta = 0
            for i in range(X.shape[1]):
                sum_WH_beta += np.sum(fast_dot(W, H[:, i]) ** beta)

        else:
            sum_WH_beta = np.sum(WH ** beta)

        sum_X_WH = np.dot(X_data, WH_data ** (beta - 1))
        res = (X_data ** beta).sum() - beta * sum_X_WH
        res += sum_WH_beta * (beta - 1)
        res /= beta * (beta - 1)

    return res


def _special_sparse_dot(W, H, X):
    """Computes np.dot(W, H), only where X is non zero."""
    if sp.issparse(X):
        ii, jj = X.nonzero()
        dot_vals = np.multiply(W[ii, :], H.T[jj, :]).sum(axis=1)
        WH = sp.coo_matrix((dot_vals, (ii, jj)), shape=X.shape)
        return WH.tocsr()
    else:
        return fast_dot(W, H)


def _compute_regularization(alpha, l1_ratio, regularization):
    """Compute L1 and L2 regularization coefficients for W and H"""
    alpha_H = 0.
    alpha_W = 0.
    if regularization in ('both', 'components'):
        alpha_H = float(alpha)
    if regularization in ('both', 'transformation'):
        alpha_W = float(alpha)

    l1_reg_W = alpha_W * l1_ratio
    l1_reg_H = alpha_H * l1_ratio
    l2_reg_W = alpha_W * (1. - l1_ratio)
    l2_reg_H = alpha_H * (1. - l1_ratio)
    return l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H


def _check_string_param(sparseness, solver, regularization, beta_loss, init):
    allowed_sparseness = (None, 'data', 'components')
    if sparseness not in allowed_sparseness:
        raise ValueError(
            'Invalid sparseness parameter: got %r instead of one of %r' %
            (sparseness, allowed_sparseness))

    allowed_solver = ('pg', 'cd', 'mu')
    if solver not in allowed_solver:
        raise ValueError(
            'Invalid solver parameter: got %r instead of one of %r' %
            (solver, allowed_solver))

    allowed_regularization = ('both', 'components', 'transformation', None)
    if regularization not in allowed_regularization:
        raise ValueError(
            'Invalid regularization parameter: got %r instead of one of %r' %
            (regularization, allowed_regularization))

    # 'mu' is the only solver that handles other beta losses than 'frobenius'
    if solver != 'mu' and beta_loss not in (2, 'frobenius'):
        raise ValueError(
            'Invalid beta_loss parameter: solver %r does not handle beta_loss'
            ' = %r' % (solver, beta_loss))

    if solver == 'mu' and init == 'nndsvd':
        warnings.warn("The multiplicative update ('mu') solver cannot update "
                      "zeros present in the initialization, and so leads to "
                      "poorer results when used jointly with init='nndsvd'. "
                      "You may try init='nndsvda' or init='nndsvdar' instead.",
                      UserWarning)

    beta_loss = _beta_loss_to_float(beta_loss)
    return beta_loss


def _beta_loss_to_float(beta_loss):
    """Convert string beta_loss to float"""
    allowed_beta_loss = {'frobenius': 2,
                         'kullback-leibler': 1,
                         'itakura-saito': 0}
    if isinstance(beta_loss, str) and beta_loss in allowed_beta_loss:
        beta_loss = allowed_beta_loss[beta_loss]

    if not isinstance(beta_loss, numbers.Number):
        raise ValueError('Invalid beta_loss parameter: got %r instead '
                         'of one of %r, or a float.' %
                         (beta_loss, allowed_beta_loss.keys()))
    return beta_loss


def _initialize_nmf(X, n_components, init=None, eps=1e-6,
                    random_state=None):
    """Algorithms for NMF initialization.

    Computes an initial guess for the non-negative
    rank k matrix approximation for X: X = WH

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The data matrix to be decomposed.

    n_components : integer
        The number of components desired in the approximation.

    init :  None | 'random' | 'nndsvd' | 'nndsvda' | 'nndsvdar'
        Method used to initialize the procedure.
        Default: 'nndsvdar' if n_components < n_features, otherwise 'random'.
        Valid options:

        - 'random': non-negative random matrices, scaled with:
            sqrt(X.mean() / n_components)

        - 'nndsvd': Nonnegative Double Singular Value Decomposition (NNDSVD)
            initialization (better for sparseness)

        - 'nndsvda': NNDSVD with zeros filled with the average of X
            (better when sparsity is not desired)

        - 'nndsvdar': NNDSVD with zeros filled with small random values
            (generally faster, less accurate alternative to NNDSVDa
            for when sparsity is not desired)

        - 'custom': use custom matrices W and H

    eps : float
        Truncate all values less then this in output to zero.

    random_state : int seed, RandomState instance, or None (default)
        Random number generator seed control, used in 'nndsvdar' and
        'random' modes.

    Returns
    -------
    W : array-like, shape (n_samples, n_components)
        Initial guesses for solving X ~= WH

    H : array-like, shape (n_components, n_features)
        Initial guesses for solving X ~= WH

    References
    ----------
    C. Boutsidis, E. Gallopoulos: SVD based initialization: A head start for
    nonnegative matrix factorization - Pattern Recognition, 2008
    http://tinyurl.com/nndsvd
    """
    check_non_negative(X, "NMF initialization")
    n_samples, n_features = X.shape

    if init is None:
        if n_components < n_features:
            init = 'nndsvd'
        else:
            init = 'random'

    # Random initialization
    if init == 'random':
        avg = np.sqrt(X.mean() / n_components)
        rng = check_random_state(random_state)
        H = avg * rng.randn(n_components, n_features)
        W = avg * rng.randn(n_samples, n_components)
        # we do not write np.abs(H, out=H) to stay compatible with
        # numpy 1.5 and earlier where the 'out' keyword is not
        # supported as a kwarg on ufuncs
        np.abs(H, H)
        np.abs(W, W)
        return W, H

    # NNDSVD initialization
    U, S, V = randomized_svd(X, n_components, random_state=random_state)
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

    if init == "nndsvd":
        pass
    elif init == "nndsvda":
        avg = X.mean()
        W[W == 0] = avg
        H[H == 0] = avg
    elif init == "nndsvdar":
        rng = check_random_state(random_state)
        avg = X.mean()
        W[W == 0] = abs(avg * rng.randn(len(W[W == 0])) / 100)
        H[H == 0] = abs(avg * rng.randn(len(H[H == 0])) / 100)
    else:
        raise ValueError(
            'Invalid init parameter: got %r instead of one of %r' %
            (init, (None, 'random', 'nndsvd', 'nndsvda', 'nndsvdar')))

    return W, H


def _nls_subproblem(V, W, H, tol, max_iter, l1_reg=0., l2_reg=0.,
                    sigma=0.01, beta=0.1):
    """Non-negative least square solver

    Solves a non-negative least squares subproblem using the projected
    gradient descent algorithm.

    Parameters
    ----------
    V : array-like, shape (n_samples, n_features)
        Constant matrix.

    W : array-like, shape (n_samples, n_components)
        Constant matrix.

    H : array-like, shape (n_components, n_features)
        Initial guess for the solution.

    tol : float
        Tolerance of the stopping condition.

    max_iter : int
        Maximum number of iterations before timing out.

    l1_reg : double, default: 0.
        L1 regularization parameter.

    l2_reg : double, default: 0.
        L2 regularization parameter.

    sigma : float
        Constant used in the sufficient decrease condition checked by the line
        search.  Smaller values lead to a looser sufficient decrease condition,
        thus reducing the time taken by the line search, but potentially
        increasing the number of iterations of the projected gradient
        procedure. 0.01 is a commonly used value in the optimization
        literature.

    beta : float
        Factor by which the step size is decreased (resp. increased) until
        (resp. as long as) the sufficient decrease condition is satisfied.
        Larger values allow to find a better step size but lead to longer line
        search. 0.1 is a commonly used value in the optimization literature.

    Returns
    -------
    H : array-like, shape (n_components, n_features)
        Solution to the non-negative least squares problem.

    grad : array-like, shape (n_components, n_features)
        The gradient.

    n_iter : int
        The number of iterations done by the algorithm.

    References
    ----------
    C.-J. Lin. Projected gradient methods for non-negative matrix
    factorization. Neural Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/
    """
    WtV = safe_sparse_dot(W.T, V)
    WtW = fast_dot(W.T, W)

    # values justified in the paper (alpha is renamed gamma)
    gamma = 1
    for n_iter in range(1, max_iter + 1):
        grad = np.dot(WtW, H) - WtV
        if l1_reg > 0:
            grad += l1_reg
        if l2_reg > 0:
            grad += l2_reg * H

        # The following multiplication with a boolean array is more than twice
        # as fast as indexing into grad.
        if norm(grad * np.logical_or(grad < 0, H > 0)) < tol:
            break

        Hp = H

        for inner_iter in range(20):
            # Gradient step.
            Hn = H - gamma * grad
            # Projection step.
            Hn *= Hn > 0
            d = Hn - H
            gradd = np.dot(grad.ravel(), d.ravel())
            dQd = np.dot(np.dot(WtW, d).ravel(), d.ravel())
            suff_decr = (1 - sigma) * gradd + 0.5 * dQd < 0
            if inner_iter == 0:
                decr_gamma = not suff_decr

            if decr_gamma:
                if suff_decr:
                    H = Hn
                    break
                else:
                    gamma *= beta
            elif not suff_decr or (Hp == Hn).all():
                H = Hp
                break
            else:
                gamma /= beta
                Hp = Hn

    if n_iter == max_iter:
        warnings.warn("Iteration limit reached in nls subproblem.")

    return H, grad, n_iter


def _update_projected_gradient_w(X, W, H, tolW, nls_max_iter, l1_reg, l2_reg,
                                 sparseness, beta, eta):
    """Helper function for _fit_projected_gradient"""
    n_samples, n_features = X.shape
    n_components_ = H.shape[0]

    if sparseness is None:
        Wt, gradW, iterW = _nls_subproblem(X.T, H.T, W.T, tolW, nls_max_iter,
                                           l1_reg=l1_reg, l2_reg=l2_reg)
    elif sparseness == 'data':
        Wt, gradW, iterW = _nls_subproblem(
            safe_vstack([X.T, np.zeros((1, n_samples))]),
            safe_vstack([H.T, np.sqrt(beta) * np.ones((1,
                         n_components_))]),
            W.T, tolW, nls_max_iter, l1_reg=l1_reg, l2_reg=l2_reg)
    elif sparseness == 'components':
        Wt, gradW, iterW = _nls_subproblem(
            safe_vstack([X.T,
                         np.zeros((n_components_, n_samples))]),
            safe_vstack([H.T,
                         np.sqrt(eta) * np.eye(n_components_)]),
            W.T, tolW, nls_max_iter, l1_reg=l1_reg, l2_reg=l2_reg)

    return Wt.T, gradW.T, iterW


def _update_projected_gradient_h(X, W, H, tolH, nls_max_iter, l1_reg, l2_reg,
                                 sparseness, beta, eta):
    """Helper function for _fit_projected_gradient"""
    n_samples, n_features = X.shape
    n_components_ = W.shape[1]

    if sparseness is None:
        H, gradH, iterH = _nls_subproblem(X, W, H, tolH, nls_max_iter,
                                          l1_reg=l1_reg, l2_reg=l2_reg)
    elif sparseness == 'data':
        H, gradH, iterH = _nls_subproblem(
            safe_vstack([X, np.zeros((n_components_, n_features))]),
            safe_vstack([W,
                         np.sqrt(eta) * np.eye(n_components_)]),
            H, tolH, nls_max_iter, l1_reg=l1_reg, l2_reg=l2_reg)
    elif sparseness == 'components':
        H, gradH, iterH = _nls_subproblem(
            safe_vstack([X, np.zeros((1, n_features))]),
            safe_vstack([W, np.sqrt(beta) * np.ones((1, n_components_))]),
            H, tolH, nls_max_iter, l1_reg=l1_reg, l2_reg=l2_reg)

    return H, gradH, iterH


def _fit_projected_gradient(X, W, H, tol, max_iter, nls_max_iter, l1_reg_W,
                            l1_reg_H, l2_reg_W, l2_reg_H,
                            sparseness, beta, eta):
    """Compute Non-negative Matrix Factorization (NMF) with Projected Gradient

    References
    ----------
    C.-J. Lin. Projected gradient methods for non-negative matrix
    factorization. Neural Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    P. Hoyer. Non-negative Matrix Factorization with Sparseness Constraints.
    Journal of Machine Learning Research 2004.
    """
    gradW = (np.dot(W, np.dot(H, H.T)) -
             safe_sparse_dot(X, H.T, dense_output=True))
    gradH = (np.dot(np.dot(W.T, W), H) -
             safe_sparse_dot(W.T, X, dense_output=True))

    init_grad = squared_norm(gradW) + squared_norm(gradH.T)
    # max(0.001, tol) to force alternating minimizations of W and H
    tolW = max(0.001, tol) * np.sqrt(init_grad)
    tolH = tolW

    for n_iter in range(1, max_iter + 1):
        # stopping condition
        # as discussed in paper
        proj_grad_W = squared_norm(gradW * np.logical_or(gradW < 0, W > 0))
        proj_grad_H = squared_norm(gradH * np.logical_or(gradH < 0, H > 0))

        if (proj_grad_W + proj_grad_H) / init_grad < tol ** 2:
            break

        # update W
        W, gradW, iterW = _update_projected_gradient_w(X, W, H, tolW,
                                                       nls_max_iter,
                                                       l1_reg_W, l2_reg_W,
                                                       sparseness, beta, eta)
        if iterW == 1:
            tolW = 0.1 * tolW

        # update H
        H, gradH, iterH = _update_projected_gradient_h(X, W, H, tolH,
                                                       nls_max_iter,
                                                       l1_reg_H, l2_reg_H,
                                                       sparseness, beta, eta)
        if iterH == 1:
            tolH = 0.1 * tolH

    H[H == 0] = 0   # fix up negative zeros

    if n_iter == max_iter:
        W, _, _ = _update_projected_gradient_w(X, W, H, tol, nls_max_iter,
                                               l1_reg_W, l2_reg_W, sparseness,
                                               beta, eta)

    return W, H, n_iter


def _update_coordinate_descent(X, W, Ht, l1_reg, l2_reg, shuffle,
                               random_state):
    """Helper function for _fit_coordinate_descent

    Update W to minimize the objective function, iterating once over all
    coordinates. By symmetry, to update H, one can call
    _update_coordinate_descent(X.T, Ht, W, ...)

    """
    n_components = Ht.shape[1]

    HHt = fast_dot(Ht.T, Ht)
    XHt = safe_sparse_dot(X, Ht)

    # L2 regularization corresponds to increase of the diagonal of HHt
    if l2_reg != 0.:
        # adds l2_reg only on the diagonal
        HHt.flat[::n_components + 1] += l2_reg
    # L1 regularization corresponds to decrease of each element of XHt
    if l1_reg != 0.:
        XHt -= l1_reg

    if shuffle:
        permutation = random_state.permutation(n_components)
    else:
        permutation = np.arange(n_components)
    # The following seems to be required on 64-bit Windows w/ Python 3.5.
    permutation = np.asarray(permutation, dtype=np.intp)
    return _update_cdnmf_fast(W, HHt, XHt, permutation)


def _fit_coordinate_descent(X, W, H, tol=1e-4, max_iter=200, l1_reg_W=0,
                            l1_reg_H=0, l2_reg_W=0, l2_reg_H=0, update_H=True,
                            verbose=0, shuffle=False, random_state=None):
    """Compute Non-negative Matrix Factorization (NMF) with Coordinate Descent

    The objective function is minimized with an alternating minimization of W
    and H. Each minimization is done with a cyclic (up to a permutation of the
    features) Coordinate Descent.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        Constant matrix.

    W : array-like, shape (n_samples, n_components)
        Initial guess for the solution.

    H : array-like, shape (n_components, n_features)
        Initial guess for the solution.

    tol : float, default: 1e-4
        Tolerance of the stopping condition.

    max_iter : integer, default: 200
        Maximum number of iterations before timing out.

    l1_reg_W : double, default: 0.
        L1 regularization parameter for W.

    l1_reg_H : double, default: 0.
        L1 regularization parameter for H.

    l2_reg_W : double, default: 0.
        L2 regularization parameter for W.

    l2_reg_H : double, default: 0.
        L2 regularization parameter for H.

    update_H : boolean, default: True
        Set to True, both W and H will be estimated from initial guesses.
        Set to False, only W will be estimated.

    verbose : integer, default: 0
        The verbosity level.

    shuffle : boolean, default: False
        If true, randomize the order of coordinates in the CD solver.

    random_state : integer seed, RandomState instance, or None (default)
        Random number generator seed control.

    Returns
    -------
    W : array-like, shape (n_samples, n_components)
        Solution to the non-negative least squares problem.

    H : array-like, shape (n_components, n_features)
        Solution to the non-negative least squares problem.

    n_iter : int
        The number of iterations done by the algorithm.

    References
    ----------
    Cichocki, Andrzej, and P. H. A. N. Anh-Huy. "Fast local algorithms for
    large scale nonnegative matrix and tensor factorizations."
    IEICE transactions on fundamentals of electronics, communications and
    computer sciences 92.3: 708-721, 2009.
    """
    # so W and Ht are both in C order in memory
    Ht = check_array(H.T, order='C')
    X = check_array(X, accept_sparse='csr')

    rng = check_random_state(random_state)

    for n_iter in range(max_iter):
        violation = 0.

        # Update W
        violation += _update_coordinate_descent(X, W, Ht, l1_reg_W,
                                                l2_reg_W, shuffle, rng)
        # Update H
        if update_H:
            violation += _update_coordinate_descent(X.T, Ht, W, l1_reg_H,
                                                    l2_reg_H, shuffle, rng)

        if n_iter == 0:
            violation_init = violation

        if violation_init == 0:
            break

        if verbose:
            print("violation:", violation / violation_init)

        if violation / violation_init <= tol:
            if verbose:
                print("Converged at iteration", n_iter + 1)
            break

    return W, Ht.T, n_iter


def _multiplicative_update_w(X, W, H, beta_loss, l1_reg_W, l2_reg_W, gamma,
                             H_sum=None, HHt=None, XHt=None, update_H=True):
    """update W in Multiplicative Update NMF"""
    if beta_loss == 2:
        # Numerator
        if XHt is None:
            XHt = safe_sparse_dot(X, H.T)
        if update_H:
            # avoid a copy of XHt, which will be re-computed (update_H=True)
            numerator = XHt
        else:
            # preserve the XHt, which is not re-computed (update_H=False)
            numerator = XHt.copy()

        # Denominator
        if HHt is None:
            HHt = fast_dot(H, H.T)
        denominator = fast_dot(W, HHt)

    else:
        # Numerator
        # if X is sparse, compute WH only where X is non zero
        WH_safe_X = _special_sparse_dot(W, H, X)
        if sp.issparse(X):
            WH_safe_X_data = WH_safe_X.data
            X_data = X.data
        else:
            WH_safe_X_data = WH_safe_X
            X_data = X
            # copy used in the Denominator
            WH = WH_safe_X.copy()

        # to avoid taking a negative power of zero
        if beta_loss - 2. < 0:
            WH_safe_X_data[WH_safe_X_data == 0] = EPSILSON

        if beta_loss == 1:
            np.divide(X_data, WH_safe_X_data, out=WH_safe_X_data)
        else:
            WH_safe_X_data **= beta_loss - 2
            # element-wise multiplication
            WH_safe_X_data *= X_data

        # here numerator = dot(X * (dot(W, H) ** (beta_loss - 2), H.T)
        numerator = safe_sparse_dot(WH_safe_X, H.T)

        # Denominator
        if beta_loss == 1:
            if H_sum is None:
                H_sum = np.sum(H, axis=1)  # shape(n_components, )
            denominator = H_sum[np.newaxis, :]

        else:
            # computation of WHHt = dot(dot(W, H) ** beta_loss - 1, H.T)
            if sp.issparse(X):
                # memory efficient computation
                # (compute row by row, avoiding the dense matrix WH)
                WHHt = np.empty(W.shape)
                for i in range(X.shape[0]):
                    WHi = fast_dot(W[i, :], H)
                    WHi **= beta_loss - 1
                    WHHt[i, :] = fast_dot(WHi, H.T)
            else:
                WH **= beta_loss - 1
                WHHt = fast_dot(WH, H.T)
            denominator = WHHt

    # Add L1 and L2 regularization
    if l1_reg_W > 0:
        denominator += l1_reg_W
    if l2_reg_W > 0:
        denominator = denominator + l2_reg_W * W
    denominator[denominator == 0] = EPSILSON

    numerator /= denominator
    delta_W = numerator

    # gamma is in ]0, 1]
    if gamma != 1:
        delta_W **= gamma

    return delta_W, H_sum, HHt, XHt


def _multiplicative_update_h(X, W, H, beta_loss, l1_reg_H, l2_reg_H, gamma):
    """update H in Multiplicative Update NMF"""
    if beta_loss == 2:
        numerator = safe_sparse_dot(W.T, X)
        denominator = fast_dot(fast_dot(W.T, W), H)

    else:
        # Numerator
        WH_safe_X = _special_sparse_dot(W, H, X)
        if sp.issparse(X):
            WH_safe_X_data = WH_safe_X.data
            X_data = X.data
        else:
            WH_safe_X_data = WH_safe_X
            X_data = X
            # copy used in the Denominator
            WH = WH_safe_X.copy()

        # to avoid division by zero
        if beta_loss - 2. < 0:
            WH_safe_X_data[WH_safe_X_data == 0] = EPSILSON

        if beta_loss == 1:
            np.divide(X_data, WH_safe_X_data, out=WH_safe_X_data)
        else:
            WH_safe_X_data **= beta_loss - 2
            # element-wise multiplication
            WH_safe_X_data *= X_data

        # here numerator = dot(W.T, (dot(W, H) ** (beta_loss - 2)) * X)
        numerator = safe_sparse_dot(W.T, WH_safe_X)

        # Denominator
        if beta_loss == 1:
            W_sum = np.sum(W, axis=0)  # shape(n_components, )
            W_sum[W_sum == 0] = 1.
            denominator = W_sum[:, np.newaxis]

        # beta_loss not in (1, 2)
        else:
            # computation of WtWH = dot(W.T, dot(W, H) ** beta_loss - 1)
            if sp.issparse(X):
                # memory efficient computation
                # (compute column by column, avoiding the dense matrix WH)
                WtWH = np.empty(H.shape)
                for i in range(X.shape[1]):
                    WHi = fast_dot(W, H[:, i])
                    WHi **= beta_loss - 1
                    WtWH[:, i] = fast_dot(W.T, WHi)
            else:
                WH **= beta_loss - 1
                WtWH = fast_dot(W.T, WH)
            denominator = WtWH

    # Add L1 and L2 regularization
    if l1_reg_H > 0:
        denominator += l1_reg_H
    if l2_reg_H > 0:
        denominator = denominator + l2_reg_H * H
    denominator[denominator == 0] = EPSILSON

    numerator /= denominator
    delta_H = numerator

    # gamma is in ]0, 1]
    if gamma != 1:
        delta_H **= gamma

    return delta_H


def _fit_multiplicative_update(X, W, H, beta_loss='frobenius',
                               max_iter=200, tol=1e-4,
                               l1_reg_W=0, l1_reg_H=0, l2_reg_W=0, l2_reg_H=0,
                               update_H=True, verbose=0):
    """Compute Non-negative Matrix Factorization with Multiplicative Update

    The objective function is beta_divergence(X, WH) and is minimized with an
    alternating minimization of W and H. Each minimization is done with a
    Multiplicative Update.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        Constant input matrix.

    W : array-like, shape (n_samples, n_components)
        Initial guess for the solution.

    H : array-like, shape (n_components, n_features)
        Initial guess for the solution.

    beta_loss : float or string, default 'frobenius'
        String must be in {'frobenius', 'kullback-leibler', 'itakura-saito'}.
        Beta divergence to be minimized, measuring the distance between X
        and the dot product WH. Note that values different from 'frobenius'
        (or 2) and 'kullback-leibler' (or 1) lead to significantly slower
        fits.

    max_iter : integer, default: 200
        Number of iterations.

    tol : float, default: 1e-4
        Tolerance of the stopping condition.

    l1_reg_W : double, default: 0.
        L1 regularization parameter for W.

    l1_reg_H : double, default: 0.
        L1 regularization parameter for H.

    l2_reg_W : double, default: 0.
        L2 regularization parameter for W.

    l2_reg_H : double, default: 0.
        L2 regularization parameter for H.

    update_H : boolean, default: True
        Set to True, both W and H will be estimated from initial guesses.
        Set to False, only W will be estimated.

    verbose : integer, default: 0
        The verbosity level.

    Returns
    -------
    W : array, shape (n_samples, n_components)
        Solution to the non-negative least squares problem.

    H : array, shape (n_components, n_features)
        Solution to the non-negative least squares problem.

    n_iter : int
        The number of iterations done by the algorithm.

    References
    ----------
    Fevotte, C., & Idier, J. (2011). Algorithms for nonnegative matrix
    factorization with the beta-divergence. Neural Computation, 23(9).
    """
    start_time = time.time()

    beta_loss = _beta_loss_to_float(beta_loss)

    # gamma for Maximization-Minimization (MM) algorithm [Fevotte 2011]
    if beta_loss < 1:
        gamma = 1. / (2. - beta_loss)
    elif beta_loss > 2:
        gamma = 1. / (beta_loss - 1.)
    else:
        gamma = 1.

    # used for the convergence criterion
    error_at_init = previous_error = beta_divergence(X, W, H, beta_loss)

    H_sum, HHt, XHt = None, None, None
    for n_iter in range(1, max_iter + 1):
        # update W
        # H_sum, HHt and XHt are saved and reused if not update_H
        delta_W, H_sum, HHt, XHt = _multiplicative_update_w(
            X, W, H, beta_loss, l1_reg_W, l2_reg_W, gamma,
            H_sum, HHt, XHt, update_H)
        W *= delta_W

        # update H
        if update_H:
            delta_H = _multiplicative_update_h(X, W, H, beta_loss, l1_reg_H,
                                               l2_reg_H, gamma)
            H *= delta_H

            # These values will be recomputed since H changed
            H_sum, HHt, XHt = None, None, None

        # test convergence criterion every 100 iterations
        if tol > 0 and n_iter % 100 == 0:
            error = beta_divergence(X, W, H, beta_loss)

            if verbose:
                iter_time = time.time()
                print("Epoch %02d reached after %.3f seconds, error: %f" %
                      (n_iter, iter_time - start_time, error))

            if (previous_error - error) / error_at_init < tol:
                break
            previous_error = error

    if verbose and n_iter % 100 != 0:
        end_time = time.time()
        print("Epoch %02d reached after %.3f seconds." %
              (n_iter, end_time - start_time))

    return W, H, n_iter


def non_negative_factorization(X, W=None, H=None, n_components=None,
                               init='random', update_H=True, solver='cd',
                               beta_loss='frobenius', tol=1e-4,
                               max_iter=200, alpha=0., l1_ratio=0.,
                               regularization=None, random_state=None,
                               verbose=0, shuffle=False, nls_max_iter=2000,
                               sparseness=None, beta=1, eta=0.1):
    """Compute Non-negative Matrix Factorization (NMF)

    Find two non-negative matrices (W, H) whose product approximates the non-
    negative matrix X. This factorization can be used for example for
    dimensionality reduction, source separation or topic extraction.

    The objective function is::

        0.5 * ||X - WH||_Fro^2
        + alpha * l1_ratio * ||vec(W)||_1
        + alpha * l1_ratio * ||vec(H)||_1
        + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2
        + 0.5 * alpha * (1 - l1_ratio) * ||H||_Fro^2

    Where::

        ||A||_Fro^2 = \sum_{i,j} A_{ij}^2 (Frobenius norm)
        ||vec(A)||_1 = \sum_{i,j} abs(A_{ij}) (Elementwise L1 norm)

    For multiplicative-update ('mu') solver, the Frobenius norm
    (0.5 * ||X - WH||_Fro^2) can be changed into another beta-divergence loss,
    by changing the beta_loss parameter.

    The objective function is minimized with an alternating minimization of W
    and H. If H is given and update_H=False, it solves for W only.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        Constant matrix.

    W : array-like, shape (n_samples, n_components)
        If init='custom', it is used as initial guess for the solution.

    H : array-like, shape (n_components, n_features)
        If init='custom', it is used as initial guess for the solution.
        If update_H=False, it is used as a constant, to solve for W only.

    n_components : integer
        Number of components, if n_components is not set all features
        are kept.

    init :  None | 'random' | 'nndsvd' | 'nndsvda' | 'nndsvdar' | 'custom'
        Method used to initialize the procedure.
        Default: 'nndsvd' if n_components < n_features, otherwise random.
        Valid options:

        - 'random': non-negative random matrices, scaled with:
            sqrt(X.mean() / n_components)

        - 'nndsvd': Nonnegative Double Singular Value Decomposition (NNDSVD)
            initialization (better for sparseness)

        - 'nndsvda': NNDSVD with zeros filled with the average of X
            (better when sparsity is not desired)

        - 'nndsvdar': NNDSVD with zeros filled with small random values
            (generally faster, less accurate alternative to NNDSVDa
            for when sparsity is not desired)

        - 'custom': use custom matrices W and H

    update_H : boolean, default: True
        Set to True, both W and H will be estimated from initial guesses.
        Set to False, only W will be estimated.

    solver : 'pg' | 'cd' | 'mu'
        Numerical solver to use:
        'pg' is a Projected Gradient solver (deprecated).
        'cd' is a Coordinate Descent solver (recommended).
        'mu' is a Multiplicative Update solver.

        .. versionadded:: 0.17
           Coordinate Descent solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver.

        .. versionadded:: 0.19
           Multiplicative Update solver.

    beta_loss : float or string, default 'frobenius'
        String must be in {'frobenius', 'kullback-leibler', 'itakura-saito'}.
        Beta divergence to be minimized, measuring the distance between X
        and the dot product WH. Note that values different from 'frobenius'
        (or 2) and 'kullback-leibler' (or 1) lead to significantly slower
        fits. Used only in 'mu' solver.

    tol : float, default: 1e-4
        Tolerance of the stopping condition.

    max_iter : integer, default: 200
        Maximum number of iterations before timing out.
        'mu' solver always does the maximum number of iteration.

    alpha : double, default: 0.
        Constant that multiplies the regularization terms.

    l1_ratio : double, default: 0.
        The regularization mixing parameter, with 0 <= l1_ratio <= 1.
        For l1_ratio = 0 the penalty is an elementwise L2 penalty
        (aka Frobenius Norm).
        For l1_ratio = 1 it is an elementwise L1 penalty.
        For 0 < l1_ratio < 1, the penalty is a combination of L1 and L2.

    regularization : 'both' | 'components' | 'transformation' | None
        Select whether the regularization affects the components (H), the
        transformation (W), both or none of them.

    random_state : integer seed, RandomState instance, or None (default)
        Random number generator seed control.

    verbose : integer, default: 0
        The verbosity level.

    shuffle : boolean, default: False
        If true, randomize the order of coordinates in the CD solver.

    nls_max_iter : integer, default: 2000
        Number of iterations in NLS subproblem.
        Used only in the deprecated 'pg' solver.

    sparseness : 'data' | 'components' | None, default: None
        Where to enforce sparsity in the model.
        Used only in the deprecated 'pg' solver.

    beta : double, default: 1
        Degree of sparseness, if sparseness is not None. Larger values mean
        more sparseness. Used only in the deprecated 'pg' solver.

    eta : double, default: 0.1
        Degree of correctness to maintain, if sparsity is not None. Smaller
        values mean larger error. Used only in the deprecated 'pg' solver.

    Returns
    -------
    W : array-like, shape (n_samples, n_components)
        Solution to the non-negative least squares problem.

    H : array-like, shape (n_components, n_features)
        Solution to the non-negative least squares problem.

    n_iter : int
        Actual number of iterations.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[1,1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from sklearn.decomposition import non_negative_factorization
    >>> W, H, n_iter = non_negative_factorization(X, n_components=2, \
        init='random', random_state=0)

    References
    ----------
    C.-J. Lin. Projected gradient methods for non-negative matrix
    factorization. Neural Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    Cichocki, Andrzej, and P. H. A. N. Anh-Huy. "Fast local algorithms for
    large scale nonnegative matrix and tensor factorizations."
    IEICE transactions on fundamentals of electronics, communications and
    computer sciences 92.3: 708-721, 2009.

    Fevotte, C., & Idier, J. (2011). Algorithms for nonnegative matrix
    factorization with the beta-divergence. Neural Computation, 23(9).
    """

    X = check_array(X, accept_sparse=('csr', 'csc'))
    check_non_negative(X, "NMF (input X)")
    beta_loss = _check_string_param(sparseness, solver, regularization,
                                    beta_loss, init)

    n_samples, n_features = X.shape
    if n_components is None:
        n_components = n_features

    if not isinstance(n_components, INTEGER_TYPES) or n_components <= 0:
        raise ValueError("Number of components must be a positive integer;"
                         " got (n_components=%r)" % n_components)
    if not isinstance(max_iter, INTEGER_TYPES) or max_iter < 0:
        raise ValueError("Maximum number of iterations must be a positive "
                         "integer; got (max_iter=%r)" % max_iter)
    if not isinstance(tol, numbers.Number) or tol < 0:
        raise ValueError("Tolerance for stopping criteria must be "
                         "positive; got (tol=%r)" % tol)

    # check W and H, or initialize them
    if init == 'custom' and update_H:
        _check_init(H, (n_components, n_features), "NMF (input H)")
        _check_init(W, (n_samples, n_components), "NMF (input W)")
    elif not update_H:
        _check_init(H, (n_components, n_features), "NMF (input H)")
        # 'mu' solver should not be initialized by zeros
        if solver == 'mu':
            avg = np.sqrt(X.mean() / n_components)
            W = avg * np.ones((n_samples, n_components))
        else:
            W = np.zeros((n_samples, n_components))
    else:
        W, H = _initialize_nmf(X, n_components, init=init,
                               random_state=random_state)

    l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H = _compute_regularization(
        alpha, l1_ratio, regularization)

    if solver == 'pg':
        warnings.warn("'pg' solver will be removed in release 0.19."
                      " Use 'cd' solver instead.", DeprecationWarning)
        if update_H:  # fit_transform
            W, H, n_iter = _fit_projected_gradient(X, W, H, tol,
                                                   max_iter,
                                                   nls_max_iter,
                                                   l1_reg_W, l1_reg_H,
                                                   l2_reg_W, l2_reg_H,
                                                   sparseness,
                                                   beta, eta)
        else:  # transform
            W, H, n_iter = _update_projected_gradient_w(X, W, H,
                                                        tol, nls_max_iter,
                                                        l1_reg_W, l2_reg_W,
                                                        sparseness, beta,
                                                        eta)
    elif solver == 'cd':
        W, H, n_iter = _fit_coordinate_descent(X, W, H, tol, max_iter,
                                               l1_reg_W, l1_reg_H,
                                               l2_reg_W, l2_reg_H,
                                               update_H=update_H,
                                               verbose=verbose,
                                               shuffle=shuffle,
                                               random_state=random_state)
    elif solver == 'mu':
        W, H, n_iter = _fit_multiplicative_update(X, W, H, beta_loss, max_iter,
                                                  tol, l1_reg_W, l1_reg_H,
                                                  l2_reg_W, l2_reg_H, update_H,
                                                  verbose)

    else:
        raise ValueError("Invalid solver parameter '%s'." % solver)

    if n_iter == max_iter:
        warnings.warn("Maximum number of iteration %d reached. Increase it to"
                      " improve convergence." % max_iter, ConvergenceWarning)

    return W, H, n_iter


class NMF(BaseEstimator, TransformerMixin):
    """Non-Negative Matrix Factorization (NMF)

    Find two non-negative matrices (W, H) whose product approximates the non-
    negative matrix X. This factorization can be used for example for
    dimensionality reduction, source separation or topic extraction.

    The objective function is::

        0.5 * ||X - WH||_Fro^2
        + alpha * l1_ratio * ||vec(W)||_1
        + alpha * l1_ratio * ||vec(H)||_1
        + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2
        + 0.5 * alpha * (1 - l1_ratio) * ||H||_Fro^2

    Where::

        ||A||_Fro^2 = \sum_{i,j} A_{ij}^2 (Frobenius norm)
        ||vec(A)||_1 = \sum_{i,j} abs(A_{ij}) (Elementwise L1 norm)

    For multiplicative-update ('mu') solver, the Frobenius norm
    (0.5 * ||X - WH||_Fro^2) can be changed into another beta-divergence loss,
    by changing the beta_loss parameter.

    The objective function is minimized with an alternating minimization of W
    and H.

    Read more in the :ref:`User Guide <NMF>`.

    Parameters
    ----------
    n_components : int or None
        Number of components, if n_components is not set all features
        are kept.

    init :  'random' | 'nndsvd' |  'nndsvda' | 'nndsvdar' | 'custom'
        Method used to initialize the procedure.
        Default: 'nndsvdar' if n_components < n_features, otherwise random.
        Valid options:

        - 'random': non-negative random matrices, scaled with:
            sqrt(X.mean() / n_components)

        - 'nndsvd': Nonnegative Double Singular Value Decomposition (NNDSVD)
            initialization (better for sparseness)

        - 'nndsvda': NNDSVD with zeros filled with the average of X
            (better when sparsity is not desired)

        - 'nndsvdar': NNDSVD with zeros filled with small random values
            (generally faster, less accurate alternative to NNDSVDa
            for when sparsity is not desired)

        - 'custom': use custom matrices W and H

    solver : 'pg' | 'cd' | 'mu'
        Numerical solver to use:
        'pg' is a Projected Gradient solver (deprecated).
        'cd' is a Coordinate Descent solver (recommended).
        'mu' is a Multiplicative Update solver.

        .. versionadded:: 0.17
           Coordinate Descent solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver.

        .. versionadded:: 0.19
           Multiplicative Update solver.


    beta_loss : float or string, default 'frobenius'
        String must be in {'frobenius', 'kullback-leibler', 'itakura-saito'}.
        Beta divergence to be minimized, measuring the distance between X
        and the dot product WH. Note that values different from 'frobenius'
        (or 2) and 'kullback-leibler' (or 1) lead to significantly slower
        fits. Used only in 'mu' solver.

    tol : float, default: 1e-4
        Tolerance of the stopping condition.

    max_iter : integer, default: 200
        Maximum number of iterations before timing out.
        'mu' solver always does the maximum number of iteration.

    random_state : integer seed, RandomState instance, or None (default)
        Random number generator seed control.

    alpha : double, default: 0.
        Constant that multiplies the regularization terms. Set it to zero to
        have no regularization.

        .. versionadded:: 0.17
           *alpha* used in the Coordinate Descent solver.

    l1_ratio : double, default: 0.
        The regularization mixing parameter, with 0 <= l1_ratio <= 1.
        For l1_ratio = 0 the penalty is an elementwise L2 penalty
        (aka Frobenius Norm).
        For l1_ratio = 1 it is an elementwise L1 penalty.
        For 0 < l1_ratio < 1, the penalty is a combination of L1 and L2.

        .. versionadded:: 0.17
           Regularization parameter *l1_ratio* used in the Coordinate Descent
           solver.

    shuffle : boolean, default: False
        If true, randomize the order of coordinates in the CD solver.

        .. versionadded:: 0.17
           *shuffle* parameter used in the Coordinate Descent solver.

    nls_max_iter : integer, default: 2000
        Number of iterations in NLS subproblem.
        Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    sparseness : 'data' | 'components' | None, default: None
        Where to enforce sparsity in the model.
        Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    beta : double, default: 1
        Degree of sparseness, if sparseness is not None. Larger values mean
        more sparseness. Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    eta : double, default: 0.1
        Degree of correctness to maintain, if sparsity is not None. Smaller
        values mean larger error. Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    Attributes
    ----------
    components_ : array, [n_components, n_features]
        Non-negative components of the data.

    reconstruction_err_ : number
        Frobenius norm of the matrix difference between
        the training data and the reconstructed data from
        the fit produced by the model. ``|| X - WH ||_2``

    n_iter_ : int
        Actual number of iterations.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[1, 1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from sklearn.decomposition import NMF
    >>> model = NMF(n_components=2, init='random', random_state=0)
    >>> W = model.fit_transform(X)
    >>> H = model.components_

    References
    ----------
    C.-J. Lin. Projected gradient methods for non-negative matrix
    factorization. Neural Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    Cichocki, Andrzej, and P. H. A. N. Anh-Huy. "Fast local algorithms for
    large scale nonnegative matrix and tensor factorizations."
    IEICE transactions on fundamentals of electronics, communications and
    computer sciences 92.3: 708-721, 2009.

    Fevotte, C., & Idier, J. (2011). Algorithms for nonnegative matrix
    factorization with the beta-divergence. Neural Computation, 23(9).
    """
    def __init__(self, n_components=None, init=None, solver='cd',
                 beta_loss='frobenius', tol=1e-4, max_iter=200,
                 random_state=None, alpha=0., l1_ratio=0., verbose=0,
                 shuffle=False, nls_max_iter=2000, sparseness=None,
                 beta=1, eta=0.1):
        self.n_components = n_components
        self.init = init
        self.solver = solver
        self.beta_loss = beta_loss
        self.tol = tol
        self.max_iter = max_iter
        self.random_state = random_state
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.verbose = verbose
        self.shuffle = shuffle

        if sparseness is not None:
            warnings.warn("Controlling regularization through the sparseness,"
                          " beta and eta arguments is only available"
                          " for 'pg' solver, which will be removed"
                          " in release 0.19. Use another solver with L1 or L2"
                          " regularization instead.", DeprecationWarning)
        self.nls_max_iter = nls_max_iter
        self.sparseness = sparseness
        self.beta = beta
        self.eta = eta

    def fit_transform(self, X, y=None, W=None, H=None):
        """Learn a NMF model for the data X and returns the transformed data.

        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Data matrix to be decomposed

        W : array-like, shape (n_samples, n_components)
            If init='custom', it is used as initial guess for the solution.

        H : array-like, shape (n_components, n_features)
            If init='custom', it is used as initial guess for the solution.

        Attributes
        ----------
        components_ : array-like, shape (n_components, n_features)
            Factorization matrix, sometimes called 'dictionary'.

        n_iter_ : int
            Actual number of iterations for the transform.

        Returns
        -------
        W: array, shape (n_samples, n_components)
            Transformed data.
        """
        X = check_array(X, accept_sparse=('csr', 'csc'))

        W, H, n_iter_ = non_negative_factorization(
            X=X, W=W, H=H, n_components=self.n_components, init=self.init,
            update_H=True, solver=self.solver, beta_loss=self.beta_loss,
            tol=self.tol, max_iter=self.max_iter, alpha=self.alpha,
            l1_ratio=self.l1_ratio, regularization='both',
            random_state=self.random_state, verbose=self.verbose,
            shuffle=self.shuffle,
            nls_max_iter=self.nls_max_iter, sparseness=self.sparseness,
            beta=self.beta, eta=self.eta)

        if self.solver == 'pg':
            self.comp_sparseness_ = _sparseness(H.ravel())
            self.data_sparseness_ = _sparseness(W.ravel())

        if self.solver == 'mu':
            self.reconstruction_err_ = beta_divergence(X, W, H,
                                                       self.beta_loss)
        else:
            self.reconstruction_err_ = beta_divergence(X, W, H, 2.)

        self.n_components_ = H.shape[0]
        self.components_ = H
        self.n_iter_ = n_iter_

        return W

    def fit(self, X, y=None, **params):
        """Learn a NMF model for the data X.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Data matrix to be decomposed

        Attributes
        ----------
        components_ : array-like, shape (n_components, n_features)
            Factorization matrix, sometimes called 'dictionary'.

        n_iter_ : int
            Actual number of iterations for the transform.

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
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Data matrix to be transformed by the model

        Attributes
        ----------
        n_iter_ : int
            Actual number of iterations for the transform.

        Returns
        -------
        W: array, shape (n_samples, n_components)
            Transformed data
        """
        check_is_fitted(self, 'n_components_')

        W, _, n_iter_ = non_negative_factorization(
            X=X, W=None, H=self.components_, n_components=self.n_components_,
            init=self.init, update_H=False, solver=self.solver,
            beta_loss=self.beta_loss, tol=self.tol, max_iter=self.max_iter,
            alpha=self.alpha, l1_ratio=self.l1_ratio, regularization='both',
            random_state=self.random_state, verbose=self.verbose,
            shuffle=self.shuffle, nls_max_iter=self.nls_max_iter,
            sparseness=self.sparseness, beta=self.beta, eta=self.eta)

        self.n_iter_ = n_iter_
        return W

    def inverse_transform(self, W):
        """Transform data back to its original space.

        Parameters
        ----------
        W: {array-like, sparse matrix}, shape (n_samples, n_components)
            Transformed data matrix

        Returns
        -------
        X: {array-like, sparse matrix}, shape (n_samples, n_features)
            Data matrix of original shape

        .. versionadded:: 0.18
        """
        check_is_fitted(self, 'n_components_')
        return np.dot(W, self.components_)


@deprecated("It will be removed in release 0.19. Use NMF instead."
            "'pg' solver is still available until release 0.19.")
class ProjectedGradientNMF(NMF):
    """Non-Negative Matrix Factorization (NMF)

    Find two non-negative matrices (W, H) whose product approximates the non-
    negative matrix X. This factorization can be used for example for
    dimensionality reduction, source separation or topic extraction.

    The objective function is::

        0.5 * ||X - WH||_Fro^2
        + alpha * l1_ratio * ||vec(W)||_1
        + alpha * l1_ratio * ||vec(H)||_1
        + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2
        + 0.5 * alpha * (1 - l1_ratio) * ||H||_Fro^2

    Where::

        ||A||_Fro^2 = \sum_{i,j} A_{ij}^2 (Frobenius norm)
        ||vec(A)||_1 = \sum_{i,j} abs(A_{ij}) (Elementwise L1 norm)

    The objective function is minimized with an alternating minimization of W
    and H.

    Read more in the :ref:`User Guide <NMF>`.

    Parameters
    ----------
    n_components : int or None
        Number of components, if n_components is not set all features
        are kept.

    init :  'random' | 'nndsvd' |  'nndsvda' | 'nndsvdar' | 'custom'
        Method used to initialize the procedure.
        Default: 'nndsvdar' if n_components < n_features, otherwise random.
        Valid options:

        - 'random': non-negative random matrices, scaled with:
            sqrt(X.mean() / n_components)

        - 'nndsvd': Nonnegative Double Singular Value Decomposition (NNDSVD)
            initialization (better for sparseness)

        - 'nndsvda': NNDSVD with zeros filled with the average of X
            (better when sparsity is not desired)

        - 'nndsvdar': NNDSVD with zeros filled with small random values
            (generally faster, less accurate alternative to NNDSVDa
            for when sparsity is not desired)

        - 'custom': use custom matrices W and H

    solver : 'pg' | 'cd' | 'mu'
        Numerical solver to use:
        'pg' is a Projected Gradient solver (deprecated).
        'cd' is a Coordinate Descent solver (recommended).
        'mu' is a Multiplicative Update solver.

        .. versionadded:: 0.17
           Coordinate Descent solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver.

        .. versionadded:: 0.19
           Multiplicative Update solver.

    tol : double, default: 1e-4
        Tolerance value used in stopping conditions.

    max_iter : integer, default: 200
        Number of iterations to compute.

    random_state : integer seed, RandomState instance, or None (default)
        Random number generator seed control.

    alpha : double, default: 0.
        Constant that multiplies the regularization terms. Set it to zero to
        have no regularization.

        .. versionadded:: 0.17
           *alpha* used in the Coordinate Descent solver.

    l1_ratio : double, default: 0.
        The regularization mixing parameter, with 0 <= l1_ratio <= 1.
        For l1_ratio = 0 the penalty is an elementwise L2 penalty
        (aka Frobenius Norm).
        For l1_ratio = 1 it is an elementwise L1 penalty.
        For 0 < l1_ratio < 1, the penalty is a combination of L1 and L2.

        .. versionadded:: 0.17
           Regularization parameter *l1_ratio* used in the Coordinate Descent
           solver.

    shuffle : boolean, default: False
        If true, randomize the order of coordinates in the CD solver.

        .. versionadded:: 0.17
           *shuffle* parameter used in the Coordinate Descent solver.

    nls_max_iter : integer, default: 2000
        Number of iterations in NLS subproblem.
        Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    sparseness : 'data' | 'components' | None, default: None
        Where to enforce sparsity in the model.
        Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    beta : double, default: 1
        Degree of sparseness, if sparseness is not None. Larger values mean
        more sparseness. Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    eta : double, default: 0.1
        Degree of correctness to maintain, if sparsity is not None. Smaller
        values mean larger error. Used only in the deprecated 'pg' solver.

        .. versionchanged:: 0.17
           Deprecated Projected Gradient solver. Use Coordinate Descent solver
           instead.

    Attributes
    ----------
    components_ : array, [n_components, n_features]
        Non-negative components of the data.

    reconstruction_err_ : number
        Frobenius norm of the matrix difference between
        the training data and the reconstructed data from
        the fit produced by the model. ``|| X - WH ||_2``

    n_iter_ : int
        Actual number of iterations.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[1, 1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from sklearn.decomposition import NMF
    >>> model = NMF(n_components=2, init='random', random_state=0)
    >>> W = model.fit_transform(X)
    >>> H = model.components_

    References
    ----------
    C.-J. Lin. Projected gradient methods for non-negative matrix
    factorization. Neural Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    Cichocki, Andrzej, and P. H. A. N. Anh-Huy. "Fast local algorithms for
    large scale nonnegative matrix and tensor factorizations."
    IEICE transactions on fundamentals of electronics, communications and
    computer sciences 92.3: 708-721, 2009.
    """

    def __init__(self, n_components=None, solver='pg', beta_loss='frobenius',
                 init=None, tol=1e-4, max_iter=200, random_state=None,
                 alpha=0., l1_ratio=0., verbose=0, nls_max_iter=2000,
                 sparseness=None, beta=1, eta=0.1):
        super(ProjectedGradientNMF, self).__init__(
            n_components=n_components, init=init, solver='pg',
            beta_loss=beta_loss, tol=tol, max_iter=max_iter,
            random_state=random_state, alpha=alpha, l1_ratio=l1_ratio,
            verbose=verbose, nls_max_iter=nls_max_iter, sparseness=sparseness,
            beta=beta, eta=eta)
