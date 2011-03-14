""" Non-negative matrix factorization
"""
# Author: Chih-Jen Lin, National Taiwan University
# Python/numpy translation: Anthony Di Franco
# scikit.learn integration: Vlad Niculae
# License: BSD


import warnings

import numpy as np

from .base import BaseEstimator
from .utils.extmath import fast_svd
from numpy.linalg import norm

_pos_ = lambda x: (x >= 0) * x
_neg_ = lambda x: (x < 0) * (-x)

def _initialize_nmf_(X, n_comp):
    """
    Computes a good initial guess for the non-negative
    rank k matrix approximation for X: X = WH

    Parameters
    ----------

    X:
        The data matrix to be decomposed.

    n_comp:
        The number of components desired in the
        approximation.

    Returns
    -------

    (W, H):
        Initial guesses for solving X ~= WH such that
        the number of columns in W is n_comp.

    Remarks
    -------

    This implements the algorithm described in
    C. Boutsidis, E. Gallopoulos: SVD based
    initialization: A head start for nonnegative
    matrix factorization - Pattern Recognition, 2008

    http://www.cs.rpi.edu/~boutsc/files/nndsvd.pdf
    """
    if (X < 0).any():
        raise ValueError("Negative values in data passed to initialization")

    U, S, V = fast_svd(X, n_comp)
    W, H = np.zeros(U.shape), np.zeros(V.shape)

    ## The leading singular triplet is non-negative
    ## so it can be used as is for initialization.
    # !!!!! apparently not!
    # W[:, 0] = np.sqrt(S[0]) * U[:, 0]
    # H[0, :] = np.sqrt(S[0]) * V[0, :]

    for j in xrange(0, n_comp):
        x, y = U[:, j], V[j, :]

        # extract positive and negative parts of column vectors
        x_p, y_p = _pos_(x), _pos_(y)
        x_n, y_n = _neg_(x), _neg_(y)

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

    return W, H


def _nls_subproblem_(V, W, Hinit, tolerance, max_iter):
    """
    Solves a non-negative least squares subproblem using the
    projected gradient descent algorithm.
    min || WH - V ||_2

    Parameters
    ----------
    V, W:
        Constant matrices

    Hinit:
        Initial guess for the solution

    tolerance:
        Tolerance of the stopping condition.

    max_iter:
        Maximum number of iterations before
        timing out.

    Returns
    -------
    H:
        Solution to the non-negative least squares problem

    grad:
        The gradient.

    iter:
        The number of iterations done by the algorithm.

    """
    if (Hinit < 0).any():
        raise ValueError("Negative values in Hinit passed to NLS solver.")

    H = Hinit
    WtV = np.dot(W.T, V)
    WtW = np.dot(W.T, W)

    # values justified in the paper
    alpha = 1
    beta = 0.1
    for iter in xrange(1, max_iter + 1):
        grad = np.dot(WtW, H) - WtV
        proj_gradient = norm(grad[np.logical_or(grad < 0, H > 0)])
        if proj_gradient < tolerance:
            break

        for inner_iter in xrange(1, 20):
            Hn = H - alpha * grad
            # Hn = np.where(Hn > 0, Hn, 0)
            Hn = _pos_(Hn)
            d = Hn - H
            gradd = np.sum(grad * d)
            dQd = np.sum(np.dot(WtW, d) * d)
            # magic numbers whoa
            suff_decr = 0.99 * gradd + 0.5 * dQd < 0
            if inner_iter == 1:
                decr_alpha = not suff_decr
                Hp = H

            if decr_alpha:
                if suff_decr:
                    H = Hn
                    break
                else:
                    alpha = alpha * beta
            else:
                if not suff_decr or (Hp == Hn).all():
                    H = Hp
                    break
                else:
                    alpha = alpha / beta
                    Hp = Hn

    if iter == max_iter:
        warnings.warn("Iteration limit reached in nls subproblem.")

    return H, grad, iter


class NMF(BaseEstimator):
    """
    Non-Negative matrix factorization (NMF, NNMF)

    Parameters
    ----------
    X: array, [n_samples, n_features]
        Data the model will be fit to.

    n_comp: int or None
        Number of components
        if n_comp is not set all components are kept

    initial: "fast_svd", int or RandomState
        Method used to initialize the procedure.
        Default: fast_svd
        Valid options:
            "fast_svd" for initialization based on fast SVD,
            "cro" for CRO-based initialization,
            int seed or RandomState for non-negative random matrices

    tolerance: double
        Tolerance value used in stopping conditions.
        Default: 0.001

    max_iter: int
        Number of iterations to compute.
        Default: 100

    nls_max_iter: int
        Number of iterations in NLS subproblem
        Default: 2000

    Attributes
    ----------
    components_: array, [n_features, n_comp]
        Non-negative components of the data
    data_: array, [n_comp, n_samples]
        Projection of the data used to fit the model onto
        the non-negative components.
    reconstruction_err_: number
        Frobenius norm of the matrix difference between the
        training data and the reconstructed data from the
        fit produced by the model. || X - WH ||_2

    Examples
    --------

    >>> import numpy as np
    >>> X = np.array([[1,1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from scikits.learn.nmf import NMF
    >>> from scipy.optimize import nnls
    >>> model = NMF(n_comp=2, initial=0)
    >>> model.fit(X) #doctest: +ELLIPSIS
    NMF(nls_max_iter=1000, n_comp=2,
      initial=<mtrand.RandomState object at 0x...>, tolerance=0.001,
      max_iter=50)
    >>> model.components_
    array([[ 0.76551648,  0.11408118],
           [ 0.4167371 ,  0.36376783]])
    >>> model.reconstruction_err_
    0.11760278071752467
    >>> np.linalg.norm(model.data_[0, :] - model.transform(X[0, :]))
    0.023723012161826595

    Notes
    -----
    This implements C.-J. Lin. Projected gradient methods
    for non-negative matrix factorization. Neural
    Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    """

    def __init__(self, n_comp=None, initial="fast_svd",
                tolerance=0.001, max_iter=100, nls_max_iter=2000):
        self.n_comp = n_comp
        self.initial = initial
        self.tolerance = tolerance
        self.max_iter = max_iter
        self.nls_max_iter = nls_max_iter

    def fit(self, X):
        """ Fit the model to the data
        """
        X = np.atleast_2d(X)
        if (X < 0).any():
            raise ValueError("Negative data passed to NMF.fit.")

        n_features, n_samples = X.shape

        if not self.n_comp:
            self.n_comp = n_features

        if self.initial == None:
            self.initial = np.random.RandomState()
        elif isinstance(self.initial, int):
            self.initial = np.random.RandomState(self.initial)

        if isinstance(self.initial, np.random.RandomState):
            W = np.abs(self.initial.randn(n_features, self.n_comp))
            H = np.abs(self.initial.randn(self.n_comp, n_samples))
        elif self.initial == "fast_svd":
            W, H = _initialize_nmf_(X, self.n_comp)
        elif self.initial == "cro":
            from cro import CRO
            m = CRO(self.n_comp)
            m.fit(X.T)
            W, H = np.abs(m.components_.T), np.abs(m.data_.T)
        else:
            raise ValueError("Invalid value for initial parameter.")

        gradW = np.dot(W, np.dot(H, H.T)) - np.dot(X, H.T)
        gradH = np.dot(np.dot(W.T, W), H) - np.dot(W.T, X)
        init_grad = norm(np.r_[gradW, gradH.T])
        tolW = max(0.001, self.tolerance) * init_grad  # why max?
        tolH = tolW

        for iter in xrange(1, self.max_iter + 1):
            # stopping condition
            # as discussed in paper
            proj_norm = norm(np.r_[gradW[np.logical_or(gradW < 0, W > 0)],
                                   gradH[np.logical_or(gradH < 0, H > 0)]])
            if proj_norm < self.tolerance * init_grad:
                break

            # update W
            W, gradW, iterW = _nls_subproblem_(X.T, H.T, W.T, tolW,
                                               self.nls_max_iter)
            W = W.T
            gradW = gradW.T
            if iterW == 1:
                tolW = 0.1 * tolW

            # update H
            H, gradH, iterH = _nls_subproblem_(X, W, H, tolH,
                                               self.nls_max_iter)
            if iterH == 1:
                tolH = 0.1 * tolH

            self.components_ = H
            self.data_ = W
            self.reconstruction_err_ = norm(X - np.dot(W, H))
        if iter == self.max_iter:
            warnings.warn("Iteration limit reached during fit")
        return self

    def transform(self, X):
        """ Transform the data X according to the model
        """
        from scipy.optimize import nnls
        X = np.atleast_2d(X)
        H = np.zeros((X.shape[0], self.n_comp))
        for j in xrange(0, X.shape[0]):
            H[j, :], _ = nnls(self.components_.T, X[j, :])
        return H
