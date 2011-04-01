""" Non-negative matrix factorization
"""
# Author: Chih-Jen Lin, National Taiwan University
# Python/numpy translation: Anthony Di Franco
# scikit.learn integration: Vlad Niculae
# License: BSD


from __future__ import division
import warnings

import numpy as np
from numpy import r_, zeros, ones, eye, sqrt
from .base import BaseEstimator
from .utils.extmath import fast_svd
from numpy.linalg import norm

_pos_ = lambda x: (x >= 0) * x
_neg_ = lambda x: (x < 0) * (-x)


class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def _sparseness_(x):
    """
    Hoyer's measure of sparsity for a vector
    """
    n = len(x)
    return (sqrt(n) - norm(x, 1) / norm(x, 2)) / (sqrt(n) - 1)


def dict_argmax(dictionary):
    val = lambda x: x[1]
    return max(dictionary.items(), key=val)[0]


def _initialize_nmf_(X, n_comp, variant=None, eps=1e-6, seed=None):
    """
    NNDSVD algorithm for NMF initialization.
    Computes a good initial guess for the non-negative
    rank k matrix approximation for X: X = WH

    Parameters
    ----------

    X:
        The data matrix to be decomposed.

    n_comp:
        The number of components desired in the
        approximation.

    variant:
        The variant of the NNDSVD algorithm.
        Accepts None, "a", "ar"
        None: leaves the zero entries as zero
        "a": Fills the zero entries with the average of X
        "ar": Fills the zero entries with standard normal random variates.

    eps:
        Truncate all values less then this in output to zero.

    seed:
        Seed for random number generator, when using variant="ar".

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
    if variant not in (None, 'a', 'ar'):
        raise ValueError("Invalid variant name")

    U, S, V = fast_svd(X, n_comp)
    W, H = np.zeros(U.shape), np.zeros(V.shape)

    # The leading singular triplet is non-negative
    # so it can be used as is for initialization.
    W[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
    H[0, :] = np.sqrt(S[0]) * np.abs(V[0, :])

    for j in xrange(1, n_comp):
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

    W[W < eps] = 0
    H[H < eps] = 0

    if variant == "a":
        avg = X.mean()
        W[W == 0] = avg
        H[H == 0] = avg
    elif variant == "ar":
        rnd = np.random.mtrand.RandomState(seed)
        W[W == 0] = abs(rnd.randn(len(W[W == 0])))
        H[H == 0] = abs(rnd.randn(len(H[H == 0])))

    return W, H


class CRO():
    """
    Closeness to Rank One Hierarchical Clustering

    Model that clusters the columns of a matrix into a given number of clusters
    by joining at each step the clusters with the largest CRO value.

    Parameters
    ----------
        n_comp: int or None
            Target number of components (clusters) to extract.
            Defaults to 1

        epsilon: double or None
            Padding value for output matrices. The value influences sparsity
            and NMF convergence time, but does not influence the performance
            of the initialization.

    Attributes
    ----------
        components_, data_:
            Output matrices to be used for NMF initialization
        clusters:
            List of clusters extracted, each one having:
                size: int, the number of columns in the cluster
                data: array, the submatrix corresponding to the cluster
                svd: tuple (u, s, v), the rank-one approximation of the data

    Examples
    --------
    The example in the paper outputs the given result
    >>> X = [[1, 2, 0, 3, 1],
    ...      [0, 0, 1, 0, 0],
    ...      [0, 0, 1, 0, 0],
    ...      [2, 4, 2, 6, 3],
    ...      [3, 6, 4, 9, 4],
    ...      [0, 0, 2, 0, 0]]
    >>> X = np.array(X)
    >>> model = CRO(n_comp=1)
    >>> model.fit(X)
    >>> model.clusters[0].data
    array([[1, 2, 3, 1, 0],
           [0, 0, 0, 0, 1],
           [0, 0, 0, 0, 1],
           [2, 4, 6, 3, 2],
           [3, 6, 9, 4, 4],
           [0, 0, 0, 0, 2]])

    Notes
    -----
    See the paper "A method of initialization for nonnegative matrix
    factorization" by Yong-Deok Kim and Seungjin Choi, available at:
    http://www.postech.ac.kr/~seungjin/publications/icassp07_ydkim.pdf

    """
    def __init__(self, n_comp=1, epsilon=1e-5):
        """
        Initializes the CRO model
        """
        self.n_comp = n_comp
        self.clusters = []
        self.epsilon = epsilon

    def fit(self, X):
        """
        Clusters the matrix X
        """

        # Each column is an individual clusters
        n_samples, n_features = X.shape
        for col in X.T:
            norm = np.linalg.norm(col)
            svd = (col / norm, norm, np.ones(1))
            self.clusters.append(Bunch(size=1,
                                       data=col,
                                       svd=svd))

        for step in xrange(n_features, self.n_comp, -1):
            cros = {}
            for i in xrange(step - 1):
                for j in xrange(i + 1, step):
                    cro = self.pairwise_cro(self.clusters[i], self.clusters[j])
                    cros[i, j] = cro

            pair = dict_argmax(cros)
            t, s = self.clusters[pair[0]], self.clusters[pair[1]]
            self.clusters[pair[0]] = self._merge(t, s)
            del self.clusters[pair[1]]

        self.data_ = np.zeros((n_samples, 0))
        self.components_ = np.zeros((self.n_comp, n_features))
        j = 0
        for i, cl in enumerate(self.clusters):
            self.data_ = np.c_[self.data_, cl.svd[1] * cl.svd[0]]
            self.components_[i, j:j + cl.size] = cl.svd[2]
            j += cl.size

        self.data_[self.data_ == 0] += self.epsilon
        self.components_[self.components_ == 0] += self.epsilon

    def _merge(self, target, source):
        """
        Merges two clusters and updates the rank-one approximation
        """

        size = target.size + source.size
        data = np.c_[target.data, source.data]
        L = np.c_[target.svd[0] * target.svd[1],
                  source.svd[0] * source.svd[1]]
        R = np.r_[np.c_[target.svd[2], np.zeros(np.shape(target.svd[2]))],
                  np.c_[np.zeros(np.shape(source.svd[2])), source.svd[2]]]
        _, S, V = np.linalg.svd(np.dot(L.T, L))
        S = np.sqrt(S) + 1e-8
        assert (S != 0).all()
        U = np.atleast_2d(np.dot(np.dot(L, V), np.diag(1 / S)))

        svd = U[:, 0], S[0], np.dot(R, V[0])
        return Bunch(size=size, data=data, svd=svd)

    def pairwise_cro(self, u, v):
        """
        CRO between two clusters
        """
        #assert(u.size == v.size == 1)
        #X = np.c_[u.data, v.data]
        #sigma = max(np.linalg.eigvals(np.dot(X.T, X)))
        #return sigma / np.linalg.norm(X) ** 2
        result = self._merge(u, v)
        return (result.svd[1] / np.linalg.norm(result.data)) ** 2


def _nls_subproblem_(V, W, Hinit, tol, max_iter):
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

    tol:
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
        if proj_gradient < tol:
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

    initial: 'nndsvd', int or RandomState
        Method used to initialize the procedure.
        Default: 'nndsvd'
        Valid options:
            'nndsvd' for SVD-based initialization,
            'cro' for CRO-based initialization,
            int seed or RandomState for non-negative random matrices

    sparseness: string or None
        'data' or 'components', where to enforce sparsity
        Default: None

    beta: double
        Degree of sparseness, if sparseness is not None
        Default: 1

    eta: double
        Degree of correctness to mantain, if sparsity is not None
        Default: 0.1

    tol: double
        Tolerance value used in stopping conditions.
        Default: 1e-4

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
    >>> model = NMF(n_comp=2, initial=0)
    >>> model.fit(X) #doctest: +ELLIPSIS
    NMF(nls_max_iter=2000, n_comp=2, max_iter=100, sparseness=None,
      initial=<mtrand.RandomState object at 0x...>, beta=1, eta=0.1,
      tol=0.0001)
    >>> model.components_
    array([[ 0.77032744,  0.38526873],
           [ 0.11118662,  0.38228063]])
    >>> model.reconstruction_err_
    0.0074679497834630824
    >>> model = NMF(n_comp=2, initial=0, sparseness='components')
    >>> model.fit(X) #doctest: +ELLIPSIS
    NMF(nls_max_iter=2000, n_comp=2, max_iter=100, sparseness='components',
      initial=<mtrand.RandomState object at 0x...>, beta=1, eta=0.1,
      tol=0.0001)
    >>> model.components_
    array([[ 1.67481991, -0.        ],
           [ 0.29614922,  0.4681982 ]])
    >>> model.reconstruction_err_
    0.51328426689311935

    Notes
    -----
    This implements C.-J. Lin. Projected gradient methods
    for non-negative matrix factorization. Neural
    Computation, 19(2007), 2756-2779.
    http://www.csie.ntu.edu.tw/~cjlin/nmf/

    """

    def __init__(self, n_comp=None, initial="nndsvd", sparseness=None, beta=1,
                 eta=0.1, tol=1e-4, max_iter=100, nls_max_iter=2000):
        self.n_comp = n_comp
        self.initial = initial
        self.tol = tol
        if sparseness not in (None, 'data', 'components'):
            raise ValueError('Invalid sparsity target')
        self.sparseness = sparseness
        self.beta = beta
        self.eta = eta
        self.max_iter = max_iter
        self.nls_max_iter = nls_max_iter

    def fit(self, X):
        """
        Fit the model to the data
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
        elif self.initial == "nndsvd":
            W, H = _initialize_nmf_(X, self.n_comp)
        elif self.initial == "cro":
            m = CRO(self.n_comp)
            m.fit(X.T)
            W, H = np.abs(m.components_.T), np.abs(m.data_.T)
        else:
            raise ValueError("Invalid value for initial parameter.")

        gradW = np.dot(W, np.dot(H, H.T)) - np.dot(X, H.T)
        gradH = np.dot(np.dot(W.T, W), H) - np.dot(W.T, X)
        init_grad = norm(np.r_[gradW, gradH.T])
        tolW = max(0.001, self.tol) * init_grad  # why max?
        tolH = tolW

        for iter in xrange(1, self.max_iter + 1):
            # stopping condition
            # as discussed in paper
            proj_norm = norm(np.r_[gradW[np.logical_or(gradW < 0, W > 0)],
                                   gradH[np.logical_or(gradH < 0, H > 0)]])
            if proj_norm < self.tol * init_grad:
                break

            # update W
            if self.sparseness == None:
                W, gradW, iterW = _nls_subproblem_(X.T, H.T, W.T, tolW,
                                                   self.nls_max_iter)
            elif self.sparseness == 'data':
                W, gradW, iterW = _nls_subproblem_(
                        r_[X.T, zeros((1, n_features))],
                        r_[H.T, sqrt(self.beta) * ones((1, self.n_comp))],
                        W.T, tolW, self.nls_max_iter)
            elif self.sparseness == 'components':
                W, gradW, iterW = _nls_subproblem_(
                        r_[X.T, zeros((self.n_comp, n_features))],
                        r_[H.T, sqrt(self.eta) * eye(self.n_comp)],
                        W.T, tolW, self.nls_max_iter)

            W = W.T
            gradW = gradW.T
            if iterW == 1:
                tolW = 0.1 * tolW

            # update H
            if self.sparseness == None:
                H, gradH, iterH = _nls_subproblem_(X, W, H, tolH,
                                                   self.nls_max_iter)
            elif self.sparseness == 'data':
                H, gradH, iterH = _nls_subproblem_(
                        r_[X, zeros((self.n_comp, n_samples))],
                        r_[W, sqrt(self.eta) * eye(self.n_comp)],
                        H, tolH, self.nls_max_iter)
            elif self.sparseness == 'components':
                H, gradH, iterH = _nls_subproblem_(
                        r_[X, zeros((1, n_samples))],
                        r_[W, sqrt(self.beta) * ones((1, self.n_comp))],
                        H, tolH, self.nls_max_iter)
            if iterH == 1:
                tolH = 0.1 * tolH
            sparH = _sparseness_(H.flatten())
            sparW = _sparseness_(W.flatten())
            self.sparseness_ = (sparH, sparW)
            self.components_ = H.T
            self.data_ = W
            self.reconstruction_err_ = norm(X - np.dot(W, H))
        if iter == self.max_iter:
            warnings.warn("Iteration limit reached during fit")
        return self

    def transform(self, X):
        """
        Transform the data X according to the model
        """
        from scipy.optimize import nnls
        X = np.atleast_2d(X)
        H = np.zeros((X.shape[0], self.n_comp))
        for j in xrange(0, X.shape[0]):
            H[j, :], _ = nnls(self.components_, X[j, :])
        return H
