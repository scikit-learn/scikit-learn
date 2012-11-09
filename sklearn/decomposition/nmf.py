""" Non-negative matrix factorization
"""
# Author: Olivier Mangin <olivier.mangin@inria.fr> (KL-NMF implementation)
# Author: Vlad Niculae
#         Lars Buitinck <L.J.Buitinck@uva.nl>
# Author: Chih-Jen Lin, National Taiwan University (original projected gradient
#     NMF implementation)
# Author: Anthony Di Franco (original Python and NumPy port)
# License: BSD


from __future__ import division

import warnings

import numpy as np
from scipy.optimize import nnls
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin
from ..utils import atleast2d_or_csr, check_random_state
from ..utils.extmath import randomized_svd, safe_sparse_dot


def safe_vstack(Xs):
    if any(sp.issparse(X) for X in Xs):
        return sp.vstack(Xs)
    else:
        return np.vstack(Xs)


def _pos(x):
    """Positive part of a vector / matrix"""
    return (x >= 0) * x


def _neg(x):
    """Negative part of a vector / matrix"""
    neg_x = -x
    neg_x *= x < 0
    return neg_x


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

    X: array, [n_samples, n_features]
        The data matrix to be decomposed.

    n_components:
        The number of components desired in the
        approximation.

    variant: None | 'a' | 'ar'
        The variant of the NNDSVD algorithm.
        Accepts None, 'a', 'ar'
        None: leaves the zero entries as zero
        'a': Fills the zero entries with the average of X
        'ar': Fills the zero entries with standard normal random variates.
        Default: None

    eps:
        Truncate all values less then this in output to zero.

    random_state: numpy.RandomState | int, optional
        The generator used to fill in the zeros, when using variant='ar'
        Default: numpy.random

    Returns
    -------

    (W, H):
        Initial guesses for solving X ~= WH such that
        the number of columns in W is n_components.

    Remarks
    -------

    This implements the algorithm described in
    C. Boutsidis, E. Gallopoulos: SVD based
    initialization: A head start for nonnegative
    matrix factorization - Pattern Recognition, 2008

    http://www.cs.rpi.edu/~boutsc/files/nndsvd.pdf
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

    for j in xrange(1, n_components):
        x, y = U[:, j], V[j, :]

        # extract positive and negative parts of column vectors
        x_p, y_p = _pos(x), _pos(y)
        x_n, y_n = _neg(x), _neg(y)

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


def _nls_subproblem(V, W, H_init, tol, max_iter):
    """Non-negative least square solver

    Solves a non-negative least squares subproblem using the
    projected gradient descent algorithm.
    min || WH - V ||_2

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

    Returns
    -------
    H : array-like
        Solution to the non-negative least squares problem.

    grad : array-like
        The gradient.

    n_iter : int
        The number of iterations done by the algorithm.

    """
    if (H_init < 0).any():
        raise ValueError("Negative values in H_init passed to NLS solver.")

    H = H_init
    WtV = safe_sparse_dot(W.T, V, dense_output=True)
    WtW = safe_sparse_dot(W.T, W, dense_output=True)

    # values justified in the paper
    alpha = 1
    beta = 0.1
    for n_iter in xrange(1, max_iter + 1):
        grad = np.dot(WtW, H) - WtV
        proj_gradient = norm(grad[np.logical_or(grad < 0, H > 0)])
        if proj_gradient < tol:
            break

        for inner_iter in xrange(1, 20):
            Hn = H - alpha * grad
            # Hn = np.where(Hn > 0, Hn, 0)
            Hn = _pos(Hn)
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
                    alpha *= beta
            elif not suff_decr or (Hp == Hn).all():
                H = Hp
                break
            else:
                alpha /= beta
                Hp = Hn

    if n_iter == max_iter:
        warnings.warn("Iteration limit reached in nls subproblem.")

    return H, grad, n_iter


class BaseNMF(BaseEstimator, TransformerMixin):

    def __init__(self, n_components=None, init=None, random_state=None):
        self.n_components = n_components
        self.init = init
        self.random_state = random_state

    def _init(self, X):
        n_samples, n_features = X.shape
        init = self.init
        if init is None:
            if self.n_components < n_features:
                init = 'nndsvd'
            else:
                init = 'random'

        if isinstance(init, (int, np.integer, np.random.RandomState)):
            random_state = check_random_state(init)
            init = "random"
            warnings.warn("Passing a random seed or generator as init "
                "is deprecated and will be removed in 0.15. Use "
                "init='random' and random_state instead.", DeprecationWarning)
        else:
            random_state = self.random_state

        if init == 'nndsvd':
            W, H = _initialize_nmf(X, self.n_components)
        elif init == 'nndsvda':
            W, H = _initialize_nmf(X, self.n_components, variant='a')
        elif init == 'nndsvdar':
            W, H = _initialize_nmf(X, self.n_components, variant='ar')
        elif init == "random":
            rng = check_random_state(random_state)
            W = rng.randn(n_samples, self.n_components)
            # we do not write np.abs(W, out=W) to stay compatible with
            # numpy 1.5 and earlier where the 'out' keyword is not
            # supported as a kwarg on ufuncs
            np.abs(W, W)
            H = rng.randn(self.n_components, n_features)
            np.abs(H, H)
        else:
            raise ValueError(
                'Invalid init parameter: got %r instead of one of %r' %
                (init, (None, 'nndsvd', 'nndsvda', 'nndsvdar', 'random')))
        return W, H


class ProjectedGradientNMF(BaseNMF):
    """Non-Negative matrix factorization by Projected Gradient (NMF)

    Parameters
    ----------
    X: {array-like, sparse matrix}, shape = [n_samples, n_features]
        Data the model will be fit to.

    n_components: int or None
        Number of components, if n_components is not set all components
        are kept

    init:  'nndsvd' |  'nndsvda' | 'nndsvdar' | 'random'
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

    sparseness: 'data' | 'components' | None, default: None
        Where to enforce sparsity in the model.

    beta: double, default: 1
        Degree of sparseness, if sparseness is not None. Larger values mean
        more sparseness.

    eta: double, default: 0.1
        Degree of correctness to mantain, if sparsity is not None. Smaller
        values mean larger error.

    tol: double, default: 1e-4
        Tolerance value used in stopping conditions.

    max_iter: int, default: 200
        Number of iterations to compute.

    nls_max_iter: int, default: 2000
        Number of iterations in NLS subproblem.

    random_state : int or RandomState
        Random number generator seed control.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Non-negative components of the data

    `reconstruction_err_` : number
        Frobenius norm of the matrix difference between the
        training data and the reconstructed data from the
        fit produced by the model. ``|| X - WH ||_2``
        Not computed for sparse input matrices because it is
        too expensive in terms of memory.

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
           [-0.        ,  0.4681982 ]])
    >>> model.reconstruction_err_ #doctest: +ELLIPSIS
    0.513...

    Notes
    -----
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
    http://www.cs.rpi.edu/~boutsc/files/nndsvd.pdf

    """

    def __init__(self, n_components=None, init=None, sparseness=None, beta=1,
            eta=0.1, tol=1e-4, max_iter=200, nls_max_iter=2000,
            random_state=None):
        BaseNMF.__init__(self, init=init, n_components=n_components,
                random_state=random_state)
        self.tol = tol
        if sparseness not in (None, 'data', 'components'):
            raise ValueError(
                'Invalid sparseness parameter: got %r instead of one of %r' %
                (sparseness, (None, 'data', 'components')))
        self.sparseness = sparseness
        self.beta = beta
        self.eta = eta
        self.max_iter = max_iter
        self.nls_max_iter = nls_max_iter

    def _update_W(self, X, H, W, tolW):
        n_samples, n_features = X.shape

        if self.sparseness == None:
            W, gradW, iterW = _nls_subproblem(X.T, H.T, W.T, tolW,
                                              self.nls_max_iter)
        elif self.sparseness == 'data':
            W, gradW, iterW = _nls_subproblem(
                    safe_vstack([X.T, np.zeros((1, n_samples))]),
                    safe_vstack([H.T, np.sqrt(self.beta) * np.ones((1,
                                 self.n_components))]),
                    W.T, tolW, self.nls_max_iter)
        elif self.sparseness == 'components':
            W, gradW, iterW = _nls_subproblem(
                    safe_vstack([X.T,
                                 np.zeros((self.n_components, n_samples))]),
                    safe_vstack([H.T, np.sqrt(self.eta)
                                      * np.eye(self.n_components)]),
                    W.T, tolW, self.nls_max_iter)

        return W, gradW, iterW

    def _update_H(self, X, H, W, tolH):
        n_samples, n_features = X.shape

        if self.sparseness == None:
            H, gradH, iterH = _nls_subproblem(X, W, H, tolH,
                                              self.nls_max_iter)
        elif self.sparseness == 'data':
            H, gradH, iterH = _nls_subproblem(
                    safe_vstack([X,
                                 np.zeros((self.n_components, n_features))]),
                    safe_vstack([W, np.sqrt(self.eta)
                                    * np.eye(self.n_components)]),
                    H, tolH, self.nls_max_iter)
        elif self.sparseness == 'components':
            H, gradH, iterH = _nls_subproblem(
                    safe_vstack([X, np.zeros((1, n_features))]),
                    safe_vstack([W, np.sqrt(self.beta) *
                          np.ones((1, self.n_components))]),
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

        n_samples, n_features = X.shape

        if not self.n_components:
            self.n_components = n_features

        W, H = self._init(X)

        gradW = (np.dot(W, np.dot(H, H.T))
                 - safe_sparse_dot(X, H.T, dense_output=True))
        gradH = (np.dot(np.dot(W.T, W), H)
                 - safe_sparse_dot(W.T, X, dense_output=True))
        init_grad = norm(np.r_[gradW, gradH.T])
        tolW = max(0.001, self.tol) * init_grad  # why max?
        tolH = tolW

        for n_iter in xrange(1, self.max_iter + 1):
            # stopping condition
            # as discussed in paper
            proj_norm = norm(np.r_[gradW[np.logical_or(gradW < 0, W > 0)],
                                   gradH[np.logical_or(gradH < 0, H > 0)]])
            if proj_norm < self.tol * init_grad:
                break

            # update W
            W, gradW, iterW = self._update_W(X, H, W, tolW)

            W = W.T
            gradW = gradW.T
            if iterW == 1:
                tolW = 0.1 * tolW

            # update H
            H, gradH, iterH = self._update_H(X, H, W, tolH)

            if iterH == 1:
                tolH = 0.1 * tolH

            self.comp_sparseness_ = _sparseness(H.ravel())
            self.data_sparseness_ = _sparseness(W.ravel())

            if not sp.issparse(X):
                self.reconstruction_err_ = norm(X - np.dot(W, H))

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
        H = np.zeros((X.shape[0], self.n_components))
        for j in xrange(0, X.shape[0]):
            H[j, :], _ = nnls(self.components_.T, X[j, :])
        return H


def _normalize_sum(a, axis=0):
    if axis >= len(a.shape):
        raise ValueError
    return a / np.expand_dims(np.sum(a, axis=axis), axis)


def _scale(matrix, factors, axis=0):
    """Scales line or columns of a matrix.

    Parameters
    ----------
    :param matrix: 2-dimensional array
    :param factors: 1-dimensional array
    :param axis: 0: columns are scaled, 1: lines are scaled
    """
    if not (len(matrix.shape) == 2):
        raise ValueError(
                "Wrong array shape: %s, should have only 2 dimensions."
                % str(matrix.shape)
                )
    if not (axis == 0 or axis == 1):
        raise ValueError('Wrong axis, should be 0 (scaling lines)\
                or 1 (scaling columns).')
    # Transform factors given as columne shaped matrices
    factors = np.squeeze(np.asarray(factors))
    if axis == 1:
        factors = factors[:, np.newaxis]
    return np.multiply(matrix, factors)


def _generalized_KL(x, y, eps=1.e-8):
    return (np.multiply(x, np.log(np.divide(x + eps, y + eps))) - x + y
            ).sum()


def _sparse_dot(a, b, refmat):
    """Computes dot product of a and b on indices where refmat is nonnzero
    and returns sparse csr matrix with same structure than refmat.

    First calls to eliminate_zeros on refmat which might modify the structure
    of refmat.

    Params
    ------
    a, b: dense arrays
    refmat: sparse matrix

    Dot product of a and b must have refmat's shape.
    """
    refmat.eliminate_zeros()
    ii, jj = refmat.nonzero()
    dot_vals = np.multiply(a[ii, :], b.T[jj, :]).sum(axis=1)
    c = sp.coo_matrix((dot_vals, (ii, jj)), shape=refmat.shape)
    return c.tocsr()


class KLdivNMF(BaseNMF):
    """Non negative factorization with Kullback Leibler divergence cost.

    Parameters
    ----------
    X: {array-like, sparse matrix}, shape = [n_samples, n_features]
        Data the model will be fit to.

    n_components: int or None
        Number of components, if n_components is not set all components
        are kept

    init:  'nndsvd' |  'nndsvda' | 'nndsvdar' | 'random'
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

    tol: double, default: 1e-4
        Tolerance value used in stopping conditions.

    max_iter: int, default: 200
        Number of iterations to compute.

    subit: int, default: 10
        Number of sub-iterations to perform on W (resp. H) before switching
        to H (resp. W) update.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Non-negative components of the data

    random_state : int or RandomState
        Random number generator seed control.

    Examples
    --------

    >>> import numpy as np
    >>> X = np.array([[1,1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from sklearn.decomposition import KLdivNMF
    >>> model = KLdivNMF(n_components=2, init='random', random_state=0)
    >>> model.fit(X) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    KLdivNMF(eps=1e-08, init='random', max_iter=200, n_components=2,
            random_state=0, subit=10, tol=1e-06)
    >>> model.components_
    array([[ 0.50328164,  0.49671836],
           [ 0.93609814,  0.06390186]])

    Notes
    -----
    This implements

    Lee D. D., Seung H. S., Learning the parts of objects by non-negative
      matrix factorization. Nature, 1999
    """

    def __init__(self, n_components=None, init='random', tol=1e-6,
            max_iter=200, eps=1.e-8, subit=10, random_state=None):
        BaseNMF.__init__(self, init=init, n_components=n_components,
                random_state=random_state)
        self.tol = tol
        self.max_iter = max_iter
        self.eps = eps
        # Only for gradient updates
        self.subit = subit

    def fit_transform(self, X, y=None, weights=1., _fit=True,
            return_errors=False, scale_W=False):
        """Learn a NMF model for the data X and returns the transformed data.

        This is more efficient than calling fit followed by transform.

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be decomposed

        weights: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Weights on the cost function used as coefficients on each
            element of the data. If smaller dimension is provided, standard
            numpy broadcasting is used.

        return_errors: boolean
            if True, the list of reconstruction errors along iterations is
            returned

        scale_W: boolean (default: False)
            Whether to force scaling of W during updates. This is only relevant
            if components are normalized.

        _fit: if True (default), update the model, else only compute transform

        Returns
        -------
        data: array, [n_samples, n_components]
            Transformed data

        or (data, errors) if return_errors
        """
        X = atleast2d_or_csr(X)
        check_non_negative(X, "NMF.fit")

        n_samples, n_features = X.shape

        if not self.n_components:
            self.n_components = n_features

        W, H = self._init(X)

        if _fit:
            self.components_ = H

        prev_error = np.Inf
        tol = self.tol * n_samples * n_features

        if return_errors:
            errors = []

        for n_iter in xrange(1, self.max_iter + 1):
            # Stopping condition
            error = self.error(X, W, self.components_, weights=weights)
            if prev_error - error < tol:
                break
            prev_error = error

            if return_errors:
                errors.append(error)

            W = self._update(X, W)

        if n_iter == self.max_iter:
            warnings.warn("Iteration limit reached during fit")

        if return_errors:
            return W, errors
        else:
            return W

    def _update(self, X, W, _fit=True, scale_W=False, eps=1.e-8):
        """Perform one update iteration.

        Updates components if _fit and returns updated coefficients.

        Params:
        -------
            _fit: boolean (default: True)
                Whether to update components.

            scale_W: boolean (default: False)
                Whether to force scaling of W. This is only relevant if
                components are normalized.
        """
        if scale_W:
            # This is only relevant if components are normalized.
            # Not always usefull but might improve convergence speed:
            # Scale W lines to have same sum than X lines
            W = _scale(_normalize_sum(W, axis=1), X.sum(axis=1), axis=1)
        Q = self._Q(X, W, self.components_, eps=eps)
        # update W
        W = self._updated_W(X, W, self.components_, Q=Q)
        if _fit:
            # update H
            self.components_ = self._updated_H(X, W, self.components_, Q=Q)
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

    def transform(self, X, **params):
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
        params['_fit'] = False
        return self.fit_transform(X, **params)

    # Helpers for beta divergence and related updates

    # Errors and performance estimations

    def error(self, X, W, H=None, weights=1., eps=1.e-8):
        X = atleast2d_or_csr(X)
        if H is None:
            H = self.components_
        if sp.issparse(X):
            WH = _sparse_dot(W, H, X)
            # Avoid computing all values of WH to get their sum
            WH_sum = np.sum(np.multiply(np.sum(W, axis=0), np.sum(H, axis=1)))
            return (np.multiply(
                X.data,
                np.log(np.divide(X.data + eps, WH.data + eps))
                )).sum() - X.data.sum() + WH_sum
        else:
            return _generalized_KL(X, np.dot(W, H))

    # Projections

    def scale(self, W, H, factors):
        """Scale W columns and H rows inversely, according to the given
        coefficients.
        """
        factors = np.array(factors)[np.newaxis, :]
        s_W = W * (factors.T + self.eps)
        s_H = H / (factors.T + self.eps)
        # TODO switch to scale method (11/09/2012)
        return s_W, s_H

    # Update rules

    @classmethod
    def _Q(cls, X, W, H, eps=1.e-8):
        """Computes X / (WH) where / is elementwise and WH is a matrix product.
        """
        # X should be at least 2D or csr
        if sp.issparse(X):
            WH = _sparse_dot(W, H, X)
            WH.data = (X.data + eps) / (WH.data + eps)
            return WH
        else:
            return np.divide(X + eps, np.dot(W, H) + eps)

    @classmethod
    def _updated_W(cls, X, W, H, weights=1., Q=None, eps=1.e-8):
        if Q is None:
            Q = cls._Q(X, W, H, eps=eps)
        W = np.multiply(W, safe_sparse_dot(Q, H.T))
        return W

    @classmethod
    def _updated_H(cls, X, W, H, weights=1., Q=None, eps=1.e-8):
        if Q is None:
            Q = cls._Q(X, W, H, eps=eps)
        H = np.multiply(H, safe_sparse_dot(W.T, Q))
        H = _normalize_sum(H, axis=1)
        return H


class NMF(ProjectedGradientNMF):
    __doc__ = ProjectedGradientNMF.__doc__
    pass
