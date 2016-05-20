import numpy as np
from ..utils.validation import check_is_fitted
from ..utils.extmath import randomized_svd
from ..utils.validation import check_array
from ..base import BaseEstimator, TransformerMixin
from numpy.linalg import norm
import warnings


def _soft_thresh(X, threshold):
    "Apply soft thresholding to an array"

    sign = np.sign(X)
    return np.multiply(sign, np.maximum(np.abs(X) - threshold, 0))


def norm_1(X):
    return np.sum(np.abs(X))


def _sv_thresh(X, threshold, num_svalue):
    """
    Perform singular value thresholding.

    Parameters
    ---------
    X : array of shape [n_samples, n_features]
        The input array.
    threshold : float
        The threshold for the singualar values.
    num_svalue : int
        The number of singular values to compute.

    Returns
    -------
    X_thresh : array of shape [n_samples, n_features]
        The output after performing singular value thresholding.
    grater_sv : int
        The number of singular values of `X` which were greater than
        `threshold`
    (U, s, V): tuple
        The singular value decomposition
    """
    m, n = X.shape
    U, s, V = randomized_svd(X, num_svalue)
    greater_sv = np.count_nonzero(s > threshold)
    s = _soft_thresh(s, threshold)
    S = np.diag(s)
    X_thresh = np.dot(U, np.dot(S, V))
    return X_thresh, greater_sv, (U, s, V)


def rpca(M, lam=None, mu=None, max_iter=1000, eps_primal=1e-7, eps_dual=1e-5,
         rho=1.6, initial_sv=10, max_mu=1e6, verbose=False):
    """Implements the Robust PCA algorithm via Principal Component Pursuit [1]_

    The Robust PCA algorithm minimizes

    .. math:: \\lVert L \\rVert_* + \\lambda \\lVert S \\rVert_1
    subject to

    .. math:: M = L + S
    where :math:`\\lVert X \\rVert_1` is the sum of absolute values of the
    matrix `X`.

    The algorithm used for optimization is the "Inexact ALM" method specified
    in [2]_

    Parameters
    ----------
    M : array-like, shape (n_samples, n_features)
        The input matrix.

    lam : float, optional
        The importance given to sparsity. Increasing this parameter will yeild
        a sparser `S`. If not given it is set to :math:`\\frac{1}{\\sqrt{n}}`
        where ``n = max(n_samples, n_features)``.

    mu : float, optional
        The initial value of the penalty parameter in the Augmented Lagrangian
        Multiplier (ALM) algorithm. This controls how much attention is given
        to the constraint in each iteration of the optimization problem.

    max_iter : int, optional
        The maximum number of iterations the optimization algortihm will run
        for.

    eps_primal : float, optional
        The threshold for the primal error in the convex optimization problem.
        If the primal and the dual error fall below ``eps_primal`` and
        ``eps_dual`` respectively, the algorithm converges.

    eps_dual :  float, optinal
         The theshold for the dual error in the convex optimzation problem.

    rho : float, optional
        The ratio of the paramter ``mu`` between two successive iterations.
        For each iteration ``mu`` is updated as ``mu = mu*rho``.

    initial_sv : int, optional
        The number of singular values to compute during the first iteration.

    rho_max : float, optional
        The maximum value that ``rho`` is allowed to take.

    verbose : bool, optional
        Whether to print convergence statistics during each iteration.


    Returns
    -------
    L : array, shape (n_samples, n_features)
        The low rank component.

    S : array, shape (n_samples, n_features)
        The sparse component.

    (U, s, Vt) : tuple of arrays
        The singular value decomposition of the ``L``

    n_iter : int
        The number of iterations taken to converge.

    References
    ----------
    .. [1] : Emmanuel J. Cand`es 1,2, Xiaodong Li, Yi Ma, John Wright4, 2009:
             Robust Principal Component Analysis?
    .. [2] : Zhouchen Lin, Minming Chen, Yi Ma, 2013 : The Augmented Lagrange
             Multiplier Method for Exact Recovery of Corrupted Low-Rank
             Matrices
    """

    # See http://arxiv.org/pdf/1009.5055v3.pdf

    # This implementation follows Algorithm 5 from the paper with minor
    # modifications

    if lam is None:
        lam = 1.0/np.sqrt(max(M.shape))

    d = min(M.shape)

    # See "Choosing Parameters" paragraph in section 4
    mu = 1.25/norm(M, 2)

    # The sparse matrix
    S = np.zeros_like(M)

    # The low rank matrix
    L = np.zeros_like(M)

    # See equation 10
    J = min(norm(M, 2), np.max(np.abs(M)))
    Y = M/J

    M_fro_norm = norm(M, 'fro')

    # This variable tried to predict how many singular values will be required.
    sv = initial_sv

    for iter_ in range(max_iter):
        # See Section 4, paragraph "Order of Updating A and E" to see why
        # `S` iterate is computed before `L` ierate.
        S_old = S
        S = _soft_thresh(M - L + (Y/mu), lam/mu)
        L, svp, (U, s, V) = _sv_thresh(M - S + (Y/mu), 1/mu, sv)
        Y = Y + mu*(M - L - S)

        mu_old = mu
        mu = rho*mu
        mu = min(mu, max_mu)

        # See Equation 18
        if svp < sv:
            sv = svp + 1
        else:
            sv = svp + int(round(0.05*d))

        sv = min(sv, M.shape[0], M.shape[1])

        primal_error = norm(M - L - S, 'fro')/M_fro_norm
        dual_error = mu_old*norm(S - S_old, 'fro')/M_fro_norm

        if verbose:
            print('rpca: Iteration %d - Primal Error = %e Dual Error = %e' %
                  (iter_, primal_error, dual_error))
        if primal_error < eps_primal and dual_error < eps_dual:
            break

    if iter_ >= max_iter:
        warnings.warn('rpca: Failed to converge within %d iterations' %
                      max_iter)

    n_iter = iter_
    return L, S, (U, s, V), n_iter


class RobustPCA(BaseEstimator, TransformerMixin):
    """Implements the Robust PCA algorithm via Principal Component Pursuit [1]_

    Robust PCA is designed to overcome the susceptibility to outliers of the
    classical :py:class:`sklearn.decomposition.PCA` algorithm. The Robust PCA
    algorithm tries to seperate out the outliers in the data into the ``S``
    matrix while ``L`` contains the low-rank approximation.

    The Robust PCA algorithm minimizes

    .. math:: \\lVert L \\rVert_* + \\lambda \\lVert S \\rVert_1
    subject to

    .. math:: M = L + S
    where :math:`\\lVert X \\rVert_1` is the sum of absolute values of the
    matrix `X`.

    The algorithm used for optimization is the "Inexact ALM" method specified
    in [2]_

    Parameters
    ----------
    M : array-like, shape (n_samples, n_features)
        The input matrix.

    lam : float, optional
        The importance given to sparsity. Increasing this parameter will yeild
        a sparser `S`. If not given it is set to :math:`\\frac{1}{\\sqrt{n}}`
        where ``n = max(n_samples, n_features)``.

    mu : float, optional
        The initial value of the penalty parameter in the Augmented Lagrangian
        Multiplier (ALM) algorithm. This controls how much attention is given
        to the constraint in each iteration of the optimization problem.

    max_iter : int, optional
        The maximum number of iterations the optimization algortihm will run
        for.

    eps_primal : float, optional
        The threshold for the primal error in the convex optimization problem.
        If the primal and the dual error fall below ``eps_primal`` and
        ``eps_dual`` respectively, the algorithm converges.

    eps_dual :  float, optinal
         The theshold for the dual error in the convex optimzation problem.

    rho : float, optional
        The ratio of the paramter ``mu`` between two successive iterations.
        For each iteration ``mu`` is updated as ``mu = mu*rho``.

    initial_sv : int, optional
        The number of singular values to compute during the first iteration.

    rho_max : float, optional
        The maximum value that ``rho`` is allowed to take.

    verbose : bool, optional
        Whether to print convergence statistics during each iteration.

    Attributes
    ----------
    n_components_ : int
        The rank of the low-rank approximation of the data. This is same as the
        rank of `low_rank_`.

    components_ : array, [n_components, n_features]
        The principal axes of the low rank subspace that the algorithm chose.

    low_rank_ : array, [n_samples, n_features]
        The low rank approximation of the fitted data. Same as ``L`` in the
        optimization problem.

    n_ter_ : int
        The number of iterations taken to converge.

    """
    def __init__(self, lam=None, mu=None, max_iter=1000, eps_primal=1e-7,
                 eps_dual=1e-5, rho=1.6, initial_sv=10, max_mu=1e6,
                 verbose=False):

        self.lam = lam
        self.mu = mu
        self.max_iter = max_iter
        self.eps_primal = eps_primal
        self.eps_dual = eps_dual
        self.rho = rho
        self.initial_sv = initial_sv
        self.max_mu = max_mu
        self.verbose = verbose

    def fit(self, X, y=None):
        """Fit the model with ``X``.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """

        X = check_array(X, dtype=np.float)
        L, S, (U, s, Vt), self.n_iter_ = rpca(X, self.lam, self.mu,
                                              self.max_iter, self.eps_primal,
                                              self.eps_dual, self.rho,
                                              self.initial_sv, self.max_mu,
                                              self.verbose)
        self.low_rank_ = L
        r = np.count_nonzero(s)
        self.n_components_ = r
        self.components_ = Vt[:r]

        return self

    def transform(self, X):
        """ Reduce dimensions of ``X`` by projecting into a low rank subspace

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
            The projection of ``X`` along the low-rank subspace
        """

        check_is_fitted(self, 'components_')
        X = check_array(X, dtype=np.float)
        return np.dot(X, self.components_.T)

    def inverse_transform(self, X):
        """Transform data back to its original space.

        Returns an array X_original whose transform would be X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            New data, where n_samples in the number of samples
            and n_components is the number of components.

        Returns
        -------
        X_original: array, [n_samples, n_features]
            ``X`` projected onto the original space.

        Notes
        -----
        Applying :py:func:`transform` and :py:func:`inverse_transform`
        successively to data should return a low-rank approximation of
        the data in the original feature space.
        """

        check_is_fitted(self, 'components_')
        X = check_array(X, dtype=np.float)
        return np.dot(X, self.components_)
