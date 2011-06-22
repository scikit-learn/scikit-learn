""" Matrix factorization with Sparse PCA
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD

import time
import sys
import itertools
from math import sqrt, floor, ceil

import numpy as np
from numpy.lib.stride_tricks import as_strided
from scipy import linalg

from ..utils.extmath import fast_svd
from ..linear_model import Lasso, lars_path
from ..externals.joblib import Parallel, delayed
from ..base import BaseEstimator, TransformerMixin


##################################
# Utilities to spread load on CPUs
# XXX: where should this be?
def _gen_even_slices(n, n_packs):
    """Generator to create n_packs slices going up to n.

    Examples
    ========

    >>> list(_gen_even_slices(10, 1))
    [slice(0, 10, None)]
    >>> list(_gen_even_slices(10, 10)) #doctest: +ELLIPSIS
    [slice(0, 1, None), slice(1, 2, None), ..., slice(9, 10, None)]
    >>> list(_gen_even_slices(10, 5)) #doctest: +ELLIPSIS
    [slice(0, 2, None), slice(2, 4, None), ..., slice(8, 10, None)]
    >>> list(_gen_even_slices(10, 3))
    [slice(0, 4, None), slice(4, 7, None), slice(7, 10, None)]

    """
    start = 0
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            yield slice(start, end, None)
            start = end


def cpu_count():
    """ Return the number of CPUs.
    """
    # XXX: should be in joblib
    try:
        import multiprocessing
    except ImportError:
        return 1
    return multiprocessing.cpu_count()


# a short preview of what will be in fabian's pull request
def _ridge_regression(X, y, alpha):
    n_samples, n_features = np.shape(X)
    if n_features > n_samples:
        # kernel ridge
        # w = X.T * inv(X X^t + alpha*Id) y
        A = np.dot(X, X.T)
        A.flat[::n_samples + 1] += alpha
        return np.dot(X.T, linalg.solve(A, y, sym_pos=True, overwrite_a=True))
    else:
        # ridge
        # w = inv(X^t X + alpha*Id) * X.T y
        A = np.dot(X.T, X)
        A.flat[::n_features + 1] += alpha
        return linalg.solve(A, np.dot(X.T, y), sym_pos=True, overwrite_a=True)


def _update_V(U, Y, alpha, V=None, Gram=None, method='lars', tol=1e-8):
    """ Update the sparse factor V in sparse_pca loop.
    Each column of V is the solution to a Lasso problem.

    Parameters
    ----------
    U: array of shape (n_samples, n_components)
        previous iteration of U

    Y: array of shape (n_samples, n_features)
        data matrix

    alpha: float
        regularization parameter for the Lasso problem

    V: array of shape (n_components, n_features)
        previous iteration of V

    Gram: array of shape (n_features, n_features)
        precomputed Gram matrix, (Y^T * Y)

    method: 'lars' | 'lasso'
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso: uses the stochastic gradient descent method to compute the
            lasso solution (linear_model.Lasso)

    tol: float
        numerical tolerance for Lasso convergence.
        Ignored if `method='lars'`

    """
    n_features = Y.shape[1]
    n_atoms = U.shape[1]
    coef = np.empty((n_atoms, n_features))
    if method == 'lars':
        if Gram is None:
            Gram = np.dot(U.T, U)
        err_mgt = np.seterr()
        np.seterr(all='ignore')
        #alpha = alpha * n_samples
        XY = np.dot(U.T, Y)
        for k in range(n_features):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            _, _, coef_path_ = lars_path(U, Y[:, k], Xy=XY[:, k], Gram=Gram,
                                         alpha_min=alpha, method='lasso')
            coef[:, k] = coef_path_[:, -1]
        np.seterr(**err_mgt)
    else:
        clf = Lasso(alpha=alpha, fit_intercept=False)
        for k in range(n_features):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            if V is not None:
                clf.coef_ = V[:, k]  # Init with previous value of Vk
            clf.fit(U, Y[:, k], max_iter=1000, tol=tol)
            coef[:, k] = clf.coef_
    return coef


def _update_V_parallel(U, Y, alpha, V=None, Gram=None, method='lars',
                       n_jobs=1, tol=1e-8):
    """ Update the sparse factor V in sparse_pca loop by efficiently
    spreading the load over the available cores.

    Parameters
    ----------
    U: array of shape (n_samples, n_components)
        previous iteration of U

    Y: array of shape (n_samples, n_features)
        data matrix

    alpha: float
        regularization parameter for the Lasso problem

    V: array of shape (n_components, n_features)
        previous iteration of V
    Gram: array of shape (n_features, n_features)
        precomputed Gram matrix, (Y^T * Y)

    method: 'lars' | 'lasso'
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso: uses the stochastic gradient descent method to compute the
            lasso solution (linear_model.Lasso)

    n_jobs: int
        number of parallel jobs to run

    tol: float
        numerical tolerance for Lasso convergence.
        Ignored if `method='lars'`

    """
    n_samples, n_features = Y.shape
    if Gram is None:
        Gram = np.dot(U.T, U)
    if n_jobs == 1:
        return _update_V(U, Y, alpha, V=V, Gram=Gram, method=method)
    n_atoms = U.shape[1]
    if V is None:
        V = np.empty((n_atoms, n_features))
    slices = list(_gen_even_slices(n_features, n_jobs))
    V_views = Parallel(n_jobs=n_jobs)(
                delayed(_update_V)(U, Y[:, this_slice], V=V[:, this_slice],
                                alpha=alpha, Gram=Gram, method=method, tol=tol)
                for this_slice in slices)
    for this_slice, this_view in zip(slices, V_views):
        V[:, this_slice] = this_view
    return V


def _update_U(U, Y, V, verbose=False, return_r2=False):
    """ Update the dense factor U in sparse_pca loop in place.

    Parameters
    ----------
    U: array of shape (n_samples, n_components)
        previous iteration of U

    Y: array of shape (n_samples, n_features)
        data matrix

    V: array of shape (n_components, n_features)
        previous iteration of V

    verbose:
        degree of output the procedure will print

    return_r2: bool
        compute and return the residual sum of squares corresponding
        to the computed solution

    """
    n_atoms = len(V)
    n_samples = Y.shape[0]
    R = -np.dot(U, V)  # Residuals, computed 'in-place' for efficiency
    R += Y
    R = np.asfortranarray(R)
    ger, = linalg.get_blas_funcs(('ger',), (U, V))
    for k in xrange(n_atoms):
        # R <- 1.0 * U_k * V_k^T + R
        R = ger(1.0, U[:, k], V[k, :], a=R, overwrite_a=True)
        U[:, k] = np.dot(R, V[k, :].T)
        # Scale Uk
        norm_square_U = np.dot(U[:, k], U[:, k])
        if norm_square_U < 1e-20:
            if verbose == 1:
                sys.stdout.write("+")
                sys.stdout.flush()
            elif verbose:
                print "Adding new random atom"
            U[:, k] = np.random.randn(n_samples)
            # Setting corresponding coefs to 0
            V[k, :] = 0.0
            U[:, k] /= sqrt(np.dot(U[:, k], U[:, k]))
        else:
            U[:, k] /= sqrt(norm_square_U)
            # R <- -1.0 * U_k * V_k^T + R
            R = ger(-1.0, U[:, k], V[k, :], a=R, overwrite_a=True)
    if return_r2:
        R **= 2
        # R is fortran-ordered. For numpy version < 1.6, sum does not
        # follow the quick striding first, and is thus inefficient on
        # fortran ordered data. We take a flat view of the data with no
        # striding
        R = as_strided(R, shape=(R.size, ), strides=(R.dtype.itemsize,))
        R = np.sum(R)
        return U, R
    return U


def sparse_pca(Y, n_atoms, alpha, max_iter=100, tol=1e-8, method='lars',
        n_jobs=1, U_init=None, V_init=None, callback=None, verbose=False):
    """
    Compute sparse matrix decomposition (PCA) with n_atoms components.

    (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || V ||_1
                 (U,V)
                with || U_k ||_2 = 1 for all  0<= k < n_atoms

    Parameters
    ----------
    Y: array of shape (n_samples, n_features)
        data matrix

    n_atoms: int,
        number of sparse atoms to extract

    alpha: int,
        sparsity controlling parameter

    max_iter: int,
        maximum number of iterations to perform

    tol: float,
        tolerance for numerical error

    method: 'lars' | 'lasso',
        method to use for solving the lasso problem

    n_jobs: int,
        number of parallel jobs to run

    U_init: array of shape (n_samples, n_atoms),
    V_init: array of shape (n_atoms, n_features),
        initial values for the decomposition for warm restart scenarios

    callback:
        callable that gets invoked every five iterations

    verbose:
        degree of output the procedure will print

    Returns
    -------
    U: array of shape (n_samples, n_atoms),
    V: array of shape (n_atoms, n_features),
        the solutions to the sparse PCA decomposition
    E: array
        vector of errors at each iteration
    """
    t0 = time.time()
    n_samples = Y.shape[0]
    # Avoid integer division problems
    alpha = float(alpha)

    if n_jobs == -1:
        n_jobs = cpu_count()

    # Init U and V with SVD of Y
    if U_init is not None and V_init is not None:
        U = np.array(U_init, order='F')
        # Don't copy V, it will happen below
        V = V_init
    else:
        U, S, V = linalg.svd(Y, full_matrices=False)
        V = S[:, np.newaxis] * V
    U = U[:, :n_atoms]
    V = V[:n_atoms, :]

    # Fortran-order V, as we are going to access its row vectors
    V = np.array(V, order='F')

    residuals = 0

    def cost_function():
        return 0.5 * residuals + alpha * np.sum(np.abs(V))

    E = []
    current_cost = np.nan

    if verbose == 1:
        print '[sparse_pca]',

    for ii in xrange(max_iter):
        dt = (time.time() - t0)
        if verbose == 1:
            sys.stdout.write(".")
            sys.stdout.flush()
        elif verbose:
            print ("Iteration % 3i "
                "(elapsed time: % 3is, % 4.1fmn, current cost % 7.3f)" %
                    (ii, dt, dt / 60, current_cost))

        # Update V
        V = _update_V_parallel(U, Y, alpha / n_samples, V=V, method=method,
                               n_jobs=n_jobs)

        # Update U
        U, residuals = _update_U(U, Y, V, verbose=verbose, return_r2=True)

        current_cost = cost_function()
        E.append(current_cost)

        if ii > 0:
            dE = E[-2] - E[-1]
            assert(dE >= -tol * E[-1])
            if dE < tol * E[-1]:
                if verbose == 1:
                    # A line return
                    print ""
                elif verbose:
                    print "--- Convergence reached after %d iterations" % ii
                break
        if ii % 5 == 0 and callback is not None:
            callback(locals())

    return U, V, E


def dict_learning(Y, n_atoms, alpha, n_iter=100, return_code=True,
        dict_init=None, callback=None, chunk_size=3, verbose=False,
        shuffle=True, n_jobs=1, coding_method='lars'):
    """
    Online dictionary learning for sparse coding 

    (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || V ||_1
                 (U,V)
                with || U_k ||_2 = 1 for all  0<= k < n_atoms

    """
    t0 = time.time()
    n_samples, n_features = Y.shape
    # Avoid integer division problems
    alpha = float(alpha)

    if n_jobs == -1:
        n_jobs = cpu_count()

    # Init V with SVD of Y
    if dict_init is not None:
        V = dict_init
    else:
        _, S, V = fast_svd(Y, n_atoms)
        V = S[:, np.newaxis] * V
    V = V[:n_atoms, :]
    V = np.ascontiguousarray(V.T)

    if verbose == 1:
        print '[dict_learning]',

    n_batches = floor(float(len(Y)) / chunk_size)
    if shuffle:
        Y_train = Y.copy()
        np.random.shuffle(Y_train)
    else:
        Y_train = Y
    batches = np.array_split(Y_train, n_batches)
    batches = itertools.cycle(batches)

    # The covariance of the dictionary
    A = np.zeros((n_atoms, n_atoms))
    # The data approximation
    B = np.zeros((n_features, n_atoms))

    for ii, this_Y in itertools.izip(xrange(n_iter), batches):
        #this_Y = this_Y.squeeze()
        dt = (time.time() - t0)
        if verbose == 1:
            sys.stdout.write(".")
            sys.stdout.flush()
        elif verbose:
            if verbose > 10 or ii % ceil(100./verbose) == 0:
                print ("Iteration % 3i (elapsed time: % 3is, % 4.1fmn)" % 
                    (ii, dt, dt/60))
        
        this_U = _update_V(V, this_Y.T, alpha)

        # Update the auxiliary variables
        if ii < chunk_size - 1:
            theta = float((ii+1)*chunk_size)
        else:
            theta = float(chunk_size**2 + ii + 1 - chunk_size)
        beta = (theta + 1 - chunk_size)/(theta + 1)

        A *= beta
        A += np.dot(this_U, this_U.T)
        B *= beta
        B += np.dot(this_Y.T, this_U.T)

        # Update dictionary
        V = _update_U(V, B, A, verbose=verbose)
        # XXX: Can the residuals be of any use?

        # Maybe we need a stopping criteria based on the amount of
        # modification in the dictionary
        if callback is not None:
            callback(locals())

    if return_code:
        if verbose > 1:
            print 'Learning code...',
        elif verbose == 1:
            print '|', 
        code = _update_V_parallel(V, Y.T, alpha, n_jobs=n_jobs, 
                    method=coding_method)
        if verbose > 1:
            dt = (time.time() - t0)
            print 'done (total time: % 3is, % 4.1fmn)' % (dt, dt/60)
        return V.T, code.T
        
    return V.T


class SparsePCA(BaseEstimator, TransformerMixin):
    """Sparse Principal Components Analysis (SparsePCA)

    Finds the set of sparse components that can optimally reconstruct the data.
    The amount of sparseness is controllable by the coefficient of the \ell_1
    penalty, given by the parameter alpha.

    Parameters
    ----------
    n_components: int,
        number of sparse atoms to extract

    alpha: int,
        sparsity controlling parameter

    max_iter: int,
        maximum number of iterations to perform

    tol: float,
        tolerance for numerical error

    method: 'lars' | 'lasso',
        method to use for solving the lasso problem

    n_jobs: int,
        number of parallel jobs to run

    U_init: array of shape (n_samples, n_atoms),
    V_init: array of shape (n_atoms, n_features),
        initial values for the decomposition for warm restart scenarios

    verbose:
        degree of verbosity of the printed output

    Attributes
    ----------
    components_: array, [n_components, n_features]
        sparse components extracted from the data

    error_: array
        vector of errors at each iteration

    See also
    --------
    PCA

    """
    def __init__(self, n_components=None, alpha=1, max_iter=1000, tol=1e-8,
                 method='lars', n_jobs=1, U_init=None, V_init=None,
                 verbose=False):
        self.n_components = n_components
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.n_jobs = n_jobs
        self.U_init = U_init
        self.V_init = V_init
        self.verbose = verbose

    def fit_transform(self, X, y=None, **params):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new array-like, shape (n_samples, n_components)
        """
        self._set_params(**params)
        X = np.asanyarray(X)

        U, V, E = sparse_pca(X, self.n_components, self.alpha, tol=self.tol,
                             max_iter=self.max_iter, method=self.method,
                             n_jobs=self.n_jobs, verbose=self.verbose)
        self.components_ = V
        self.error_ = E
        return U

    def fit(self, X, y=None, **params):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self.fit_transform(X, y, **params)
        return self

    def transform(self, X, alpha=0):
        """Apply the projection onto the learned sparse components
        to new data.

        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        alpha: float
            Amount of ridge shrinkage to apply in order to improve conditioning

        Returns
        -------
        X_new array, shape (n_samples, n_components)
            Transformed data
        """

#        if alpha != 0:
#            raise NotImplemented('SparsePCA.transform only does OLS for now')
        # TODO: Ridge regression with controllable shrinkage
        n_samples = len(X)
        U = np.zeros((n_samples, self.n_components))
        for k in xrange(n_samples):
            U[k, :] = _ridge_regression(self.components_.T, X[k, :], alpha)
        # U = linalg.lstsq(self.components_.T, X.T)[0].T
        U /= np.sqrt((U ** 2).sum(axis=0))
        return U
