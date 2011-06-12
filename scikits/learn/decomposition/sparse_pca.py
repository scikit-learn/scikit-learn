import time
import sys
from math import sqrt

import numpy as np
from numpy.lib.stride_tricks import as_strided
from scipy import linalg

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


###########
# sparsePCA
def _update_V(U, Y, V, alpha, Gram=None, method='lars', tol=1e-8):
    """ Update V (dictionary) in sparse_pca loop.
    """
    coef = np.empty_like(V)
    if method == 'lars':
        err_mgt = np.seterr()
        np.seterr(all='ignore')
        #alpha = alpha * n_samples
        XY = np.dot(U.T, Y)
        for k in range(V.shape[1]):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            _, _, coef_path_ = lars_path(U, Y[:, k], Xy=XY[:, k], Gram=Gram,
                                         alpha_min=alpha, method='lasso')
            coef[:, k] = coef_path_[:, -1]
        np.seterr(**err_mgt)
    else:
        clf = Lasso(alpha=alpha, fit_intercept=False)
        for k in range(V.shape[1]):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            clf.coef_ = V[:, k]  # Init with previous value of Vk
            clf.fit(U, Y[:, k], max_iter=1000, tol=tol)
            coef[:, k] = clf.coef_
    return coef


def _update_U(U, Y, V, verbose=False, return_r2=False):
    """ Update U (data) in sparse_pca loop in place.
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
    Compute sparse PCA with n_atoms components.

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

    slices = list(_gen_even_slices(V.shape[1], n_jobs))
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
        Gram = np.dot(U.T, U)
        slices = list(_gen_even_slices(V.shape[1], n_jobs))
        V_views = Parallel(n_jobs=n_jobs)(
                    delayed(_update_V)(U, Y[:, this_slice], V[:, this_slice],
                                    alpha=alpha / n_samples,
                                    Gram=Gram, method=method, tol=tol)
                    for this_slice in slices)
        for this_slice, this_view in zip(slices, V_views):
            V[:, this_slice] = this_view

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


class SparsePCA(BaseEstimator, TransformerMixin):
    """Sparse Principal Components Analysis (SparsePCA)

    Finds the best decomposition of the data matrix with sparse components.

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

    Attributes
    ----------
    components_: array, [n_components, n_features]
        sparse components extracted from the data
    
    error_: array
        vector of errors at each iteration

    """
    def __init__(self, n_components=None, alpha=1, max_iter=1000, tol=1e-8,
                 method='lars', n_jobs=1, U_init=None, V_init=None):
        self.n_components = n_components
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.n_jobs = n_jobs
        self.U_init = U_init
        self.V_init = V_init

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
                             n_jobs=self.n_jobs)
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

    def transform(self, X):
        """Apply the projection onto the learned sparse components
        to new data.

        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.
        
        Returns
        -------
        X_new array, shape (n_samples, n_components)
            Transformed data
        """
        # TODO: Ridge
        U = linalg.lstsq(self.components_.T, X.T)[0].T
        U /= np.sqrt((U ** 2).sum(axis=0))
        return U
