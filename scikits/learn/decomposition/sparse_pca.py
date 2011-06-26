""" Matrix factorization with Sparse PCA
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD

import time
import sys

from math import sqrt

import numpy as np
from numpy.lib.stride_tricks import as_strided
from scipy import linalg

from ..linear_model import Lasso, lars_path
from ..externals.joblib import Parallel, delayed, cpu_count
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


def _update_code(dictionary, Y, alpha, code=None, Gram=None, method='lars',
                 tol=1e-8):
    """ Update the sparse code factor in sparse_pca loop.
    Each column of the result is the solution to a Lasso problem.

    Parameters
    ----------
    dictionary: array of shape (n_samples, n_components)
        dictionary against which to optimize the sparse code

    Y: array of shape (n_samples, n_features)
        data matrix

    alpha: float
        regularization parameter for the Lasso problem

    code: array of shape (n_components, n_features)
        previous iteration of the sparse code

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
    n_atoms = dictionary.shape[1]
    new_code = np.empty((n_atoms, n_features))
    if method == 'lars':
        if Gram is None:
            Gram = np.dot(dictionary.T, dictionary)
        err_mgt = np.seterr()
        np.seterr(all='ignore')
        #alpha = alpha * n_samples
        XY = np.dot(dictionary.T, Y)
        for k in range(n_features):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            _, _, coef_path_ = lars_path(dictionary, Y[:, k], Xy=XY[:, k],
                                    Gram=Gram, alpha_min=alpha, method='lasso')
            new_code[:, k] = coef_path_[:, -1]
        np.seterr(**err_mgt)
    else:
        clf = Lasso(alpha=alpha, fit_intercept=False)
        for k in range(n_features):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            if code is not None:
                clf.coef_ = code[:, k]  # Init with previous value of Vk
            clf.fit(dictionary, Y[:, k], max_iter=1000, tol=tol)
            new_code[:, k] = clf.coef_
    return new_code


def _update_code_parallel(dictionary, Y, alpha, code=None, Gram=None,
                          method='lars', n_jobs=1, tol=1e-8):
    """ Update the sparse factor V in sparse_pca loop by efficiently
    spreading the load over the available cores.

    Parameters
    ----------
    dictionary: array of shape (n_samples, n_components)
        dictionary against which to optimize the sparse code

    Y: array of shape (n_samples, n_features)
        data matrix

    alpha: float
        regularization parameter for the Lasso problem

    code: array of shape (n_components, n_features)
        previous iteration of the sparse code

    Gram: array of shape (n_features, n_features)
        precomputed Gram matrix, (Y^T * Y)

    method: 'lars' | 'lasso'
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso: uses the stochastic gradient descent method to compute the
            lasso solution (linear_model.Lasso)

    n_jobs: int
        number of parallel jobs to run

    tol: float
        numerical tolerance for coordinate descent Lasso convergence.
        Only used if `method='lasso`.

    """
    n_samples, n_features = Y.shape
    n_atoms = dictionary.shape[1]
    if Gram is None:
        Gram = np.dot(dictionary.T, dictionary)
    if n_jobs == 1:
        return _update_code(dictionary, Y, alpha, code=code, Gram=Gram,
                            method=method)
    if code is None:
        code = np.empty((n_atoms, n_features))
    slices = list(_gen_even_slices(n_features, n_jobs))
    code_views = Parallel(n_jobs=n_jobs)(
                delayed(_update_code)(dictionary, Y[:, this_slice], 
                                      code=code[:, this_slice], alpha=alpha,
                                      Gram=Gram, method=method, tol=tol)
                for this_slice in slices)
    for this_slice, this_view in zip(slices, V_views):
        code[:, this_slice] = this_view
    return code


def _update_dict(dictionary, Y, code, verbose=False, return_r2=False):
    """ Update the dense dictionary factor in place. 

    Parameters
    ----------
    dictionary: array of shape (n_samples, n_components)
        value of the dictionary at the previous iteration

    Y: array of shape (n_samples, n_features)
        data matrix

    code: array of shape (n_components, n_features)
        sparse coding of the data against which to optimize the dictionary

    verbose:
        degree of output the procedure will print

    return_r2: bool
        whether to compute and return the residual sum of squares corresponding
        to the computed solution

    """
    n_atoms = len(code)
    n_samples = Y.shape[0]
    # Residuals, computed 'in-place' for efficiency
    R = -np.dot(dictionary, code)
    R += Y
    R = np.asfortranarray(R)
    ger, = linalg.get_blas_funcs(('ger',), (dictionary, code))
    for k in xrange(n_atoms):
        # R <- 1.0 * U_k * V_k^T + R
        R = ger(1.0, dictionary[:, k], code[k, :], a=R, overwrite_a=True)
        dictionary[:, k] = np.dot(R, code[k, :].T)
        # Scale k'th atom
        atom_norm_square = np.dot(dictionary[:, k], dictionary[:, k])
        if atom_norm_square < 1e-20:
            if verbose == 1:
                sys.stdout.write("+")
                sys.stdout.flush()
            elif verbose:
                print "Adding new random atom"
            dictionary[:, k] = np.random.randn(n_samples)
            # Setting corresponding coefs to 0
            code[k, :] = 0.0
            dictionary[:, k] /= sqrt(np.dot(dictionary[:, k],
                                            dictionary[:, k]))
        else:
            dictionary[:, k] /= sqrt(atom_norm_square)
            # R <- -1.0 * U_k * V_k^T + R
            R = ger(-1.0, dictionary[:, k], code[k, :], a=R, overwrite_a=True)
    if return_r2:
        R **= 2
        # R is fortran-ordered. For numpy version < 1.6, sum does not
        # follow the quick striding first, and is thus inefficient on
        # fortran ordered data. We take a flat view of the data with no
        # striding
        R = as_strided(R, shape=(R.size, ), strides=(R.dtype.itemsize,))
        R = np.sum(R)
        return dictionary, R
    return dictionary


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
        V = _update_code_parallel(U, Y, alpha / n_samples, code=V,
                                  method=method, n_jobs=n_jobs)

        # Update U
        U, residuals = _update_dict(U, Y, V, verbose=verbose, return_r2=True)

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
