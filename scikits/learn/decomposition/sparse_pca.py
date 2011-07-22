""" Matrix factorization with Sparse PCA
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD

import time
import sys

from math import sqrt, floor, ceil
import itertools

import numpy as np
from numpy.lib.stride_tricks import as_strided
from scipy import linalg

from ..utils import check_random_state
from ..utils import gen_even_slices
from ..utils.extmath import fast_svd
from ..linear_model import Lasso, lars_path, ridge_regression
from ..externals.joblib import Parallel, delayed, cpu_count
from ..base import BaseEstimator, TransformerMixin


def _update_code(dictionary, Y, alpha, code=None, Gram=None, method='lars',
                 tol=1e-8):
    """ Update the sparse code factor in sparse_pca loop.
    Each column of the result is the solution to a Lasso problem.

    Parameters
    ----------
    dictionary: array of shape (n_samples, n_components)
        Dictionary against which to optimize the sparse code.

    Y: array of shape (n_samples, n_features)
        Data matrix.

    alpha: float
        Regularization parameter for the Lasso problem.

    code: array of shape (n_components, n_features)
        Value of the sparse codes at the previous iteration.

    Gram: array of shape (n_features, n_features)
        Precomputed Gram matrix, (Y^T * Y).

    method: {'lars', 'cd'}
        lars: uses the least angle regression method (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    tol: float
        Numerical tolerance for coordinate descent Lasso convergence.
        Only used if `method='cd'`

    Returns
    -------
    new_code : array of shape (n_components, n_features)
        The sparse codes precomputed using this iteration's dictionary
    """
    n_features = Y.shape[1]
    n_atoms = dictionary.shape[1]
    new_code = np.empty((n_atoms, n_features))
    if Gram is None:
        Gram = np.dot(dictionary.T, dictionary)
    if method == 'lars':
        XY = np.dot(dictionary.T, Y)
        try:
            err_mgt = np.seterr(all='ignore')
            for k in range(n_features):
                # A huge amount of time is spent in this loop. It needs to be
                # tight.
                _, _, coef_path_ = lars_path(dictionary, Y[:, k], Xy=XY[:, k],
                                             Gram=Gram, alpha_min=alpha,
                                             method='lasso')
                new_code[:, k] = coef_path_[:, -1]
        finally:
            np.seterr(**err_mgt)
    elif method == 'cd':
        clf = Lasso(alpha=alpha, fit_intercept=False, precompute=Gram)
        for k in range(n_features):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            if code is not None:
                clf.coef_ = code[:, k]  # Init with previous value of Vk
            clf.fit(dictionary, Y[:, k], max_iter=1000, tol=tol)
            new_code[:, k] = clf.coef_
    else:
        raise NotImplemented("Lasso method %s is not implemented." % method)
    return new_code


def _update_code_parallel(dictionary, Y, alpha, code=None, Gram=None,
                          method='lars', n_jobs=1, tol=1e-8):
    """ Update the sparse factor V in sparse_pca loop by efficiently
    spreading the load over the available cores.

    Parameters
    ----------
    dictionary: array of shape (n_samples, n_components)
        Dictionary against which to optimize the sparse code.

    Y: array of shape (n_samples, n_features)
        Data matrix.

    alpha: float
        Regularization parameter for the Lasso problem.

    code: array of shape (n_components, n_features)
        Previous iteration of the sparse code.

    Gram: array of shape (n_features, n_features)
        Precomputed Gram matrix, (Y^T * Y).

    method: 'lars' | 'cd'
        lars: uses the least angle regression method (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        lasso solution (linear_model.Lasso). Lars will be faster if
        the components extracted are sparse.

    n_jobs: int
        Number of parallel jobs to run.

    tol: float
        Numerical tolerance for coordinate descent Lasso convergence.
        Only used if `method='cd`.

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
    slices = list(gen_even_slices(n_features, n_jobs))
    code_views = Parallel(n_jobs=n_jobs)(
                delayed(_update_code)(dictionary, Y[:, this_slice],
                                      code=code[:, this_slice], alpha=alpha,
                                      Gram=Gram, method=method, tol=tol)
                for this_slice in slices)
    for this_slice, this_view in zip(slices, code_views):
        code[:, this_slice] = this_view
    return code


def _update_dict(dictionary, Y, code, verbose=False, return_r2=False,
                 random_state=None):
    """ Update the dense dictionary factor in place.

    Parameters
    ----------
    dictionary: array of shape (n_samples, n_components)
        Value of the dictionary at the previous iteration.

    Y: array of shape (n_samples, n_features)
        Data matrix.

    code: array of shape (n_components, n_features)
        Sparse coding of the data against which to optimize the dictionary.

    verbose:
        Degree of output the procedure will print.

    return_r2: bool
        Whether to compute and return the residual sum of squares corresponding
        to the computed solution.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Returns
    -------
    dictionary: array of shape (n_samples, n_components)
        Updated dictionary.

    """
    n_atoms = len(code)
    n_samples = Y.shape[0]
    random_state = check_random_state(random_state)
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
            dictionary[:, k] = random_state.randn(n_samples)
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


def dict_learning(X, n_atoms, alpha, max_iter=100, tol=1e-8, method='lars',
                  n_jobs=1, dict_init=None, code_init=None, callback=None,
                  verbose=False, random_state=None):
    """Solves a dictionary learning matrix factorization problem.

    Finds the best dictionary and the corresponding sparse code for
    approximating the data matrix X by solving::

    (U^*, V^*) = argmin 0.5 || X - U V ||_2^2 + alpha * || U ||_1
                 (U,V)
                with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    where V is the dictionary and U is the sparse code.

    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        Data matrix.

    n_atoms: int,
        Number of dictionary atoms to extract.

    alpha: int,
        Sparsity controlling parameter.

    max_iter: int,
        Maximum number of iterations to perform.

    tol: float,
        Tolerance for the stopping condition.

    method: {'lars', 'cd'}
        lars: uses the least angle regression method (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    n_jobs: int,
        Number of parallel jobs to run, or -1 to autodetect.

    dict_init: array of shape (n_atoms, n_features),
        Initial value for the dictionary for warm restart scenarios.

    code_init: array of shape (n_samples, n_atoms),
        Initial value for the sparse code for warm restart scenarios.

    callback:
        Callable that gets invoked every five iterations.

    verbose:
        Degree of output the procedure will print.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Returns
    -------
    code: array of shape (n_samples, n_atoms)
        The sparse code factor in the matrix factorization.

    dictionary: array of shape (n_atoms, n_features),
        The dictionary factor in the matrix factorization.

    errors: array
        Vector of errors at each iteration.

    """
    t0 = time.time()
    n_features = X.shape[1]
    # Avoid integer division problems
    alpha = float(alpha)
    random_state = check_random_state(random_state)

    if n_jobs == -1:
        n_jobs = cpu_count()

    # Init U and V with SVD of Y
    if code_init is not None and code_init is not None:
        code = np.array(code_init, order='F')
        # Don't copy V, it will happen below
        dictionary = dict_init
    else:
        code, S, dictionary = linalg.svd(X, full_matrices=False)
        dictionary = S[:, np.newaxis] * dictionary
    r = len(dictionary)
    if n_atoms <= r:
        code = code[:, :n_atoms]
        dictionary = dictionary[:n_atoms, :]
    else:
        code = np.c_[code, np.zeros((len(code), n_atoms - r))]
        dictionary = np.r_[dictionary,
                           np.zeros((n_atoms - r, dictionary.shape[1]))]

    # Fortran-order dict, as we are going to access its row vectors
    #code = np.array(code, order='F')
    dictionary = np.array(dictionary, order='F')

    residuals = 0

    errors = []
    current_cost = np.nan

    if verbose == 1:
        print '[dict_learning]',

    for ii in xrange(max_iter):
        dt = (time.time() - t0)
        if verbose == 1:
            sys.stdout.write(".")
            sys.stdout.flush()
        elif verbose:
            print ("Iteration % 3i "
                "(elapsed time: % 3is, % 4.1fmn, current cost % 7.3f)" %
                    (ii, dt, dt / 60, current_cost))

        # Update code
        code = _update_code_parallel(dictionary.T, X.T, alpha / n_features,
                                     code.T, method=method, n_jobs=n_jobs)
        code = code.T
        # Update dictionary
        dictionary, residuals = _update_dict(dictionary.T, X.T, code.T,
                                             verbose=verbose, return_r2=True,
                                             random_state=random_state)
        dictionary = dictionary.T

        # Cost function
        current_cost = 0.5 * residuals + alpha * np.sum(np.abs(code))
        errors.append(current_cost)

        if ii > 0:
            dE = errors[-2] - errors[-1]
            assert(dE >= -tol * errors[-1])
            if dE < tol * errors[-1]:
                if verbose == 1:
                    # A line return
                    print ""
                elif verbose:
                    print "--- Convergence reached after %d iterations" % ii
                break
        if ii % 5 == 0 and callback is not None:
            callback(locals())

    return code, dictionary, errors


def dict_learning_online(X, n_atoms, alpha, n_iter=100, return_code=True,
                         dict_init=None, callback=None, chunk_size=3,
                         verbose=False, shuffle=True, n_jobs=1,
                         method='lars', iter_offset=0, random_state=None):
    """Solves a dictionary learning matrix factorization problem online.

    Finds the best dictionary and the corresponding sparse code for
    approximating the data matrix X by solving:

    (U^*, V^*) = argmin 0.5 || X - U V ||_2^2 + alpha * || U ||_1
                 (U,V)
                 with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    where V is the dictionary and U is the sparse code. This is accomplished
    by repeatedly iterating over mini-batches of the input data.

    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        data matrix

    n_atoms: int,
        number of dictionary atoms to extract

    alpha: int,
        sparsity controlling parameter

    n_iter: int,
        number of iterations to perform

    return_code: boolean,
        whether to also return the code U or just the dictionary V

    dict_init: array of shape (n_atoms, n_features),
        initial value for the dictionary for warm restart scenarios

    callback:
        callable that gets invoked every five iterations

    chunk_size: int,
        the number of samples to take in each batch

    verbose:
        degree of output the procedure will print

    shuffle: boolean,
        whether to shuffle the data before splitting it in batches

    n_jobs: int,
        number of parallel jobs to run, or -1 to autodetect.

    method: 'lars' | 'cd'
        lars: uses the least angle regression method (linear_model.lars_path)
        cd: uses the stochastic gradient descent method to compute the
            lasso solution (linear_model.Lasso)

    iter_offset: int, default 0
        number of previous iterations completed on the dictionary used for
        initialization

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Returns
    -------
    dictionary: array of shape (n_atoms, n_features),
        the solutions to the dictionary learning problem

    code: array of shape (n_samples, n_atoms),
        the sparse code (only returned if `return_code=True`)
    """
    t0 = time.time()
    n_samples, n_features = X.shape
    # Avoid integer division problems
    alpha = float(alpha)
    random_state = check_random_state(random_state)

    if n_jobs == -1:
        n_jobs = cpu_count()

    # Init V with SVD of X
    if dict_init is not None:
        dictionary = dict_init
    else:
        _, S, dictionary = fast_svd(X, n_atoms)
        dictionary = S[:, np.newaxis] * dictionary
    r = len(dictionary)
    if n_atoms <= r:
        dictionary = dictionary[:n_atoms, :]
    else:
        dictionary = np.r_[dictionary,
                           np.zeros((n_atoms - r, dictionary.shape[1]))]
    dictionary = np.ascontiguousarray(dictionary.T)

    if verbose == 1:
        print '[dict_learning]',

    n_batches = floor(float(len(X)) / chunk_size)
    if shuffle:
        X_train = X.copy()
        random_state.shuffle(X_train)
    else:
        X_train = X
    batches = np.array_split(X_train, n_batches)
    batches = itertools.cycle(batches)

    # The covariance of the dictionary
    A = np.zeros((n_atoms, n_atoms))
    # The data approximation
    B = np.zeros((n_features, n_atoms))

    for ii, this_X in itertools.izip(xrange(iter_offset, iter_offset + n_iter),
                                     batches):
        #this_Y = this_Y.squeeze()
        dt = (time.time() - t0)
        if verbose == 1:
            sys.stdout.write(".")
            sys.stdout.flush()
        elif verbose:
            if verbose > 10 or ii % ceil(100. / verbose) == 0:
                print ("Iteration % 3i (elapsed time: % 3is, % 4.1fmn)" %
                    (ii, dt, dt / 60))

        this_code = _update_code(dictionary, this_X.T, alpha, method=method)

        # Update the auxiliary variables
        if ii < chunk_size - 1:
            theta = float((ii + 1) * chunk_size)
        else:
            theta = float(chunk_size ** 2 + ii + 1 - chunk_size)
        beta = (theta + 1 - chunk_size) / (theta + 1)

        A *= beta
        A += np.dot(this_code, this_code.T)
        B *= beta
        B += np.dot(this_X.T, this_code.T)

        # Update dictionary
        dictionary = _update_dict(dictionary, B, A, verbose=verbose,
                                  random_state=random_state)
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
        code = _update_code_parallel(dictionary, X.T, alpha, n_jobs=n_jobs,
                    method=method)
        if verbose > 1:
            dt = (time.time() - t0)
            print 'done (total time: % 3is, % 4.1fmn)' % (dt, dt / 60)
        return dictionary.T, code.T

    return dictionary.T


class SparsePCA(BaseEstimator, TransformerMixin):
    """Sparse Principal Components Analysis (SparsePCA)

    Finds the set of sparse components that can optimally reconstruct the data.
    The amount of sparseness is controllable by the coefficient of the L1
    penalty, given by the parameter alpha.

    Parameters
    ----------
    n_components: int,
        Number of sparse atoms to extract.

    alpha: int,
        Sparsity controlling parameter.

    max_iter: int,
        Maximum number of iterations to perform.

    tol: float,
        Tolerance for the stopping condition.

    method: {'lars', 'cd'}
        lars: uses the least angle regression method (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    n_jobs: int,
        Number of parallel jobs to run.

    U_init: array of shape (n_samples, n_atoms),
        Initial values for the loadings for warm restart scenarios.

    V_init: array of shape (n_atoms, n_features),
        Initial values for the components for warm restart scenarios.

    verbose:
        Degree of verbosity of the printed output.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    components_: array, [n_components, n_features]
        Sparse components extracted from the data.

    error_: array
        Vector of errors at each iteration.

    See also
    --------
    PCA

    """
    def __init__(self, n_components=None, alpha=1, max_iter=1000, tol=1e-8,
                 method='lars', n_jobs=1, U_init=None, V_init=None,
                 verbose=False, random_state=None):
        self.n_components = n_components
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.n_jobs = n_jobs
        self.U_init = U_init
        self.V_init = V_init
        self.verbose = verbose
        self.random_state = random_state

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
        self._set_params(**params)
        self.random_state = check_random_state(self.random_state)
        X = np.asanyarray(X)
        code_init = self.V_init.T if self.V_init != None else None
        dict_init = self.U_init.T if self.U_init != None else None
        Vt, _, E = dict_learning(X.T, self.n_components, self.alpha,
                                tol=self.tol, max_iter=self.max_iter,
                                method=self.method, n_jobs=self.n_jobs,
                                verbose=self.verbose,
                                random_state=self.random_state,
                                code_init=code_init,
                                dict_init=dict_init)
        self.components_ = Vt.T
        self.error_ = E
        return self

    def transform(self, X, ridge_alpha=0.01):
        """Least Squares projection of the data onto the learned sparse
        components.

        To avoid instability issues in case the system is under-determined,
        regularization can be applied (Ridge regression) via the
        `ridge_alpha` parameter.

        Note that Sparse PCA components orthogonality is not enforced as in PCA
        hence one cannot use a simple linear projection.

        Parameters
        ----------
        X: array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        ridge_alpha: float, default: 0.01
            Amount of ridge shrinkage to apply in order to improve
            conditioning.

        Returns
        -------
        X_new array, shape (n_samples, n_components)
            Transformed data.
        """
        U = ridge_regression(self.components_.T, X.T, ridge_alpha,
                             solver='dense_cholesky')
        U /= np.sqrt((U ** 2).sum(axis=0))
        return U
