""" Dictionary learning
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD

import time
import sys
import itertools

from math import sqrt, floor, ceil

import numpy as np
from scipy import linalg
from numpy.lib.stride_tricks import as_strided

from ..base import BaseEstimator, TransformerMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..utils import check_random_state
from ..utils import gen_even_slices
from ..utils.extmath import fast_svd
from ..linear_model import Lasso, orthogonal_mp_gram, lars_path


def sparse_encode(X, Y, gram=None, cov=None, algorithm='lasso_lars',
                  n_nonzero_coefs=None, alpha=None,
                  overwrite_gram=False, overwrite_cov=False, init=None):
    """Generic sparse coding

    Each column of the result is the solution to a Lasso problem.

    Parameters
    ----------
    X: array of shape (n_samples, n_components)
        Dictionary against which to optimize the sparse code.

    Y: array of shape (n_samples, n_features)
        Data matrix.

    gram: array, shape=(n_components, n_components)
        Precomputed Gram matrix, X^T * X

    cov: array, shape=(n_components, n_features)
        Precomputed covariance, X^T * Y

    algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection X.T * Y

    n_nonzero_coefs: int, 0.1 * n_features by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    alpha: float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threhold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    init: array of shape (n_components, n_features)
        Initialization value of the sparse codes. Only used if
        `algorithm='lasso_cd'`.

    overwrite_gram: boolean,
        Whether to overwrite the precomputed Gram matrix.

    overwrite_cov: boolean,
        Whether to overwrite the precomputed covariance matrix.

    Returns
    -------
    code: array of shape (n_components, n_features)
        The sparse codes

    See also
    --------
    linear_model.lars_path
    linear_model.orthogonal_mp
    linear_model.Lasso
    """
    alpha = float(alpha) if alpha is not None else None
    X, Y = map(np.asanyarray, (X, Y))
    if Y.ndim == 1:
        Y = Y[:, np.newaxis]
    n_features = Y.shape[1]
    # This will always use Gram
    if gram is None:
        # I think it's never safe to overwrite Gram when n_features > 1
        # but I'd like to avoid the complicated logic.
        # The parameter could be removed in this case. Discuss.
        gram = np.dot(X.T, X)
    if cov is None and algorithm != 'lasso_cd':
        # overwrite_cov is safe
        overwrite_cov = True
        cov = np.dot(X.T, Y)

    if algorithm == 'lasso_lars':
        if alpha is None:
            alpha = 1.
        try:
            new_code = np.empty((X.shape[1], n_features))
            err_mgt = np.seterr(all='ignore')
            for k in range(n_features):
                # A huge amount of time is spent in this loop. It needs to be
                # tight.
                _, _, coef_path_ = lars_path(X, Y[:, k], Xy=cov[:, k],
                                             Gram=gram, alpha_min=alpha,
                                             method='lasso')
                new_code[:, k] = coef_path_[:, -1]
        finally:
            np.seterr(**err_mgt)

    elif algorithm == 'lasso_cd':
        if alpha is None:
            alpha = 1.
        new_code = np.empty((X.shape[1], n_features))
        clf = Lasso(alpha=alpha, fit_intercept=False, precompute=gram,
                    max_iter=1000)
        for k in xrange(n_features):
            # A huge amount of time is spent in this loop. It needs to be
            # tight.
            if init is not None:
                clf.coef_ = init[:, k]  # Init with previous value of Vk
            clf.fit(X, Y[:, k])
            new_code[:, k] = clf.coef_

    elif algorithm == 'lars':
        if n_nonzero_coefs is None:
            n_nonzero_coefs = n_features / 10
        try:
            new_code = np.empty((X.shape[1], n_features))
            err_mgt = np.seterr(all='ignore')
            for k in xrange(n_features):
                # A huge amount of time is spent in this loop. It needs to be
                # tight.
                _, _, coef_path_ = lars_path(X, Y[:, k], Xy=cov[:, k],
                                             Gram=gram, method='lar',
                                             max_iter=n_nonzero_coefs)
                new_code[:, k] = coef_path_[:, -1]
        finally:
            np.seterr(**err_mgt)

    elif algorithm == 'threshold':
        if alpha is None:
            alpha = 1.
        new_code = np.sign(cov) * np.maximum(np.abs(cov) - alpha, 0)

    elif algorithm == 'omp':
        if n_nonzero_coefs is None and alpha is None:
            n_nonzero_coefs = n_features / 10
        norms_squared = np.sum((Y ** 2), axis=0)
        new_code = orthogonal_mp_gram(gram, cov, n_nonzero_coefs, alpha,
                                      norms_squared, overwrite_Xy=overwrite_cov
                                      )
    else:
        raise NotImplemented('Sparse coding method %s not implemented' %
                             algorithm)
    return new_code


def sparse_encode_parallel(X, Y, gram=None, cov=None, algorithm='lasso_lars',
                  n_nonzero_coefs=None, alpha=None, overwrite_gram=False,
                  overwrite_cov=False, init=None, n_jobs=1):
    """Parallel sparse coding using joblib

    Each column of the result is the solution to a Lasso problem.

    Parameters
    ----------
    X: array of shape (n_samples, n_components)
        Dictionary against which to optimize the sparse code.

    Y: array of shape (n_samples, n_features)
        Data matrix.

    gram: array, shape=(n_components, n_components)
        Precomputed Gram matrix, X^T * X

    cov: array, shape=(n_components, n_features)
        Precomputed covariance, X^T * Y

    algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection X.T * Y

    n_nonzero_coefs: int, 0.1 * n_features by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    alpha: float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threhold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    init: array of shape (n_components, n_features)
        Initialization value of the sparse codes. Only used if
        `algorithm='lasso_cd'`.

    overwrite_gram: boolean,
        Whether to overwrite the precomputed Gram matrix.

    overwrite_cov: boolean,
        Whether to overwrite the precomputed covariance matrix.

    n_jobs: int,
        Number of parallel jobs to run.

    Returns
    -------
    code: array of shape (n_components, n_features)
        The sparse codes

    See also
    --------
    linear_model.lars_path
    linear_model.orthogonal_mp
    linear_model.Lasso
    """
    n_samples, n_features = Y.shape
    n_components = X.shape[1]
    if gram is None:
        overwrite_gram = True
        gram = np.dot(X.T, X)
    if cov is None and algorithm != 'lasso_cd':
        overwrite_cov = True
        cov = np.dot(X.T, Y)
    if n_jobs == 1 or algorithm == 'threshold':
        return sparse_encode(X, Y, gram, cov, algorithm, n_nonzero_coefs,
                             alpha, overwrite_gram, overwrite_cov, init)
    code = np.empty((n_components, n_features))
    slices = list(gen_even_slices(n_features, n_jobs))
    code_views = Parallel(n_jobs=n_jobs)(
                delayed(sparse_encode)(X, Y[:, this_slice], gram,
                                       cov[:, this_slice], algorithm,
                                       n_nonzero_coefs, alpha,
                                       overwrite_gram, overwrite_cov,
                                       init=init[:, this_slice] if init is not
                                       None else None)
                for this_slice in slices)
    for this_slice, this_view in zip(slices, code_views):
        code[:, this_slice] = this_view
    return code


def _update_dict(dictionary, Y, code, verbose=False, return_r2=False,
                 random_state=None):
    """Update the dense dictionary factor in place.

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


def dict_learning(X, n_atoms, alpha, max_iter=100, tol=1e-8,
                  method='lars', n_jobs=1, dict_init=None, code_init=None,
                  callback=None, verbose=False, random_state=None):
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
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
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
    if method not in ('lars', 'cd'):
        raise ValueError('Coding method not supported as a fit algorithm.')
    method = 'lasso_' + method

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
    if n_atoms <= r:  # True even if n_atoms=None
        code = code[:, :n_atoms]
        dictionary = dictionary[:n_atoms, :]
    else:
        code = np.c_[code, np.zeros((len(code), n_atoms - r))]
        dictionary = np.r_[dictionary,
                           np.zeros((n_atoms - r, dictionary.shape[1]))]

    # Fortran-order dict, as we are going to access its row vectors
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
        code = sparse_encode_parallel(dictionary.T, X.T, algorithm=method,
                                      alpha=alpha / n_features,
                                      init=code.T, n_jobs=n_jobs)
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
            # assert(dE >= -tol * errors[-1])
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

    where V is the dictionary and U is the sparse code. This is
    accomplished by repeatedly iterating over mini-batches by slicing
    the input data.

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

    method: {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

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
    if method not in ('lars', 'cd'):
        raise ValueError('Coding method not supported as a fit algorithm.')
    method = 'lasso_' + method

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
        dt = (time.time() - t0)
        if verbose == 1:
            sys.stdout.write(".")
            sys.stdout.flush()
        elif verbose:
            if verbose > 10 or ii % ceil(100. / verbose) == 0:
                print ("Iteration % 3i (elapsed time: % 3is, % 4.1fmn)" %
                    (ii, dt, dt / 60))

        this_code = sparse_encode(dictionary, this_X.T, algorithm=method,
                                  alpha=alpha)

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
        code = sparse_encode_parallel(dictionary, X.T, algorithm=method,
                                      alpha=alpha, n_jobs=n_jobs)
        if verbose > 1:
            dt = (time.time() - t0)
            print 'done (total time: % 3is, % 4.1fmn)' % (dt, dt / 60)
        return code.T, dictionary.T

    return dictionary.T


class BaseDictionaryLearning(BaseEstimator, TransformerMixin):
    """Dictionary learning base class"""

    def __init__(self, n_atoms, transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 split_sign=False, n_jobs=1):
        self.n_atoms = n_atoms
        self.transform_algorithm = transform_algorithm
        self.transform_n_nonzero_coefs = transform_n_nonzero_coefs
        self.transform_alpha = transform_alpha
        self.split_sign = split_sign
        self.n_jobs = n_jobs

    def transform(self, X, y=None):
        """Encode the data as a sparse combination of the dictionary atoms.

        Coding method is determined by the object parameter
        `transform_algorithm`.

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
        # XXX : kwargs is not documented
        X = np.atleast_2d(X)
        n_samples, n_features = X.shape

        code = sparse_encode_parallel(
            self.components_.T, X.T, algorithm=self.transform_algorithm,
            n_nonzero_coefs=self.transform_n_nonzero_coefs,
            alpha=self.transform_alpha, n_jobs=self.n_jobs)
        code = code.T

        if self.split_sign:
            # feature vector is split into a positive and negative side
            n_samples, n_features = code.shape
            split_code = np.empty((n_samples, 2 * n_features))
            split_code[:, :n_features] = np.maximum(code, 0)
            split_code[:, n_features:] = -np.minimum(code, 0)
            code = split_code

        return code


class DictionaryLearning(BaseDictionaryLearning):
    """ Dictionary learning

    Finds a dictionary (a set of atoms) that can best be used to represent data
    using a sparse code.

    Solves the optimization problem:
    (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
                 (U,V)
                with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    Parameters
    ----------
    n_atoms: int,
        number of dictionary elements to extract

    alpha: int,
        sparsity controlling parameter

    max_iter: int,
        maximum number of iterations to perform

    tol: float,
        tolerance for numerical error

    fit_algorithm: {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    transform_algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        Algorithm used to transform the data
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection X.T * Y

    transform_n_nonzero_coefs: int, 0.1 * n_features by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    transform_alpha: float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threhold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    n_jobs: int,
        number of parallel jobs to run

    code_init: array of shape (n_samples, n_atoms),
        initial value for the code, for warm restart

    dict_init: array of shape (n_atoms, n_features),
        initial values for the dictionary, for warm restart

    verbose:
        degree of verbosity of the printed output

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    components_: array, [n_atoms, n_features]
        dictionary atoms extracted from the data

    error_: array
        vector of errors at each iteration

    References
    ----------
    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)


    See also
    --------
    :class:`sklearn.decomposition.SparsePCA` which solves the transposed
    problem, finding sparse components to represent data.

    """
    def __init__(self, n_atoms, alpha=1, max_iter=1000, tol=1e-8,
                 fit_algorithm='lars', transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 n_jobs=1, code_init=None, dict_init=None, verbose=False,
                 split_sign=False, random_state=None):
        BaseDictionaryLearning.__init__(self, n_atoms, transform_algorithm,
                 transform_n_nonzero_coefs, transform_alpha, split_sign,
                 n_jobs)
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.fit_algorithm = fit_algorithm
        self.code_init = code_init
        self.dict_init = dict_init
        self.verbose = verbose
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self: object
            Returns the object itself
        """
        self.random_state = check_random_state(self.random_state)
        X = np.asanyarray(X)
        V, U, E = dict_learning(X, self.n_atoms, self.alpha,
                                tol=self.tol, max_iter=self.max_iter,
                                method=self.fit_algorithm,
                                n_jobs=self.n_jobs,
                                code_init=self.code_init,
                                dict_init=self.dict_init,
                                verbose=self.verbose,
                                random_state=self.random_state)
        self.components_ = U
        self.error_ = E
        return self


class MiniBatchDictionaryLearning(BaseDictionaryLearning):
    """Mini-batch dictionary learning

    Finds a dictionary (a set of atoms) that can best be used to represent data
    using a sparse code.

    Solves the optimization problem:
    (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
                 (U,V)
                with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    Parameters
    ----------
    n_atoms: int,
        number of dictionary elements to extract

    alpha: int,
        sparsity controlling parameter

    n_iter: int,
        total number of iterations to perform

    fit_algorithm: {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    transform_algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        Algorithm used to transform the data.
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection X.T * Y

    transform_n_nonzero_coefs: int, 0.1 * n_features by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    transform_alpha: float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threshold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    n_jobs: int,
        number of parallel jobs to run

    dict_init: array of shape (n_atoms, n_features),
        initial value of the dictionary for warm restart scenarios

    verbose:
        degree of verbosity of the printed output

    chunk_size: int,
        number of samples in each mini-batch

    shuffle: bool,
        whether to shuffle the samples before forming batches

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    components_: array, [n_atoms, n_features]
        components extracted from the data

    References
    ----------
    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)


    See also
    --------
    :class:`sklearn.decomposition.SparsePCA` which solves the transposed
    problem, finding sparse components to represent data.

    """
    def __init__(self, n_atoms, alpha=1, n_iter=1000,
                 fit_algorithm='lars', n_jobs=1, chunk_size=3,
                 shuffle=True, dict_init=None, transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 verbose=False, split_sign=False, random_state=None):
        BaseDictionaryLearning.__init__(self, n_atoms, transform_algorithm,
                 transform_n_nonzero_coefs, transform_alpha, split_sign,
                 n_jobs)
        self.alpha = alpha
        self.n_iter = n_iter
        self.fit_algorithm = fit_algorithm
        self.dict_init = dict_init
        self.verbose = verbose
        self.shuffle = shuffle
        self.chunk_size = chunk_size
        self.split_sign = split_sign
        self.random_state = random_state

    def fit(self, X, y=None):
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
        self.random_state = check_random_state(self.random_state)
        X = np.asanyarray(X)
        U = dict_learning_online(X, self.n_atoms, self.alpha,
                                 n_iter=self.n_iter, return_code=False,
                                 method=self.fit_algorithm,
                                 n_jobs=self.n_jobs,
                                 dict_init=self.dict_init,
                                 chunk_size=self.chunk_size,
                                 shuffle=self.shuffle, verbose=self.verbose,
                                 random_state=self.random_state)
        self.components_ = U
        return self

    def partial_fit(self, X, y=None, iter_offset=0):
        """Updates the model using the data in X as a mini-batch.

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
        self.random_state = check_random_state(self.random_state)
        X = np.atleast_2d(X)
        if hasattr(self, 'components_'):
            dict_init = self.components_
        else:
            dict_init = self.dict_init
        U = dict_learning_online(X, self.n_atoms, self.alpha,
                                 n_iter=self.n_iter,
                                 method=self.fit_algorithm,
                                 n_jobs=self.n_jobs, dict_init=dict_init,
                                 chunk_size=len(X), shuffle=False,
                                 verbose=self.verbose, return_code=False,
                                 iter_offset=iter_offset,
                                 random_state=self.random_state)
        self.components_ = U
        return self
