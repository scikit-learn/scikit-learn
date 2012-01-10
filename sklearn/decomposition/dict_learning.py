""" Dictionary learning
"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD

import time
import sys
import itertools
import warnings

from math import sqrt, floor, ceil

import numpy as np
from scipy import linalg
from numpy.lib.stride_tricks import as_strided

from ..base import BaseEstimator, TransformerMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..utils import array2d, check_random_state, gen_even_slices, deprecated
from ..utils.extmath import randomized_svd
from ..linear_model import Lasso, orthogonal_mp_gram, lars_path


def _sparse_encode(X, dictionary, gram=None, cov=None, algorithm='lasso_lars',
                  n_nonzero_coefs=None, alpha=None, copy_gram=True,
                  copy_cov=True, init=None, max_iter=1000):
    """Generic sparse coding

    Each column of the result is the solution to a Lasso problem.

    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        Data matrix.

    dictionary: array of shape (n_atoms, n_features)
        The dictionary matrix against which to solve the sparse coding of
        the data. Some of the algorithms assume normalized rows.

    gram: array, shape=(n_atoms, n_atoms)
        Precomputed Gram matrix, dictionary * dictionary'

    cov: array, shape=(n_atoms, n_samples)
        Precomputed covariance, dictionary * X'

    algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection dictionary * data'

    n_nonzero_coefs: int, 0.1 * n_features by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    alpha: float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threshold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    init: array of shape (n_samples, n_atoms)
        Initialization value of the sparse code. Only used if
        `algorithm='lasso_cd'`.

    max_iter: int, 1000 by default
        Maximum number of iterations to perform if `algorithm='lasso_cd'`.

    copy_gram: boolean, optional
        Whether to copy the precomputed Gram matrix; if False, it may be
        overwritten.

    copy_cov: boolean, optional
        Whether to copy the precomputed covariance matrix; if False, it may be
        overwritten.

    Returns
    -------
    code: array of shape (n_components, n_features)
        The sparse codes

    See also
    --------
    sklearn.linear_model.lars_path
    sklearn.linear_model.orthogonal_mp
    sklearn.linear_model.Lasso
    SparseCoder
    """
    alpha = float(alpha) if alpha is not None else None
    dictionary = np.asarray(dictionary)
    X = np.asarray(X)
    if X.ndim == 1:
        X = X[:, np.newaxis]
    n_samples, n_features = X.shape
    n_atoms = dictionary.shape[0]
    # This will always use Gram
    if gram is None:
        # I think it's never safe to overwrite Gram when n_features > 1
        # but I'd like to avoid the complicated logic.
        # The parameter could be removed in this case. Discuss.
        gram = np.dot(dictionary, dictionary.T)
    if cov is None and algorithm != 'lasso_cd':
        # overwriting cov is safe
        copy_cov = False
        cov = np.dot(dictionary, X.T)

    if algorithm == 'lasso_lars':
        if alpha is None:
            alpha = 1.
        alpha /= n_features  # account for scaling
        try:
            new_code = np.empty((n_samples, n_atoms))
            err_mgt = np.seterr(all='ignore')
            for k in range(n_samples):
                # A huge amount of time is spent in this loop. It needs to be
                # tight.
                _, _, coef_path_ = lars_path(dictionary.T, X[k], Xy=cov[:, k],
                                             Gram=gram, alpha_min=alpha,
                                             method='lasso')
                new_code[k] = coef_path_[:, -1]
        finally:
            np.seterr(**err_mgt)

    elif algorithm == 'lasso_cd':
        if alpha is None:
            alpha = 1.
        alpha /= n_features  # account for scaling
        new_code = np.empty((n_samples, n_atoms))
        clf = Lasso(alpha=alpha, fit_intercept=False, precompute=gram,
                    max_iter=1000)
        for k in xrange(n_samples):
            # A huge amount of time is spent in this loop. It needs to be
            # tight

            if init is not None:
                clf.coef_ = init[k]  # Init with previous value of the code
            clf.fit(dictionary.T, X[k])
            new_code[k] = clf.coef_

    elif algorithm == 'lars':
        if n_nonzero_coefs is None:
            n_nonzero_coefs = max(n_features / 10, 1)
        try:
            new_code = np.empty((n_samples, n_atoms))
            err_mgt = np.seterr(all='ignore')
            for k in xrange(n_samples):
                # A huge amount of time is spent in this loop. It needs to be
                # tight.
                _, _, coef_path_ = lars_path(dictionary.T, X[k], Xy=cov[:, k],
                                             Gram=gram, method='lar',
                                             max_iter=n_nonzero_coefs)
                new_code[k] = coef_path_[:, -1]
        finally:
            np.seterr(**err_mgt)

    elif algorithm == 'threshold':
        if alpha is None:
            alpha = 1.
        new_code = (np.sign(cov) * np.maximum(np.abs(cov) - alpha, 0)).T

    elif algorithm == 'omp':
        if n_nonzero_coefs is None and alpha is None:
            n_nonzero_coefs = max(n_features / 10, 1)
        norms_squared = np.sum((X ** 2), axis=1)
        new_code = orthogonal_mp_gram(gram, cov, n_nonzero_coefs, alpha,
                                      norms_squared, copy_Xy=copy_cov
                                      ).T
    else:
        raise NotImplemented('Sparse coding method %s not implemented' %
                             algorithm)
    return new_code


def sparse_encode(X, dictionary, gram=None, cov=None, algorithm='lasso_lars',
                  n_nonzero_coefs=None, alpha=None, copy_gram=True,
                  copy_cov=True, init=None, max_iter=1000, n_jobs=1):
    """Sparse coding

    Each row of the result is the solution to a sparse coding problem.
    The goal is to find a sparse array `code` such that::

        X ~= code * dictionary

    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        Data matrix

    dictionary: array of shape (n_atoms, n_features)
        The dictionary matrix against which to solve the sparse coding of
        the data. Some of the algorithms assume normalized rows for meaningful
        output.

    gram: array, shape=(n_atoms, n_atoms)
        Precomputed Gram matrix, dictionary * dictionary'

    cov: array, shape=(n_atoms, n_samples)
        Precomputed covariance, dictionary' * X

    algorithm: {'lasso_lars', 'lasso_cd', 'lars', 'omp', 'threshold'}
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection dictionary * X'

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

    init: array of shape (n_samples, n_atoms)
        Initialization value of the sparse codes. Only used if
        `algorithm='lasso_cd'`.

    max_iter: int, 1000 by default
        Maximum number of iterations to perform if `algorithm='lasso_cd'`.

    copy_gram: boolean, optional
        Whether to copy the precomputed Gram matrix; if False, it may be
        overwritten.

    copy_cov: boolean, optional
        Whether to copy the precomputed covariance matrix; if False, it may be
        overwritten.

    n_jobs: int, optional
        Number of parallel jobs to run.

    Returns
    -------
    code: array of shape (n_samples, n_atoms)
        The sparse codes

    See also
    --------
    sklearn.linear_model.lars_path
    sklearn.linear_model.orthogonal_mp
    sklearn.linear_model.Lasso
    SparseCoder
    """
    warnings.warn("Please note: the interface of sparse_encode has changed: "
                  "It now follows the dictionary learning API and it also "
                  "handles parallelization. Please read the docstring for "
                  "more information.")
    dictionary = np.asarray(dictionary)
    X = np.asarray(X)
    n_samples, n_features = X.shape
    n_atoms = dictionary.shape[0]
    if gram is None:
        copy_gram = False
        gram = np.dot(dictionary, dictionary.T)
    if cov is None and algorithm != 'lasso_cd':
        copy_cov = False
        cov = np.dot(dictionary, X.T)
    if n_jobs == 1 or algorithm == 'threshold':
        return _sparse_encode(X, dictionary, gram, cov, algorithm,
                             n_nonzero_coefs, alpha, copy_gram, copy_cov, init)
    code = np.empty((n_samples, n_atoms))
    slices = list(gen_even_slices(n_samples, n_jobs))
    code_views = Parallel(n_jobs=n_jobs)(
                delayed(sparse_encode)(X[this_slice], dictionary, gram,
                                       cov[:, this_slice], algorithm,
                                       n_nonzero_coefs, alpha,
                                       copy_gram, copy_cov,
                                       init=init[this_slice] if init is not
                                       None else None)
                for this_slice in slices)
    for this_slice, this_view in zip(slices, code_views):
        code[this_slice] = this_view
    return code


@deprecated('Use sparse_encode instead')
def sparse_encode_parallel():
    pass


def _update_dict(dictionary, Y, code, verbose=False, return_r2=False,
                 random_state=None):
    """Update the dense dictionary factor in place.

    Parameters
    ----------
    dictionary: array of shape (n_features, n_atoms)
        Value of the dictionary at the previous iteration.

    Y: array of shape (n_features, n_samples)
        Data matrix.

    code: array of shape (n_atoms, n_samples)
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
    dictionary: array of shape (n_features, n_atoms)
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

    See also
    --------
    dict_learning_online
    DictionaryLearning
    MiniBatchDictionaryLearning
    SparsePCA
    MiniBatchSparsePCA
    """
    if method not in ('lars', 'cd'):
        raise ValueError('Coding method not supported as a fit algorithm.')
    method = 'lasso_' + method

    t0 = time.time()
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
        code = sparse_encode(X, dictionary, algorithm=method, alpha=alpha,
                             init=code, n_jobs=n_jobs)
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
    approximating the data matrix X by solving::

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
    code: array of shape (n_samples, n_atoms),
        the sparse code (only returned if `return_code=True`)

    dictionary: array of shape (n_atoms, n_features),
        the solutions to the dictionary learning problem

    See also
    --------
    dict_learning
    DictionaryLearning
    MiniBatchDictionaryLearning
    SparsePCA
    MiniBatchSparsePCA

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
        _, S, dictionary = randomized_svd(X, n_atoms)
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

        this_code = sparse_encode(this_X, dictionary.T, algorithm=method,
                                  alpha=alpha).T

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
        code = sparse_encode(X, dictionary.T, algorithm=method, alpha=alpha,
                             n_jobs=n_jobs)
        if verbose > 1:
            dt = (time.time() - t0)
            print 'done (total time: % 3is, % 4.1fmn)' % (dt, dt / 60)
        return code, dictionary.T

    return dictionary.T


class SparseCodingMixin(TransformerMixin):
    """Sparse coding mixin"""

    def _set_sparse_coding_params(self, n_atoms, transform_algorithm='omp',
                                  transform_n_nonzero_coefs=None,
                                  transform_alpha=None, split_sign=False,
                                  n_jobs=1):
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
        X : array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Transformed data

        """
        # XXX : kwargs is not documented
        X = array2d(X)
        n_samples, n_features = X.shape

        code = sparse_encode(
            X, self.components_, algorithm=self.transform_algorithm,
            n_nonzero_coefs=self.transform_n_nonzero_coefs,
            alpha=self.transform_alpha, n_jobs=self.n_jobs)

        if self.split_sign:
            # feature vector is split into a positive and negative side
            n_samples, n_features = code.shape
            split_code = np.empty((n_samples, 2 * n_features))
            split_code[:, :n_features] = np.maximum(code, 0)
            split_code[:, n_features:] = -np.minimum(code, 0)
            code = split_code

        return code


class SparseCoder(BaseEstimator, SparseCodingMixin):
    """Sparse coding

    Finds a sparse representation of data against a fixed, precomputed
    dictionary.

    Each row of the result is the solution to a sparse coding problem.
    The goal is to find a sparse array `code` such that::

        X ~= code * dictionary

    Parameters
    ----------
    dictionary : array, [n_atoms, n_features]
        The dictionary atoms used for sparse coding. Lines are assumed to be
        normalized to unit norm.

    transform_algorithm : {'lasso_lars', 'lasso_cd', 'lars', 'omp', \
    'threshold'}
        Algorithm used to transform the data:
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection ``dictionary * X'``

    transform_n_nonzero_coefs : int, ``0.1 * n_features`` by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    transform_alpha : float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threshold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    split_sign : bool, False by default
        Whether to split the sparse feature vector into the concatenation of
        its negative part and its positive part. This can improve the
        performance of downstream classifiers.

    n_jobs : int,
        number of parallel jobs to run

    Attributes
    ----------
    `components_` : array, [n_atoms, n_features]
        The unchanged dictionary atoms

    See also
    --------
    DictionaryLearning
    MiniBatchDictionaryLearning
    SparsePCA
    MiniBatchSparsePCA
    sparse_encode
    """

    def __init__(self, dictionary, transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 split_sign=False, n_jobs=1):
        self._set_sparse_coding_params(dictionary.shape[0],
                                       transform_algorithm,
                                       transform_n_nonzero_coefs,
                                       transform_alpha, split_sign, n_jobs)
        self.components_ = dictionary

    def fit(self, X, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence
        work in pipelines.
        """
        return self


class DictionaryLearning(BaseEstimator, SparseCodingMixin):
    """Dictionary learning

    Finds a dictionary (a set of atoms) that can best be used to represent data
    using a sparse code.

    Solves the optimization problem::

        (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
                    (U,V)
                    with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    Parameters
    ----------
    n_atoms : int,
        number of dictionary elements to extract

    alpha : int,
        sparsity controlling parameter

    max_iter : int,
        maximum number of iterations to perform

    tol : float,
        tolerance for numerical error

    fit_algorithm : {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    transform_algorithm : {'lasso_lars', 'lasso_cd', 'lars', 'omp', \
    'threshold'}
        Algorithm used to transform the data
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection ``dictionary * X'``

    transform_n_nonzero_coefs : int, ``0.1 * n_features`` by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    transform_alpha : float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threshold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    split_sign : bool, False by default
        Whether to split the sparse feature vector into the concatenation of
        its negative part and its positive part. This can improve the
        performance of downstream classifiers.

    n_jobs : int,
        number of parallel jobs to run

    code_init : array of shape (n_samples, n_atoms),
        initial value for the code, for warm restart

    dict_init : array of shape (n_atoms, n_features),
        initial values for the dictionary, for warm restart

    verbose :
        degree of verbosity of the printed output

    random_state : int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    `components_` : array, [n_atoms, n_features]
        dictionary atoms extracted from the data

    `error_` : array
        vector of errors at each iteration

    Notes
    -----
    **References:**

    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)

    See also
    --------
    SparseCoder
    MiniBatchDictionaryLearning
    SparsePCA
    MiniBatchSparsePCA
    """
    def __init__(self, n_atoms, alpha=1, max_iter=1000, tol=1e-8,
                 fit_algorithm='lars', transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 n_jobs=1, code_init=None, dict_init=None, verbose=False,
                 split_sign=False, random_state=None):
        self._set_sparse_coding_params(n_atoms, transform_algorithm,
                                       transform_n_nonzero_coefs,
                                       transform_alpha, split_sign, n_jobs)
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
        X = np.asarray(X)
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


class MiniBatchDictionaryLearning(BaseEstimator, SparseCodingMixin):
    """Mini-batch dictionary learning

    Finds a dictionary (a set of atoms) that can best be used to represent data
    using a sparse code.

    Solves the optimization problem::

       (U^*,V^*) = argmin 0.5 || Y - U V ||_2^2 + alpha * || U ||_1
                    (U,V)
                    with || V_k ||_2 = 1 for all  0 <= k < n_atoms

    Parameters
    ----------
    n_atoms : int,
        number of dictionary elements to extract

    alpha : int,
        sparsity controlling parameter

    n_iter : int,
        total number of iterations to perform

    fit_algorithm : {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    transform_algorithm : {'lasso_lars', 'lasso_cd', 'lars', 'omp', \
    'threshold'}
        Algorithm used to transform the data.
        lars: uses the least angle regression method (linear_model.lars_path)
        lasso_lars: uses Lars to compute the Lasso solution
        lasso_cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). lasso_lars will be faster if
        the estimated components are sparse.
        omp: uses orthogonal matching pursuit to estimate the sparse solution
        threshold: squashes to zero all coefficients less than alpha from
        the projection dictionary * X'

    transform_n_nonzero_coefs : int, ``0.1 * n_features`` by default
        Number of nonzero coefficients to target in each column of the
        solution. This is only used by `algorithm='lars'` and `algorithm='omp'`
        and is overridden by `alpha` in the `omp` case.

    transform_alpha : float, 1. by default
        If `algorithm='lasso_lars'` or `algorithm='lasso_cd'`, `alpha` is the
        penalty applied to the L1 norm.
        If `algorithm='threshold'`, `alpha` is the absolute value of the
        threshold below which coefficients will be squashed to zero.
        If `algorithm='omp'`, `alpha` is the tolerance parameter: the value of
        the reconstruction error targeted. In this case, it overrides
        `n_nonzero_coefs`.

    split_sign : bool, False by default
        Whether to split the sparse feature vector into the concatenation of
        its negative part and its positive part. This can improve the
        performance of downstream classifiers.

    n_jobs : int,
        number of parallel jobs to run

    dict_init : array of shape (n_atoms, n_features),
        initial value of the dictionary for warm restart scenarios

    verbose :
        degree of verbosity of the printed output

    chunk_size : int,
        number of samples in each mini-batch

    shuffle : bool,
        whether to shuffle the samples before forming batches

    random_state : int or RandomState
        Pseudo number generator state used for random sampling.

    Attributes
    ----------
    `components_` : array, [n_atoms, n_features]
        components extracted from the data

    Notes
    -----
    **References:**

    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
    for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)

    See also
    --------
    SparseCoder
    DictionaryLearning
    SparsePCA
    MiniBatchSparsePCA

    """
    def __init__(self, n_atoms, alpha=1, n_iter=1000,
                 fit_algorithm='lars', n_jobs=1, chunk_size=3,
                 shuffle=True, dict_init=None, transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 verbose=False, split_sign=False, random_state=None):

        self._set_sparse_coding_params(n_atoms, transform_algorithm,
                                       transform_n_nonzero_coefs,
                                       transform_alpha, split_sign, n_jobs)
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
        X = np.asarray(X)
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
        X = array2d(X)
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
