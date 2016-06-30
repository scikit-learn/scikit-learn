"""K-SVD dictionary learning algorithm.
"""


import numpy as np
import scipy.sparse.linalg

from ..base import BaseEstimator
from .dict_learning import sparse_encode, SparseCodingMixin
from ..preprocessing import normalize
from ..utils import check_array, check_random_state
from ..utils.extmath import norm
from ..utils.random import choice


# atoms with norm below are considered zero
_ATOM_NORM_TOLERANCE = 1e-6


def _worst_represented_example(X, dictionary, sparse_codes):
    """Find the sample that has the most representation error."""

    residuals = X - np.dot(sparse_codes, dictionary)
    errors_squared = np.sum(residuals * residuals, axis=1)
    worst_index = np.argmax(errors_squared)
    return X[worst_index, :]


def _decompose_svd(residual):
    """Decompose the residual using 1-rank SVD."""

    # Sparse version is used because non-sparse does not allow to
    # set the desired number of singular values.
    U, S, V = scipy.sparse.linalg.svds(residual, 1)
    # U is unitary, so d is normalized
    d = U.ravel()
    g = np.dot(V.T, S).ravel()
    return d, g


def _decompose_projections(residual, g_old):
    """Decompose the residual using SVD-approximating projections."""

    d = np.dot(residual, g_old)
    d_norm = norm(d)
    if d_norm >= _ATOM_NORM_TOLERANCE:
        d /= d_norm
    g = np.dot(residual.T, d)
    return d, g


def _atom_update(X, atom_index, dictionary, sparse_codes, approximate_svd):
    """Update single dictionary atom and the corresponding codes."""

    atom_usages = np.nonzero(sparse_codes[:, atom_index])[0]
    if len(atom_usages) == 0:
        # replace with the new atom
        atom = _worst_represented_example(X, dictionary, sparse_codes)
        dictionary[atom_index, :] = atom / norm(atom)
    elif len(atom_usages) == 1:
        # update is trivial if the atom is used in only one sample
        atom = X[atom_usages[0], :]
        atom_norm = norm(atom)
        dictionary[atom_index, :] = atom / atom_norm
        sparse_codes[atom_usages[0], atom_index] = atom_norm
    else:
        dictionary[atom_index, :] = 0
        representation = np.dot(sparse_codes[atom_usages, :], dictionary)
        residual = (X[atom_usages, :] - representation).T
        if approximate_svd:
            atom, atom_codes = _decompose_projections(
                residual, sparse_codes[atom_usages, atom_index])
        else:
            atom, atom_codes = _decompose_svd(residual)
        dictionary[atom_index, :] = atom
        sparse_codes[atom_usages, atom_index] = atom_codes.T

    return dictionary, sparse_codes


def _reconstruction_error(X, dictionary, sparse_codes):
    """Calculate the total reconstruction error for data."""
    errors = X - np.dot(sparse_codes, dictionary)
    return norm(errors)


def _init_dictionary_from_samples(X, n_components, random_state):
    samples_norms_squared = np.sum(X*X, axis=1)
    nonzero_examples = list(
        np.where(samples_norms_squared > _ATOM_NORM_TOLERANCE)[0])

    if len(nonzero_examples) < n_components:
        # this is not effective, but still allowed
        replace = True
    else:
        replace = False
    chosen_examples = choice(nonzero_examples, replace=replace,
                             size=n_components, random_state=random_state)
    return X[chosen_examples, :]


def learn_dictionary_ksvd(X, n_nonzero_coefs, n_components, iteration_count,
                          init_dictionary=None, approximate_svd=True,
                          random_state=None, n_jobs=1, callback=None):
    """K-SVD dictionary learning algorithm.

    Finds the best dictionary D and the corresponding sparse code C for
    the following problem:

        (C^*, D^*) = argmin || X - C D ||_2^2  s.t.
            ||C_i||_0 <= n_nonzero_coefs

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        Matrix with training samples as columns.

    n_nonzero_coefs : int
        Maximum number of atoms in sparse code for each example.

    n_components : int
        Number of atoms (columns) in the resulting dictionary.

    iteration_count : int
        Maximum number of learning iterations.

    init_dictionary : array of shape (n_components, n_features)
        Initional dictionary for learning.  If the value is `None`,
        the random training samples are chosen as initional dictionary
        columns.

    approximate_svd : boolean, True by default
        If True, the SVD-decomposition is approximated for faster
        execution.

    random_state : int or RandomState
        Pseudo number generator state used for random sampling.

    n_jobs : int, 1 by default
        Number of parallel jobs to run.

    callback : function
        The function that is called at each iteration.  The value of
        `locals()` is passed to the function.

    Returns
    -------
    sparse_codes : array of shape (n_samples, n_components)
        Sparse codes for all the training samples.

    dictionary : array of shape (n_components, n_features)
        The learned dictionary.

    errors : list
        List of errors after each iteration.

    Notes
    -----

    """

    X = check_array(X)
    n_samples, n_features = X.shape

    random_state = check_random_state(random_state)
    if init_dictionary is None:
        dictionary = _init_dictionary_from_samples(X, n_components,
                                                   random_state)
    else:
        assert init_dictionary.shape == (n_components, n_features), \
            "Initional dictionary must have consistent shape."
        dictionary = init_dictionary.copy()

    dictionary = normalize(dictionary)
    sparse_codes = np.zeros((n_samples, n_components))

    errors = []
    for iteration in range(iteration_count):
        sparse_codes = sparse_encode(X, dictionary, algorithm='omp',
                                     n_nonzero_coefs=n_nonzero_coefs,
                                     n_jobs=n_jobs)
        for atom_index in range(n_components):
            dictionary, sparse_codes = _atom_update(
                X, atom_index, dictionary, sparse_codes, approximate_svd)

        errors.append(_reconstruction_error(X, dictionary, sparse_codes))

        if callback is not None:
            callback(locals())

    return sparse_codes, dictionary, errors


def _setup_default_sparse_parameters(X, n_components, fit_n_nonzero_coefs):
    if n_components is None:
        n_components = X.shape[1]
    if fit_n_nonzero_coefs is None:
        fit_n_nonzero_coefs = max(int(0.1 * n_components), 1)

    return n_components, fit_n_nonzero_coefs


class DictionaryLearningKSVD(BaseEstimator, SparseCodingMixin):
    """Dictionary learning using K-SVD method.

    Finds the best dictionary D and the corresponding sparse code C for
    the following problem:

        (C^*, D^*) = argmin || X - C D ||_2^2  s.t.
            ||C_i||_0 <= n_nonzero_coefs

    Parameters
    ----------

    fit_n_nonzero_coefs : int
        Maximum number of atoms in sparse code for each example when
        learning.  If None (by default) this value is set to 10% of
        `n_components`.

    n_components : int
        Number of atoms (columns) in the resulting dictionary.  If None
        (by default), the number of features in learning data is used.

    iteration_count : int, 100 by default
        Maximum number of learning iterations.

    init_dictionary : array of shape (n_components, n_features)
        Initional dictionary for learning.  If the value is None,
        the random training samples are chosen as initional dictionary
        columns.

    approximate_svd : boolean, True by default
        If True, the SVD-decomposition is approximated for faster
        execution.

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

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    n_jobs: int, 1 by default
        Number of parallel jobs to run.

    callback : function
        The function that is called at each iteration.  The value of
        `locals()` is passed to the function.

    Attributes
    ----------
    components_ : array, [n_components, n_features]
        Dictionary atoms extracted from the data.

    error_ : list
        List of errors after each iteration.

    Notes
    -----
    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

    """

    def __init__(self, fit_n_nonzero_coefs=None, n_components=None,
                 iteration_count=100, init_dictionary=None,
                 approximate_svd=True, transform_algorithm='omp',
                 transform_n_nonzero_coefs=None, transform_alpha=None,
                 split_sign=False, random_state=None, n_jobs=1, callback=None):
        self._set_sparse_coding_params(n_components, transform_algorithm,
                                       transform_n_nonzero_coefs,
                                       transform_alpha, split_sign, n_jobs)
        self.fit_n_nonzero_coefs = fit_n_nonzero_coefs
        self.iteration_count = iteration_count
        self.init_dictionary = init_dictionary
        self.approximate_svd = approximate_svd
        self.random_state = random_state
        self.callback = callback

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X : array of shape (n_samples, n_features)
            Matrix with training samples as columns.

        Returns
        -------
        self: object
            Returns the object itself

        """

        random_state = check_random_state(self.random_state)
        X = check_array(X)

        n_components, fit_n_nonzero_coefs = _setup_default_sparse_parameters(
            X, self.n_components, self.fit_n_nonzero_coefs)

        codes, dictionary, errors = learn_dictionary_ksvd(
            X, n_nonzero_coefs=fit_n_nonzero_coefs,
            n_components=n_components,
            iteration_count=self.iteration_count,
            init_dictionary=self.init_dictionary,
            approximate_svd=self.approximate_svd, random_state=random_state,
            n_jobs=self.n_jobs, callback=self.callback)
        self.components_ = dictionary
        self.error_ = errors

        return self
