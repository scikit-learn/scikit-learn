"""K-SVD dictionary learning algorithm.
"""


import numpy
import scipy.sparse.linalg

from .dict_learning import sparse_encode
from ..preprocessing import normalize
from ..utils import check_array, check_random_state
from ..utils.extmath import norm


# atoms with norm below are considered zero
_ATOM_NORM_TOLERANCE = 1e-6


def _worst_represented_example(X, dictionary, sparse_codes):
    """Find the sample that has the most representation error."""

    # TODO: [optimize] can use squared norm or extmath.norm
    errors = numpy.linalg.norm(X - numpy.dot(sparse_codes, dictionary), axis=1)
    worst_index = numpy.argmax(errors)
    return X[worst_index, :]


def _decompose_svd(residual):
    """Decompose the residual using 1-rank SVD."""

    # Sparse version is used because non-sparse does not allow to
    # set the desired number of singular values.
    U, S, V = scipy.sparse.linalg.svds(residual, 1)
    # U is unitary, so d is normalized
    d = U.ravel()
    g = numpy.dot(V.T, S).ravel()
    return d, g


def _decompose_projections(residual, g_old):
    """Decompose the residual using SVD-approximating projections."""

    d = numpy.dot(residual, g_old)
    d_norm = norm(d)
    if d_norm >= _ATOM_NORM_TOLERANCE:
        d /= d_norm
    g = numpy.dot(residual.T, d)
    return d, g


def _atom_update(X, atom_index, dictionary, sparse_codes, approximate_svd):
    """Update single dictionary atom."""

    atom_usages = numpy.nonzero(sparse_codes[:, atom_index])[0]
    if len(atom_usages) == 0:
        atom = _worst_represented_example(X, dictionary, sparse_codes)
        dictionary[atom_index, :] = atom / norm(atom)
    else:
        dictionary[atom_index, :] = 0
        # TODO: [optimize] can save memory by not-calculating this, or doing it inplace
        representation = numpy.dot(sparse_codes[atom_usages, :], dictionary)
        residual = (X[atom_usages, :] - representation).T
        if approximate_svd:
            atom, atom_codes = _decompose_projections(
                residual, sparse_codes[atom_usages, atom_index])
        else:
            atom, atom_codes = _decompose_svd(residual)
        dictionary[atom_index, :] = atom
        sparse_codes[atom_usages, atom_index] = atom_codes.T

    # TODO: [optimize] update and return R
    return dictionary, sparse_codes


def ksvd(X, n_nonzero_coefs, n_atoms, iteration_count, init_dictionary=None,
         approximate_svd=True, random_state=None, n_jobs=1):
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

    n_atoms : int
        Number of atoms (columns) in the resulting dictionary.

    iteration_count : int
        Maximum number of learning iterations.

    init_dictionary : array of shape (n_atoms, n_features)
        Initional dictionary for learning.  If the value is `None`,
        the random training samples are chosen as initional dictionary
        columns.

    approximate_svd : boolean, True by default
        If True, the SVD-decomposition is approximated for faster
        execution.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    n_jobs: int, 1 by default
        Number of parallel jobs to run.

    Returns
    -------
    sparse_codes : array of shape (n_samples, n_atoms)
        Sparse codes for all the training samples.

    dictionary : array of shape (n_atoms, n_features)
        The learned dictionary.

    Notes
    -----
    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

    """

    X = check_array(X)
    n_samples, n_features = X.shape

    random_state = check_random_state(random_state)
    if init_dictionary is None:
        # TODO: [optimize] can use squared norm or extmath.norm
        nonzero_examples = list(numpy.where(
            numpy.linalg.norm(X, axis=1) > _ATOM_NORM_TOLERANCE)[0])
        chosen_examples = random_state.choice(nonzero_examples, replace=False,
                                              size=n_atoms)
        dictionary = X[chosen_examples, :]
    else:
        assert init_dictionary.shape == (n_atoms, n_features), \
            "Initional dictionary must have consistent shape."
        dictionary = init_dictionary.copy()

    dictionary = normalize(dictionary)
    sparse_codes = numpy.zeros((n_samples, n_atoms))

    for iteration in range(iteration_count):
        sparse_codes = sparse_encode(X, dictionary, algorithm='omp',
                                     n_nonzero_coefs=n_nonzero_coefs,
                                     n_jobs=n_jobs)
        for atom_index in range(n_atoms):
            dictionary, sparse_codes = _atom_update(
                X, atom_index, dictionary, sparse_codes, approximate_svd)

    return sparse_codes, dictionary
