# Authors: Manoj Kumar
#          Thomas Unterthiner
#          Hamzeh Alsalhi

# License: BSD 3 clause
from __future__ import division
import scipy.sparse as sp
import numpy as np
import array

from numpy.random import permutation
from sklearn.utils.random import choice
from sklearn.utils.random import sample_without_replacement
from .fixes import sparse_min_max
from .sparsefuncs_fast import (csr_mean_variance_axis0,
                               csc_mean_variance_axis0)


def _raise_typeerror(X):
    """Raises a TypeError if X is not a CSR or CSC matrix"""
    input_type = X.format if sp.issparse(X) else type(X)
    err = "Expected a CSR or CSC sparse matrix, got %s." % input_type
    raise TypeError(err)


def inplace_csr_column_scale(X, scale):
    """Inplace column scaling of a CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale: float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    assert scale.shape[0] == X.shape[1]
    X.data *= scale.take(X.indices, mode='clip')


def inplace_csr_row_scale(X, scale):
    """ Inplace row scaling of a CSR matrix.

    Scale each sample of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
    matrix to be scaled.

    scale: float array with shape (n_samples,)
    Array of precomputed sample-wise values to use for scaling.
    """
    assert scale.shape[0] == X.shape[0]
    X.data *= np.repeat(scale, np.diff(X.indptr))


def mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSR or CSC matrix

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    if isinstance(X, sp.csr_matrix):
        return csr_mean_variance_axis0(X)
    elif isinstance(X, sp.csc_matrix):
        return csc_mean_variance_axis0(X)
    else:
        _raise_typeerror(X)


def inplace_column_scale(X, scale):
    """Inplace column scaling of a CSC/CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSC or CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale: float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    if isinstance(X, sp.csc_matrix):
        inplace_csr_row_scale(X.T, scale)
    elif isinstance(X, sp.csr_matrix):
        inplace_csr_column_scale(X, scale)
    else:
        _raise_typeerror(X)


def inplace_row_scale(X, scale):
    """ Inplace row scaling of a CSR or CSC matrix.

    Scale each row of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
    matrix to be scaled.

    scale: float array with shape (n_features,)
    Array of precomputed sample-wise values to use for scaling.
    """
    if isinstance(X, sp.csc_matrix):
        inplace_csr_column_scale(X.T, scale)
    elif isinstance(X, sp.csr_matrix):
        inplace_csr_row_scale(X, scale)
    else:
        _raise_typeerror(X)


def inplace_swap_row_csc(X, m, n):
    """
    Swaps two rows of a CSC matrix in-place.

    Parameters
    ----------
    X: scipy.sparse.csc_matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m: int
        Index of the row of X to be swapped.

    n: int
        Index of the row of X to be swapped.
    """
    for t in [m, n]:
        if isinstance(t, np.ndarray):
            raise TypeError("m and n should be valid integers")

    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    m_mask = X.indices == m
    X.indices[X.indices == n] = m
    X.indices[m_mask] = n


def inplace_swap_row_csr(X, m, n):
    """
    Swaps two rows of a CSR matrix in-place.

    Parameters
    ----------
    X: scipy.sparse.csr_matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m: int
        Index of the row of X to be swapped.

    n: int
        Index of the row of X to be swapped.
    """
    for t in [m, n]:
        if isinstance(t, np.ndarray):
            raise TypeError("m and n should be valid integers")

    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    # The following swapping makes life easier since m is assumed to be the
    # smaller integer below.
    if m > n:
        m, n = n, m

    indptr = X.indptr
    m_start = indptr[m]
    m_stop = indptr[m + 1]
    n_start = indptr[n]
    n_stop = indptr[n + 1]
    nz_m = m_stop - m_start
    nz_n = n_stop - n_start

    if nz_m != nz_n:
        # Modify indptr first
        X.indptr[m + 2:n] += nz_n - nz_m
        X.indptr[m + 1] = m_start + nz_n
        X.indptr[n] = n_stop - nz_m

    X.indices = np.concatenate([X.indices[:m_start],
                                X.indices[n_start:n_stop],
                                X.indices[m_stop:n_start],
                                X.indices[m_start:m_stop],
                                X.indices[n_stop:]])
    X.data = np.concatenate([X.data[:m_start],
                             X.data[n_start:n_stop],
                             X.data[m_stop:n_start],
                             X.data[m_start:m_stop],
                             X.data[n_stop:]])


def inplace_swap_row(X, m, n):
    """
    Swaps two rows of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m: int
        Index of the row of X to be swapped.

    n: int
        Index of the row of X to be swapped.
    """
    if isinstance(X, sp.csc_matrix):
        return inplace_swap_row_csc(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        return inplace_swap_row_csr(X, m, n)
    else:
        _raise_typeerror(X)


def inplace_swap_column(X, m, n):
    """
    Swaps two columns of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape=(n_samples, n_features)
        Matrix whose two columns are to be swapped.

    m: int
        Index of the column of X to be swapped.

    n : int
        Index of the column of X to be swapped.
    """
    if m < 0:
        m += X.shape[1]
    if n < 0:
        n += X.shape[1]
    if isinstance(X, sp.csc_matrix):
        return inplace_swap_row_csr(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        return inplace_swap_row_csc(X, m, n)
    else:
        _raise_typeerror(X)


def min_max_axis(X, axis):
    """Compute minimum and maximum along axis 0 on a CSR or CSC matrix

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    mins: float array with shape (n_features,)
        Feature-wise minima

    maxs: float array with shape (n_features,)
        Feature-wise maxima
    """
    if isinstance(X, sp.csr_matrix) or isinstance(X, sp.csc_matrix):
        return sparse_min_max(X, axis=axis)
    else:
        _raise_typeerror(X)


def random_choice_csc(n_samples, classes, class_probability=None,
                      random_state=None):
    """Generate a sparse random matrix given column class distributions

    Parameters
    ----------
    n_samples : int,
        Number of samples to draw in each column.

    classes : list of size n_outputs of arrays of size (n_classes, )
        List of classes for each column.

    class_probability : list of size n_outputs of arrays of size (n_classes, )
        Optional (default=None). Class distribution of each column. If None the
        uniform distribution is assumed.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Return
    ------
    random_matrix : sparse csc matrix of size (n_samples, n_outputs)

    """
    data = array.array('i')
    indices = array.array('i')
    indptr = array.array('i', [0])

    for j in range(len(classes)):
        classes[j] = classes[j].astype(int)

        # use uniform distribution if no class_probability is given
        if class_probability is None:
            class_prob_j = np.empty(shape=classes[j].shape[0])
            class_prob_j.fill(1 / classes[j].shape[0])
        else:
            class_prob_j = class_probability[j]
        class_prob_j = np.asarray(class_prob_j)

        if class_prob_j.shape[0] != classes[j].shape[0]:
            raise ValueError("classes[{0}] (length {1}) and "
                             "class_probability[{0}] (length {2}) have "
                             "different length.".format(j,
                                                        classes[j].shape[0],
                                                        class_prob_j.shape[0]))

        # If 0 is not present in the classes insert it with a probability 0.0
        if 0 not in classes[j]:
            classes[j] = np.insert(classes[j], 0, 0)
            class_prob_j = np.insert(class_prob_j, 0, 0.0)

        # If there are nonzero classes choose randomly using class_probability
        if classes[j].shape[0] > 1:
            p_nonzero = 1 - class_prob_j[classes == 0]
            nnz = int(n_samples * p_nonzero)
            ind_sample = sample_without_replacement(n_population=n_samples,
                                                    n_samples=nnz,
                                                    random_state=random_state)
            indices.extend(ind_sample)

            # Normalize probabilites for the nonzero elements
            classes_j_nonzero = classes[j] != 0
            class_probability_nz = class_prob_j[classes_j_nonzero]
            class_probability_nz_norm = (class_probability_nz /
                                         np.sum(class_probability_nz))
            data.extend(choice(classes[j][classes_j_nonzero],
                               size=nnz,
                               p=class_probability_nz_norm))
        indptr.append(len(indices))

    return sp.csc_matrix((data, indices, indptr),
                         (n_samples, len(classes)),
                         dtype=int)


def sparse_class_distribution(y, sample_weight=None):
    """Compute class priors from sparse multioutput-multiclass target data

    Parameters
    ----------
    y : sparse matrix of size (n_samples, n_outputs)
        The labels for each example.

   sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Return
    ------
    classes : list of size n_outputs of arrays of size (n_classes, )
        List of classes for each column.

    n_classes : list of integrs of size n_outputs
        Number of classes in each column

    class_prior : list of size n_outputs of arrays of size (n_classes, )
        Class distribution of each column.

    """
    classes_ = []
    n_classes = []
    class_prior_ = []

    y = y.tocsc()
    y.eliminate_zeros()
    y_nnz = np.diff(y.indptr)

    n_samples, n_outputs = y.shape

    for k in range(n_outputs):
        col_nonzero = y.indices[y.indptr[k]:y.indptr[k+1]]
        # separate sample weights for zero and non-zero elements
        
        nz_sample_weight = None
        z_sample_weight_sum = y.shape[0] - y_nnz[k]
        if sample_weight is not None:
            nz_sample_weight = sample_weight[col_nonzero]
            z_sample_weight_sum = (np.sum(sample_weight) -
                                  np.sum(nz_sample_weight))

        classes, y_k = np.unique(y.data[y.indptr[k]:y.indptr[k+1]],
                                 return_inverse=True)
        class_prior = np.bincount(y_k, weights=nz_sample_weight)
        if y_nnz[k] < y.shape[0]:
            classes = np.insert(classes, 0, 0)
            class_prior = np.insert(class_prior, 0,
                                    z_sample_weight_sum)
        classes_.append(classes)
        n_classes.append(classes.shape[0])
        class_prior_.append(class_prior / float(class_prior.sum()))

    return (classes_, n_classes, class_prior_)
