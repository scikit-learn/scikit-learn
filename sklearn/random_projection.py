"""Random Projection transformers

Random Projections are an efficient way to reduce the dimensionality
of the data by trading a controlled amout of accuracy (as additional
variance) for faster processing times and smaller model sizes.

The dimensions and distribution of Random Projections matrices are
controlled so as to preserve the pairwise distances between any two
samples of the dataset.

The main theoretical result behind their efficiency is the
Johnson-Lindenstrauss lemma (quoting Wikipedia):

  In mathematics, the Johnson-Lindenstrauss lemma is a result
  concerning low-distortion embeddings of points from high-dimensional
  into low-dimensional Euclidean space. The lemma states that a small set
  of points in a high-dimensional space can be embedded into a space of
  much lower dimension in such a way that distances between the points are
  nearly preserved. The map used for the embedding is at least Lipschitz,
  and can even be taken to be an orthogonal projection.

  http://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma

"""
# Authors: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simple BSD

import math
import warnings

import numpy as np
import scipy.sparse as sp

from sklearn.utils import check_random_state
from sklearn.base import BaseEstimator
from sklearn.base import TransformerMixin


def johnson_lindenstrauss_bound(n_samples, eps=0.1):
    """Find a 'safe' number of components to randomly project to

    Minimum dimensionality of a random projection eps-embedding such
    that with good probability:

      (1 - eps) ||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps) ||u - v||^2

    Where u and v are any rows taken from a dataset of shape [n_samples,
    n_features] and p is a projection by a random gaussian N(0, 1) matrix
    with shape [n_components, n_features] (or a sparse Achlioptas matrix).

    The minimum number of components to garantee the embedding is given by:

      n_components >= 4 log(n_samples) / (eps^2 / 2 - eps^3 / 3)

    Note that the number of dimensions is independent of the original
    number of features but instead depends on the size of the dataset:
    the larger the dataset, the higher is the minimal dimensionality of
    random projection embedding.

    Examples
    --------

    >>> johnson_lindenstrauss_bound(1e6, eps=0.5)
    663

    >>> johnson_lindenstrauss_bound(1e6, eps=0.1)
    11841

    >>> johnson_lindenstrauss_bound(1e6, eps=0.01)
    1112658

    References
    ----------
    - http://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma

    - An elementary proof of the Johnson-Lindenstrauss Lemma.
      Sanjoy Dasgupta and Anupam Gupta, 1999
      http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.45.3654

    """
    denominator = (eps ** 2 / 2) - (eps ** 3 / 3)
    return int(4 * math.log(n_samples) / denominator)


def sparse_random_matrix(n_components, n_features, density='auto',
                         random_state=None):
    """Generalized Achlioptas random sparse matrix for random projection

    Setting density to 1/3 will yield the original matrix by Dimitris
    Achlioptas while setting a lower value will yield the generalization
    by Ping Li et al:

    If we note `s = 1 / density` the components of the random matrix are:

      - -sqrt(s) / sqrt(n_components)   with probability 1 / 2s
      -  0                              with probability 1 - 1 / s
      - +sqrt(s) / sqrt(n_components)   with probability 1 / 2s

    Parameters
    ----------
    n_components: int
        Dimensionality of the target projection space.

    density: float in range (0, 1/3], optional
        Ratio of non-zero component in the random projection matrix.

        By default the value is set to the minimum density as recommended
        by Ping Li et al.: 1 / sqrt(n_features)

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

    Examples
    --------

      >>> import numpy as np
      >>> from sklearn.random_projection import sparse_random_matrix

      >>> n_components, n_features = 10, 10000

      >>> r = sparse_random_matrix(n_components, n_features, random_state=0)
      >>> r                                   # doctest: +NORMALIZE_WHITESPACE
      <10x10000 sparse matrix of type '<type 'numpy.float64'>'
          with 1002 stored elements in Compressed Sparse Row format>

    The random matrix has only two possible non-zero values::

      >>> np.unique(r.data)                              # doctest: +ELLIPSIS
      array([-3.16...,  3.16...])

    The density is adjusted based on the number of features::

      >>> expected_density = 1 / math.sqrt(n_features)
      >>> actual_density = float(r.nnz) / (n_components * n_features)
      >>> np.abs(actual_density - expected_density) < 0.0001
      True

    The matrix is centered on zero::

      >>> np.abs(r.mean())                                # doctest: +ELLIPSIS
      0.00...

    References
    ----------

    - Very Sparse Random Projections. Ping Li, Trevor Hastie
      and Kenneth Church, 2006
      http://www.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf

    - Database-friendly random projections, Dimitris Achlioptas, 2001
      http://www.cs.ucsc.edu/~optas/papers/jl.pdf

    """
    random_state = check_random_state(random_state)

    if density is 'auto':
        density = min(1 / math.sqrt(n_features), 1 / 3.)
    elif density <= 0 or density > 1 / float(3):
        raise ValueError("Expected density in range (0, 1/3], got: %r"
                         % density)

    # placeholders for the CSR datastructure
    indices = []
    data = []
    offset = 0
    indptr = [offset]

    prob_nonzero = density
    for i in xrange(n_components):
        # find the indices of the non-zero components for row i
        u = random_state.uniform(size=n_features)
        indices_i = np.arange(n_features)[u < prob_nonzero].copy()
        indices.append(indices_i)

        # among non zero components the probability of the sign is 50%/50%
        n_nonzero_i = indices_i.shape[0]
        data_i = np.ones(n_nonzero_i)
        u = random_state.uniform(size=n_nonzero_i)
        data_i[u < 0.5] *= -1
        data.append(data_i)
        offset += n_nonzero_i
        indptr.append(offset)

    # build the CSR structure by concatenating the rows
    r = sp.csr_matrix(
        (np.concatenate(data), np.concatenate(indices), np.array(indptr)),
        shape=(n_components, n_features))

    return math.sqrt(1 / density) / math.sqrt(n_components) * r


class SparseRandomProjection(BaseEstimator, TransformerMixin):
    """Transformer to reduce the dimensionality with sparse random projection

    Alternative to the dense Gaussian Random matrix that garantees
    similar embedding quality while being much more memory efficient
    and allowing faster computation of the projected data.

    The implementation uses a CSR matrix internally.

    If we note `s = 1 / density` the components of the random matrix are:

      - -sqrt(s) / sqrt(n_components)   with probability 1 / 2s
      -  0                              with probability 1 - 1 / s
      - +sqrt(s) / sqrt(n_components)   with probability 1 / 2s

    Parameters
    ----------
    n_components: int, optional
        Dimensionality of the target projection space.

        By default n_components is automatically adjusted according to
        the number of samples in the dataset and the bound given by the
        Johnson Lindenstrauss lemma.

    density: float in range (0, 1/3], optional
        Ratio of non-zero component in the random projection matrix.

        By default the value is set to the minimum density as recommended
        by Ping Li et al.: 1 / sqrt(n_features)

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

    eps: strictly positive float, optional, default 0.1
        Parameter to control the quality of the embedding according to
        the Johnson-Lindenstrauss lemma when n_components is set to
        'auto'.

        Smaller values lead to better embedding and higher number of
        dimensions (n_components) in the target projection space.

    random_state: integer, RandomState instance or None (default)
        Control the pseudo random number generator used to generate the
        matrix at fit time.

    Attributes
    ----------
    components_: CSR matrix with shape [n_components, n_features]
        Random matrix used for the projection.

    References
    ----------

    - Very Sparse Random Projections. Ping Li, Trevor Hastie
      and Kenneth Church, 2006
      http://www.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf

    - Database-friendly random projections, Dimitris Achlioptas, 2001
      http://www.cs.ucsc.edu/~optas/papers/jl.pdf

    """

    def __init__(self, n_components='auto', density='auto', eps=0.1,
                 random_state=None):
        self.n_components = n_components
        self.density = density
        self.eps = eps
        self.random_state = random_state

    def fit(self, X, y=None):
        """Generate a sparse random projection matrix

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            Training set: only the shape is used to find optimal random
            matrix dimensions based on the theory referenced in the
            afore mentioned papers.

        y : is not used: placehold to allow for usage in a Pipeline.

        Returns
        -------
        The fitted estimator.

        """
        if not sp.issparse(X):
            X = np.atleast_2d(X)
        self.random_state = check_random_state(self.random_state)
        n_samples, n_features = X.shape

        if self.n_components == 'auto':
            self.n_components = johnson_lindenstrauss_bound(
                n_samples, eps=self.eps)
            if self.n_components > n_features:
                warnings.warn(
                    'eps=%f and n_samples=%d lead to a target dimension of '
                    '%d which is larger than the original space with '
                    'n_features=%d' % (self.eps, n_samples, self.n_components,
                                       n_features))
                self.n_components = n_features

        if self.density == 'auto':
            self.density = min(1 / math.sqrt(n_features), 1 / 3.)

        self.components_ = sparse_random_matrix(
            self.n_components, n_features, density=self.density,
            random_state=self.random_state)
        return self

    def transform(self, X, y=None):
        """Project the data by using matrix product with the random matrix

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            The input data to project into a smaller dimensional space.

        y : is not used: placehold to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.

        """
        if not sp.issparse(X):
            X = np.atleast_2d(X)
        return X * self.components_.T
