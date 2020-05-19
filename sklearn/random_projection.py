# -*- coding: utf8
"""Random Projection transformers

Random Projections are a simple and computationally efficient way to
reduce the dimensionality of the data by trading a controlled amount
of accuracy (as additional variance) for faster processing times and
smaller model sizes.

The dimensions and distribution of Random Projections matrices are
controlled so as to preserve the pairwise distances between any two
samples of the dataset.

The main theoretical result behind the efficiency of random projection is the
`Johnson-Lindenstrauss lemma (quoting Wikipedia)
<https://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma>`_:

  In mathematics, the Johnson-Lindenstrauss lemma is a result
  concerning low-distortion embeddings of points from high-dimensional
  into low-dimensional Euclidean space. The lemma states that a small set
  of points in a high-dimensional space can be embedded into a space of
  much lower dimension in such a way that distances between the points are
  nearly preserved. The map used for the embedding is at least Lipschitz,
  and can even be taken to be an orthogonal projection.

"""
# Authors: Olivier Grisel <olivier.grisel@ensta.org>,
#          Arnaud Joly <a.joly@ulg.ac.be>
# License: BSD 3 clause

import warnings
from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp

from .base import BaseEstimator, TransformerMixin

from .utils import check_random_state
from .utils.extmath import safe_sparse_dot
from .utils.random import sample_without_replacement
from .utils.validation import check_array, check_is_fitted
from .utils.validation import _deprecate_positional_args
from .exceptions import DataDimensionalityWarning
from .utils import deprecated


__all__ = ["SparseRandomProjection",
           "GaussianRandomProjection",
           "johnson_lindenstrauss_min_dim"]


@_deprecate_positional_args
def johnson_lindenstrauss_min_dim(n_samples, *, eps=0.1):
    """Find a 'safe' number of components to randomly project to

    The distortion introduced by a random projection `p` only changes the
    distance between two points by a factor (1 +- eps) in an euclidean space
    with good probability. The projection `p` is an eps-embedding as defined
    by:

      (1 - eps) ||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps) ||u - v||^2

    Where u and v are any rows taken from a dataset of shape [n_samples,
    n_features], eps is in ]0, 1[ and p is a projection by a random Gaussian
    N(0, 1) matrix with shape [n_components, n_features] (or a sparse
    Achlioptas matrix).

    The minimum number of components to guarantee the eps-embedding is
    given by:

      n_components >= 4 log(n_samples) / (eps^2 / 2 - eps^3 / 3)

    Note that the number of dimensions is independent of the original
    number of features but instead depends on the size of the dataset:
    the larger the dataset, the higher is the minimal dimensionality of
    an eps-embedding.

    Read more in the :ref:`User Guide <johnson_lindenstrauss>`.

    Parameters
    ----------
    n_samples : int or numpy array of int greater than 0,
        Number of samples. If an array is given, it will compute
        a safe number of components array-wise.

    eps : float or numpy array of float in ]0,1[, optional (default=0.1)
        Maximum distortion rate as defined by the Johnson-Lindenstrauss lemma.
        If an array is given, it will compute a safe number of components
        array-wise.

    Returns
    -------
    n_components : int or numpy array of int,
        The minimal number of components to guarantee with good probability
        an eps-embedding with n_samples.

    Examples
    --------

    >>> johnson_lindenstrauss_min_dim(1e6, eps=0.5)
    663

    >>> johnson_lindenstrauss_min_dim(1e6, eps=[0.5, 0.1, 0.01])
    array([    663,   11841, 1112658])

    >>> johnson_lindenstrauss_min_dim([1e4, 1e5, 1e6], eps=0.1)
    array([ 7894,  9868, 11841])

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma

    .. [2] Sanjoy Dasgupta and Anupam Gupta, 1999,
           "An elementary proof of the Johnson-Lindenstrauss Lemma."
           http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.45.3654

    """
    eps = np.asarray(eps)
    n_samples = np.asarray(n_samples)

    if np.any(eps <= 0.0) or np.any(eps >= 1):
        raise ValueError(
            "The JL bound is defined for eps in ]0, 1[, got %r" % eps)

    if np.any(n_samples) <= 0:
        raise ValueError(
            "The JL bound is defined for n_samples greater than zero, got %r"
            % n_samples)

    denominator = (eps ** 2 / 2) - (eps ** 3 / 3)
    return (4 * np.log(n_samples) / denominator).astype(np.int)


def _check_density(density, n_features):
    """Factorize density check according to Li et al."""
    if density == 'auto':
        density = 1 / np.sqrt(n_features)

    elif density <= 0 or density > 1:
        raise ValueError("Expected density in range ]0, 1], got: %r"
                         % density)
    return density


def _check_input_size(n_components, n_features):
    """Factorize argument checking for random matrix generation"""
    if n_components <= 0:
        raise ValueError("n_components must be strictly positive, got %d" %
                         n_components)
    if n_features <= 0:
        raise ValueError("n_features must be strictly positive, got %d" %
                         n_features)


# TODO: remove in 0.24
@deprecated("gaussian_random_matrix is deprecated in "
            "0.22 and will be removed in version 0.24.")
def gaussian_random_matrix(n_components, n_features, random_state=None):
    return _gaussian_random_matrix(n_components, n_features, random_state)


def _gaussian_random_matrix(n_components, n_features, random_state=None):
    """Generate a dense Gaussian random matrix.

    The components of the random matrix are drawn from

        N(0, 1.0 / n_components).

    Read more in the :ref:`User Guide <gaussian_random_matrix>`.

    Parameters
    ----------
    n_components : int,
        Dimensionality of the target projection space.

    n_features : int,
        Dimensionality of the original source space.

    random_state : int, RandomState instance or None, optional (default=None)
        Controls the pseudo random number generator used to generate the matrix
        at fit time.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    components : numpy array of shape [n_components, n_features]
        The generated Gaussian random matrix.

    See Also
    --------
    GaussianRandomProjection
    """
    _check_input_size(n_components, n_features)
    rng = check_random_state(random_state)
    components = rng.normal(loc=0.0,
                            scale=1.0 / np.sqrt(n_components),
                            size=(n_components, n_features))
    return components


# TODO: remove in 0.24
@deprecated("gaussian_random_matrix is deprecated in "
            "0.22 and will be removed in version 0.24.")
def sparse_random_matrix(n_components, n_features, density='auto',
                         random_state=None):
    return _sparse_random_matrix(n_components, n_features, density,
                                 random_state)


def _sparse_random_matrix(n_components, n_features, density='auto',
                          random_state=None):
    """Generalized Achlioptas random sparse matrix for random projection

    Setting density to 1 / 3 will yield the original matrix by Dimitris
    Achlioptas while setting a lower value will yield the generalization
    by Ping Li et al.

    If we note :math:`s = 1 / density`, the components of the random matrix are
    drawn from:

      - -sqrt(s) / sqrt(n_components)   with probability 1 / 2s
      -  0                              with probability 1 - 1 / s
      - +sqrt(s) / sqrt(n_components)   with probability 1 / 2s

    Read more in the :ref:`User Guide <sparse_random_matrix>`.

    Parameters
    ----------
    n_components : int,
        Dimensionality of the target projection space.

    n_features : int,
        Dimensionality of the original source space.

    density : float in range ]0, 1] or 'auto', optional (default='auto')
        Ratio of non-zero component in the random projection matrix.

        If density = 'auto', the value is set to the minimum density
        as recommended by Ping Li et al.: 1 / sqrt(n_features).

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

    random_state : int, RandomState instance or None, optional (default=None)
        Controls the pseudo random number generator used to generate the matrix
        at fit time.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    components : array or CSR matrix with shape [n_components, n_features]
        The generated Gaussian random matrix.

    See Also
    --------
    SparseRandomProjection

    References
    ----------

    .. [1] Ping Li, T. Hastie and K. W. Church, 2006,
           "Very Sparse Random Projections".
           https://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf

    .. [2] D. Achlioptas, 2001, "Database-friendly random projections",
           http://www.cs.ucsc.edu/~optas/papers/jl.pdf

    """
    _check_input_size(n_components, n_features)
    density = _check_density(density, n_features)
    rng = check_random_state(random_state)

    if density == 1:
        # skip index generation if totally dense
        components = rng.binomial(1, 0.5, (n_components, n_features)) * 2 - 1
        return 1 / np.sqrt(n_components) * components

    else:
        # Generate location of non zero elements
        indices = []
        offset = 0
        indptr = [offset]
        for _ in range(n_components):
            # find the indices of the non-zero components for row i
            n_nonzero_i = rng.binomial(n_features, density)
            indices_i = sample_without_replacement(n_features, n_nonzero_i,
                                                   random_state=rng)
            indices.append(indices_i)
            offset += n_nonzero_i
            indptr.append(offset)

        indices = np.concatenate(indices)

        # Among non zero components the probability of the sign is 50%/50%
        data = rng.binomial(1, 0.5, size=np.size(indices)) * 2 - 1

        # build the CSR structure by concatenating the rows
        components = sp.csr_matrix((data, indices, indptr),
                                   shape=(n_components, n_features))

        return np.sqrt(1 / density) / np.sqrt(n_components) * components


class BaseRandomProjection(TransformerMixin, BaseEstimator, metaclass=ABCMeta):
    """Base class for random projections.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """

    @abstractmethod
    def __init__(self, n_components='auto', *, eps=0.1, dense_output=False,
                 random_state=None):
        self.n_components = n_components
        self.eps = eps
        self.dense_output = dense_output
        self.random_state = random_state

    @abstractmethod
    def _make_random_matrix(self, n_components, n_features):
        """ Generate the random projection matrix

        Parameters
        ----------
        n_components : int,
            Dimensionality of the target projection space.

        n_features : int,
            Dimensionality of the original source space.

        Returns
        -------
        components : numpy array or CSR matrix [n_components, n_features]
            The generated random matrix.

        """

    def fit(self, X, y=None):
        """Generate a sparse random projection matrix

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            Training set: only the shape is used to find optimal random
            matrix dimensions based on the theory referenced in the
            afore mentioned papers.

        y
            Ignored

        Returns
        -------
        self

        """
        X = self._validate_data(X, accept_sparse=['csr', 'csc'])

        n_samples, n_features = X.shape

        if self.n_components == 'auto':
            self.n_components_ = johnson_lindenstrauss_min_dim(
                n_samples=n_samples, eps=self.eps)

            if self.n_components_ <= 0:
                raise ValueError(
                    'eps=%f and n_samples=%d lead to a target dimension of '
                    '%d which is invalid' % (
                        self.eps, n_samples, self.n_components_))

            elif self.n_components_ > n_features:
                raise ValueError(
                    'eps=%f and n_samples=%d lead to a target dimension of '
                    '%d which is larger than the original space with '
                    'n_features=%d' % (self.eps, n_samples, self.n_components_,
                                       n_features))
        else:
            if self.n_components <= 0:
                raise ValueError("n_components must be greater than 0, got %s"
                                 % self.n_components)

            elif self.n_components > n_features:
                warnings.warn(
                    "The number of components is higher than the number of"
                    " features: n_features < n_components (%s < %s)."
                    "The dimensionality of the problem will not be reduced."
                    % (n_features, self.n_components),
                    DataDimensionalityWarning)

            self.n_components_ = self.n_components

        # Generate a projection matrix of size [n_components, n_features]
        self.components_ = self._make_random_matrix(self.n_components_,
                                                    n_features)

        # Check contract
        assert self.components_.shape == (self.n_components_, n_features), (
                'An error has occurred the self.components_ matrix has '
                ' not the proper shape.')

        return self

    def transform(self, X):
        """Project the data by using matrix product with the random matrix

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            The input data to project into a smaller dimensional space.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.
        """
        X = check_array(X, accept_sparse=['csr', 'csc'])

        check_is_fitted(self)

        if X.shape[1] != self.components_.shape[1]:
            raise ValueError(
                'Impossible to perform projection:'
                'X at fit stage had a different number of features. '
                '(%s != %s)' % (X.shape[1], self.components_.shape[1]))

        X_new = safe_sparse_dot(X, self.components_.T,
                                dense_output=self.dense_output)
        return X_new


class GaussianRandomProjection(BaseRandomProjection):
    """Reduce dimensionality through Gaussian random projection

    The components of the random matrix are drawn from N(0, 1 / n_components).

    Read more in the :ref:`User Guide <gaussian_random_matrix>`.

    .. versionadded:: 0.13

    Parameters
    ----------
    n_components : int or 'auto', optional (default = 'auto')
        Dimensionality of the target projection space.

        n_components can be automatically adjusted according to the
        number of samples in the dataset and the bound given by the
        Johnson-Lindenstrauss lemma. In that case the quality of the
        embedding is controlled by the ``eps`` parameter.

        It should be noted that Johnson-Lindenstrauss lemma can yield
        very conservative estimated of the required number of components
        as it makes no assumption on the structure of the dataset.

    eps : strictly positive float, optional (default=0.1)
        Parameter to control the quality of the embedding according to
        the Johnson-Lindenstrauss lemma when n_components is set to
        'auto'.

        Smaller values lead to better embedding and higher number of
        dimensions (n_components) in the target projection space.

    random_state : int, RandomState instance or None, optional (default=None)
        Controls the pseudo random number generator used to generate the
        projection matrix at fit time.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    n_components_ : int
        Concrete number of components computed when n_components="auto".

    components_ : numpy array of shape [n_components, n_features]
        Random matrix used for the projection.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.random_projection import GaussianRandomProjection
    >>> rng = np.random.RandomState(42)
    >>> X = rng.rand(100, 10000)
    >>> transformer = GaussianRandomProjection(random_state=rng)
    >>> X_new = transformer.fit_transform(X)
    >>> X_new.shape
    (100, 3947)

    See Also
    --------
    SparseRandomProjection

    """
    @_deprecate_positional_args
    def __init__(self, n_components='auto', *, eps=0.1, random_state=None):
        super().__init__(
            n_components=n_components,
            eps=eps,
            dense_output=True,
            random_state=random_state)

    def _make_random_matrix(self, n_components, n_features):
        """ Generate the random projection matrix

        Parameters
        ----------
        n_components : int,
            Dimensionality of the target projection space.

        n_features : int,
            Dimensionality of the original source space.

        Returns
        -------
        components : numpy array or CSR matrix [n_components, n_features]
            The generated random matrix.

        """
        random_state = check_random_state(self.random_state)
        return _gaussian_random_matrix(n_components,
                                       n_features,
                                       random_state=random_state)


class SparseRandomProjection(BaseRandomProjection):
    """Reduce dimensionality through sparse random projection

    Sparse random matrix is an alternative to dense random
    projection matrix that guarantees similar embedding quality while being
    much more memory efficient and allowing faster computation of the
    projected data.

    If we note `s = 1 / density` the components of the random matrix are
    drawn from:

      - -sqrt(s) / sqrt(n_components)   with probability 1 / 2s
      -  0                              with probability 1 - 1 / s
      - +sqrt(s) / sqrt(n_components)   with probability 1 / 2s

    Read more in the :ref:`User Guide <sparse_random_matrix>`.

    .. versionadded:: 0.13

    Parameters
    ----------
    n_components : int or 'auto', optional (default = 'auto')
        Dimensionality of the target projection space.

        n_components can be automatically adjusted according to the
        number of samples in the dataset and the bound given by the
        Johnson-Lindenstrauss lemma. In that case the quality of the
        embedding is controlled by the ``eps`` parameter.

        It should be noted that Johnson-Lindenstrauss lemma can yield
        very conservative estimated of the required number of components
        as it makes no assumption on the structure of the dataset.

    density : float in range ]0, 1], optional (default='auto')
        Ratio of non-zero component in the random projection matrix.

        If density = 'auto', the value is set to the minimum density
        as recommended by Ping Li et al.: 1 / sqrt(n_features).

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

    eps : strictly positive float, optional, (default=0.1)
        Parameter to control the quality of the embedding according to
        the Johnson-Lindenstrauss lemma when n_components is set to
        'auto'.

        Smaller values lead to better embedding and higher number of
        dimensions (n_components) in the target projection space.

    dense_output : boolean, optional (default=False)
        If True, ensure that the output of the random projection is a
        dense numpy array even if the input and random projection matrix
        are both sparse. In practice, if the number of components is
        small the number of zero components in the projected data will
        be very small and it will be more CPU and memory efficient to
        use a dense representation.

        If False, the projected data uses a sparse representation if
        the input is sparse.

    random_state : int, RandomState instance or None, optional (default=None)
        Controls the pseudo random number generator used to generate the
        projection matrix at fit time.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    n_components_ : int
        Concrete number of components computed when n_components="auto".

    components_ : CSR matrix with shape [n_components, n_features]
        Random matrix used for the projection.

    density_ : float in range 0.0 - 1.0
        Concrete density computed from when density = "auto".

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.random_projection import SparseRandomProjection
    >>> rng = np.random.RandomState(42)
    >>> X = rng.rand(100, 10000)
    >>> transformer = SparseRandomProjection(random_state=rng)
    >>> X_new = transformer.fit_transform(X)
    >>> X_new.shape
    (100, 3947)
    >>> # very few components are non-zero
    >>> np.mean(transformer.components_ != 0)
    0.0100...

    See Also
    --------
    GaussianRandomProjection

    References
    ----------

    .. [1] Ping Li, T. Hastie and K. W. Church, 2006,
           "Very Sparse Random Projections".
           https://web.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf

    .. [2] D. Achlioptas, 2001, "Database-friendly random projections",
           https://users.soe.ucsc.edu/~optas/papers/jl.pdf

    """
    @_deprecate_positional_args
    def __init__(self, n_components='auto', *, density='auto', eps=0.1,
                 dense_output=False, random_state=None):
        super().__init__(
            n_components=n_components,
            eps=eps,
            dense_output=dense_output,
            random_state=random_state)

        self.density = density

    def _make_random_matrix(self, n_components, n_features):
        """ Generate the random projection matrix

        Parameters
        ----------
        n_components : int,
            Dimensionality of the target projection space.

        n_features : int,
            Dimensionality of the original source space.

        Returns
        -------
        components : numpy array or CSR matrix [n_components, n_features]
            The generated random matrix.

        """
        random_state = check_random_state(self.random_state)
        self.density_ = _check_density(self.density, n_features)
        return _sparse_random_matrix(n_components,
                                     n_features,
                                     density=self.density_,
                                     random_state=random_state)
