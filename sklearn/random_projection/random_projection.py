"""Random Projection transformers

Random Projections are a simple and computationally efficient way to
reduce the dimensionality of the data by trading a controlled amount
of accuracy (as additional variance) for faster processing times and
smaller model sizes.

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
# Authors: Olivier Grisel <olivier.grisel@ensta.org>,
#          Arnaud Joly <a.joly@ulg.ac.be>
# License: Simple BSD

from __future__ import division
import warnings

import numpy as np
from numpy.testing import assert_equal
import scipy.sparse as sp

from sklearn.utils import check_random_state
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.random_projection._random_projection import sample_int


__all__ = ["BernouilliRandomProjection",
           "GaussianRandomProjection",
           "johnson_lindenstrauss_min_dim"]


def johnson_lindenstrauss_min_dim(n_samples, eps=0.1):
    """Find a 'safe' number of components to randomly project to

    The distortion introduced by a random projection `p` is asserted by
    the fact that `p` is defining an eps-embedding with good probability
    as defined by:

      (1 - eps) ||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps) ||u - v||^2

    Where u and v are any rows taken from a dataset of shape [n_samples,
    n_features] and p is a projection by a random gaussian N(0, 1) matrix
    with shape [n_components, n_features] (or a sparse Achlioptas matrix).

    The minimum number of components to guarantees the eps-embedding is
    given by:

      n_components >= 4 log(n_samples) / (eps^2 / 2 - eps^3 / 3)

    Note that the number of dimensions is independent of the original
    number of features but instead depends on the size of the dataset:
    the larger the dataset, the higher is the minimal dimensionality of
    random projection embedding.

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
    - http://en.wikipedia.org/wiki/Johnson%E2%80%93Lindenstrauss_lemma

    - An elementary proof of the Johnson-Lindenstrauss Lemma.
      Sanjoy Dasgupta and Anupam Gupta, 1999
      http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.45.3654

    """
    eps = np.asarray(eps)
    if np.any(eps <= 0.0) or np.any(eps > 1):
        raise ValueError(
            "The JL bound is defined for eps in (0, 1]: got %r" % eps)
    denominator = (eps ** 2 / 2) - (eps ** 3 / 3)
    return (4 * np.log(n_samples) / denominator).astype(np.int)


def _check_density(density, n_features):
    """Factorize density check according to Li et al."""
    if density is 'auto':
        density = min(1 / np.sqrt(n_features), 1 / 3.)

    elif density <= 0 or density > 1:
        raise ValueError("Expected density in range (0, 1], got: %r"
                         % density)
    return density


def _check_input_size(n_components, n_features):
    """Factorize argument checking for random matrix generation"""
    if n_components <= 0:
        raise ValueError("n_components must be strictly positive, got %d" %
                         n_components)
    if n_features <= 0:
        raise ValueError("n_features must be strictly positive, got %d" %
                         n_components)


def gaussian_random_matrix(n_components, n_features, random_state=None):
    """ Generate a dense gaussian random matrix.

    The components of the random matrix are drawn from
    N(0, 1.0 / n_components).

    Parameters
    ----------
    n_components: int,
        Dimensionality of the target projection space.

    n_features: int,
        Dimensionality of the original source space.

    random_state : int, RandomState instance or None (default)
        Control the pseudo random number generator used to generate the
        matrix at fit time.

    See Also
    --------
    bernouilli_random_matrix
    """
    _check_input_size(n_components, n_features)
    rng = check_random_state(random_state)
    rp_matrix = rng.normal(loc=0.0,
                           scale=1.0 / np.sqrt(n_components),
                           size=(n_components, n_features))
    return rp_matrix


def bernouilli_random_matrix(n_components, n_features, density='auto',
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
    n_components: int,
        Dimensionality of the target projection space.

    n_features: int,
        Dimensionality of the original source space.

    density: float in range (0, 1/3], optional
        Ratio of non-zero component in the random projection matrix.

        By default the value is set to the minimum density as recommended
        by Ping Li et al.: 1 / sqrt(n_features)

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

        Use density = 1 if you want a dense bernouilli matrix.

    random_state : integer, RandomState instance or None (default)
        Control the pseudo random number generator used to generate the
        matrix at fit time.

    See Also
    --------
    gaussian_random_matrix

    References
    ----------

    - Very Sparse Random Projections. Ping Li, Trevor Hastie
      and Kenneth Church, 2006
      http://www.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf

    - Database-friendly random projections, Dimitris Achlioptas, 2001
      http://www.cs.ucsc.edu/~optas/papers/jl.pdf

    """
    _check_input_size(n_components, n_features)
    density = _check_density(density, n_features)
    rng = check_random_state(random_state)

    if density == 1:
        # efficient implementation for dense bernouilli projection
        rp_matrix = rng.binomial(1, 0.5, (n_components, n_features)) * 2 - 1
        return 1 / np.sqrt(n_components) * rp_matrix

    else:
        # Generate location of non zero element
        indices = []
        offset = 0
        indptr = [offset]
        for i in xrange(n_components):
            # find the indices of the non-zero components for row i
            n_nonzero_i = rng.binomial(n_features, density)
            indices_i = sample_int(n_features, n_nonzero_i,
                                   random_state=rng)
            indices.append(indices_i)

            # among non zero components the probability of the sign is
            # 50%/50%
            offset += n_nonzero_i
            indptr.append(offset)

        indices = np.concatenate(indices)

        # Generate data
        data = rng.binomial(1, 0.5, size=len(indices)) * 2 - 1

        # build the CSR structure by concatenating the rows
        rp_matrix = sp.csr_matrix((data, indices, indptr),
                                  shape=(n_components, n_features))

        return np.sqrt(1 / density) / np.sqrt(n_components) * rp_matrix


class BaseRandomProjection(BaseEstimator, TransformerMixin):
    """Base class for random projections.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """
    def __init__(self, n_components='auto', density='auto', eps=0.5,
                 dense_output=False, random_state=None,
                 distribution="bernouilli"):
        self.n_components = n_components
        self.density = density
        self.eps = eps
        self.dense_output = dense_output
        self.random_state = random_state
        self.distribution = distribution

        self.components_ = None
        self.n_components_ = None
        self.density_ = None

    def fit(self, X, y=None):
        """Generate a sparse random projection matrix

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            Training set: only the shape is used to find optimal random
            matrix dimensions based on the theory referenced in the
            afore mentioned papers.

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        self

        """
        if not sp.issparse(X):
            X = np.atleast_2d(X)

        self.random_state = check_random_state(self.random_state)
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
                                 % self.n_components_)

            elif self.n_components > n_features:
                warnings.warn(
                    "The number of components is higher than the number of"
                    " features: n_features > n_components (%s > %s)."
                    "The dimensionality of the problem will not be reduced."
                    % (n_features, self.n_components))

            self.n_components_ = self.n_components

        self.density_ = _check_density(self.density, n_features)

        if self.distribution == "bernouilli":
            self.components_ = bernouilli_random_matrix(
                self.n_components_, n_features, density=self.density,
                random_state=self.random_state)

        elif self.distribution == "gaussian":
            if self.density != 1.0:
                raise NotImplementedError(
                    "Sparse gaussion projection not implemented."
                    "Density should be equal to 1.0 for Gaussian random "
                    "projection.")

            self.components_ = gaussian_random_matrix(
                self.n_components_, n_features, random_state=self.random_state)

        else:
            raise ValueError('Unknown specified distribution.')

        assert_equal(
            self.components_.shape,
            (self.n_components_, n_features),
            err_msg=('An error has occured the self.components_ matrix has not'
                     ' the proper shape.'))

        return self

    def transform(self, X, y=None):
        """Project the data by using matrix product with the random matrix

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape [n_samples, n_features]
            The input data to project into a smaller dimensional space.

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.

        """
        if self.components_ is None:
            raise ValueError('No random projection matrix had been fit.')

        if X.shape[1] != self.components_.shape[1]:
            raise ValueError(
                'Impossible to perform projection:'
                'X at fit stage had a different number of features.'
                '(%s != %s)' % (X.shape[1], self.components_.shape[1]))

        if not sp.issparse(X):
            X = np.atleast_2d(X)

        return safe_sparse_dot(X, self.components_.T,
                               dense_output=self.dense_output)


class GaussianRandomProjection(BaseRandomProjection):
    """Transformer to reduce the dimensionality with Gaussian random
    projection.

    The components of the random matrix are drawn from
        N(0, 1.0 / n_components).

    Parameters
    ----------
    n_components : int or 'auto', optional
        Dimensionality of the target projection space.

        n_components can be automatically adjusted according to the
        number of samples in the dataset and the bound given by the
        Johnson Lindenstrauss lemma. In that case the quality of the
        embedding is controlled by the ``eps`` parameter.

        It should be noted that Johnson Lindenstrauss lemma can yield
        very conservative estimated of the required number of components
        as it makes no assumption on the structure of the dataset.

    density : float in range (0, 1], optional
        Ratio of non-zero component in the random projection matrix.

        By default the value is set to the minimum density as recommended
        by Ping Li et al.: 1 / sqrt(n_features)

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

        Use density = 1 if you want to get a dense Bernouilli random matrix.

    eps : strictly positive float, optional, default 0.1
        Parameter to control the quality of the embedding according to
        the Johnson-Lindenstrauss lemma when n_components is set to
        'auto'.

        Smaller values lead to better embedding and higher number of
        dimensions (n_components) in the target projection space.

    random_state : integer, RandomState instance or None (default)
        Control the pseudo random number generator used to generate the
        matrix at fit time.

    Attributes
    ----------
    n_component_: int
        Concrete number of components computed when n_components="auto".

    components_: CSR matrix with shape [n_components, n_features]
        Random matrix used for the projection.

    """
    def __init__(self, n_components='auto', eps=0.5, random_state=None):
        super(GaussianRandomProjection, self).__init__(
            n_components=n_components,
            density=1.0,
            eps=eps,
            dense_output=True,
            random_state=random_state,
            distribution="gaussian")


class BernouilliRandomProjection(BaseRandomProjection):
    """Transformer to reduce the dimensionality with bernouilli random
    projection.

    Bernouilli random matrix is an alternative to Gaussian projection.

    Sparse bernouilli random matrix is an alternative to dense random
    projection matrix that guarantees similar embedding quality while being
    much more memory efficient and allowing faster computation of the
    projected data.

    If we note `s = 1 / density` the components of the random matrix are
    drawn from:

      - -sqrt(s) / sqrt(n_components)   with probability 1 / 2s
      -  0                              with probability 1 - 1 / s
      - +sqrt(s) / sqrt(n_components)   with probability 1 / 2s

    Parameters
    ----------
    n_components : int or 'auto', optional
        Dimensionality of the target projection space.

        n_components can be automatically adjusted according to the
        number of samples in the dataset and the bound given by the
        Johnson Lindenstrauss lemma. In that case the quality of the
        embedding is controlled by the ``eps`` parameter.

        It should be noted that Johnson Lindenstrauss lemma can yield
        very conservative estimated of the required number of components
        as it makes no assumption on the structure of the dataset.

    density : float in range (0, 1], optional
        Ratio of non-zero component in the random projection matrix.

        By default the value is set to the minimum density as recommended
        by Ping Li et al.: 1 / sqrt(n_features)

        Use density = 1 / 3.0 if you want to reproduce the results from
        Achlioptas, 2001.

        Use density = 1 if you want dense bernouilli random projection.

    eps : strictly positive float, optional, default 0.1
        Parameter to control the quality of the embedding according to
        the Johnson-Lindenstrauss lemma when n_components is set to
        'auto'.

        Smaller values lead to better embedding and higher number of
        dimensions (n_components) in the target projection space.

    dense_output : boolean, False by default
        If True, ensure that the output of the random projection is a
        dense numpy array even if the input and random projection matrix
        are both sparse. In practice, if the number of components is
        small the number of zero components in the projected data will
        be very small and it will be more CPU and memory efficient to
        use a dense representation.

        If False, the projected data uses a sparse representation if
        the input is sparse.

    random_state : integer, RandomState instance or None (default)
        Control the pseudo random number generator used to generate the
        matrix at fit time.

    Attributes
    ----------
    n_component_: int
        Concrete number of components computed when n_components="auto".

    components_: CSR matrix with shape [n_components, n_features]
        Random matrix used for the projection.

    density_: float in range 0.0 - 1.0
        Concrete density computed from when density="auto".

    References
    ----------
    - Very Sparse Random Projections. Ping Li, Trevor Hastie
      and Kenneth Church, 2006
      http://www.stanford.edu/~hastie/Papers/Ping/KDD06_rp.pdf

    - Database-friendly random projections, Dimitris Achlioptas, 2001
      http://www.cs.ucsc.edu/~optas/papers/jl.pdf

    """
    def __init__(self, n_components='auto', density='auto', eps=0.5,
                 dense_output=False, random_state=None):
        super(BernouilliRandomProjection, self).__init__(
            n_components=n_components,
            density=density,
            eps=eps,
            dense_output=dense_output,
            random_state=random_state,
            distribution="bernouilli")
