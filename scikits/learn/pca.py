""" Principal Component Analysis
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD Style.

import numpy as np
from scipy import linalg

from .base import BaseEstimator
from .utils.extmath import fast_logdet
from .utils.extmath import fast_svd
from .utils.extmath import safe_sparse_dot


def _assess_dimension_(spectrum, rank, n_samples, dim):
    """Compute the likelihood of a rank rank dataset

    The dataset is assumed to be embedded in gaussian noise of shape(n,
    dimf) having spectrum spectrum.

    Parameters
    ----------
    spectrum: array of shape (n)
        data spectrum
    rank: int,
        tested rank value
    n_samples: int,
        number of samples
    dim: int,
        embedding/empirical dimension

    Returns
    -------
    ll: float,
        The log-likelihood

    Notes
    -----
    This implements the method of Thomas P. Minka:
    Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
    """
    if rank > dim:
        raise ValueError("the dimension cannot exceed dim")
    from scipy.special import gammaln

    pu = -rank * np.log(2)
    for i in range(rank):
        pu += gammaln((dim - i) / 2) - np.log(np.pi) * (dim - i) / 2

    pl = np.sum(np.log(spectrum[:rank]))
    pl = -pl * n_samples / 2

    if rank == dim:
        pv = 0
        v = 1
    else:
        v = np.sum(spectrum[rank:dim]) / (dim - rank)
        pv = -np.log(v) * n_samples * (dim - rank) / 2

    m = dim * rank - rank * (rank + 1) / 2
    pp = np.log(2 * np.pi) * (m + rank + 1) / 2

    pa = 0
    spectrum_ = spectrum.copy()
    spectrum_[rank:dim] = v
    for i in range(rank):
        for j in range(i + 1, dim):
            pa += (np.log((spectrum[i] - spectrum[j])
                          * (1. / spectrum_[j] - 1. / spectrum_[i]))
                   + np.log(n_samples))

    ll = pu + pl + pv + pp - pa / 2 - rank * np.log(n_samples) / 2

    return ll


def _infer_dimension_(spectrum, n, p):
    """This method infers the dimension of a dataset of shape (n, p)

    The dataset is described by its spectrum `spectrum`.
    """
    ll = []
    for rank in range(min(n, p, len(spectrum))):
        ll.append(_assess_dimension_(spectrum, rank, n, p))
    ll = np.array(ll)
    return ll.argmax()


class PCA(BaseEstimator):
    """Principal component analysis (PCA)

    Linear dimensionality reduction using Singular Value Decomposition of the
    data and keeping only the most significant singular vectors to project the
    data to a lower dimensional space.

    This implementation uses the scipy.linalg implementation of the singular
    value decomposition. It only works for dense arrays and is not scalable to
    large dimensional data.

    The time complexity of this implementation is O(n ** 3) assuming
    n ~ n_samples ~ n_features.

    Parameters
    ----------
    n_components: int, none or string
        Number of components to keep.
        if n_components is not set all components are kept:
            n_components == min(n_samples, n_features)

        if n_components == 'mle', Minka's MLE is used to guess the dimension

        if 0 < n_components < 1, select the number of components such that
                                 the explained variance ratio is greater
                                 than n_components

    copy: bool
        If False, data passed to fit are overwritten

    whiten: bool, optional
        When True (False by default) the components_ vectors are divided
        by n_samples times singular values to ensure uncorrelated outputs
        with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making there data respect some hard-wired assumptions.

    Attributes
    ----------
    components_: array, [n_features, n_components]
        Components with maximum variance.

    explained_variance_ratio_: array, [n_components]
        Percentage of variance explained by each of the selected components.
        k is not set then all components are stored and the sum of
        explained variances is equal to 1.0

    Notes
    -----
    For n_components='mle', this class uses the method of Thomas P. Minka:
    Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.pca import PCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = PCA(n_components=2)
    >>> pca.fit(X)
    PCA(copy=True, n_components=2, whiten=False)
    >>> print pca.explained_variance_ratio_
    [ 0.99244289  0.00755711]

    See also
    --------
    ProbabilisticPCA
    RandomizedPCA

    """
    def __init__(self, n_components=None, copy=True, whiten=False):
        self.n_components = n_components
        self.copy = copy
        self.whiten = whiten

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
        X = np.atleast_2d(X)
        n_samples, n_features = X.shape
        if self.copy:
            X = X.copy()
        # Center data
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_
        U, S, V = linalg.svd(X, full_matrices=False)
        self.explained_variance_ = (S ** 2) / n_samples
        self.explained_variance_ratio_ = self.explained_variance_ / \
                                        self.explained_variance_.sum()

        if self.whiten:
            n = X.shape[0]
            self.components_ = np.dot(V.T, np.diag(1.0 / S)) * np.sqrt(n)
        else:
            self.components_ = V.T

        if self.n_components == 'mle':
            self.n_components = _infer_dimension_(self.explained_variance_,
                                            n_samples, X.shape[1])

        elif 0 < self.n_components and self.n_components < 1.0:
            # number of components for which the cumulated explained variance
            # percentage is superior to the desired threshold
            n_remove = np.sum(self.explained_variance_ratio_.cumsum() >=
                              self.n_components) - 1
            self.n_components = n_features - n_remove

        if self.n_components is not None:
            self.components_ = self.components_[:, :self.n_components]
            self.explained_variance_ = \
                    self.explained_variance_[:self.n_components]
            self.explained_variance_ratio_ = \
                    self.explained_variance_ratio_[:self.n_components]

        return self

    def transform(self, X):
        """Apply the dimension reduction learned on the train data."""
        X_transformed = X - self.mean_
        X_transformed = np.dot(X_transformed, self.components_)
        return X_transformed

    def inverse_transform(self, X):
        """Return an input X_original whose transform would be X

        Note: if whitening is enabled, inverse_transform does not compute the
        exact inverse operation as transform.
        """
        return np.dot(X, self.components_.T) + self.mean_


class ProbabilisticPCA(PCA):
    """Additional layer on top of PCA that add a probabilistic evaluation

    """ + PCA.__doc__

    def fit(self, X, y=None, homoscedastic=True):
        """Additionally to PCA.fit, learns a covariance model

        Parameters
        ----------
        X: array of shape(n_samples, n_dim)
            The data to fit
        homoscedastic: bool, optional,
            If True, average variance across remaining dimensions
        """
        PCA.fit(self, X)
        self.dim = X.shape[1]
        Xr = X - self.mean_
        Xr -= np.dot(np.dot(Xr, self.components_), self.components_.T)
        n_samples = X.shape[0]
        if self.dim <= self.n_components:
            delta = np.zeros(self.dim)
        elif homoscedastic:
            delta = (Xr ** 2).sum() * np.ones(self.dim) \
                    / (n_samples * self.dim)
        else:
            delta = (Xr ** 2).mean(0) / (self.dim - self.n_components)
        self.covariance_ = np.diag(delta)
        for k in range(self.n_components):
            add_cov = np.dot(
                self.components_[:, k:k + 1], self.components_[:, k:k + 1].T)
            self.covariance_ += self.explained_variance_[k] * add_cov
        return self

    def score(self, X):
        """Return a score associated to new data

        Parameters
        ----------
        X: array of shape(n_samples, n_dim)
            The data to test

        Returns
        -------
        ll: array of shape (n_samples),
            log-likelihood of each row of X under the current model
        """
        Xr = X - self.mean_
        log_like = np.zeros(X.shape[0])
        self.precision_ = np.linalg.inv(self.covariance_)
        for i in range(X.shape[0]):
            log_like[i] = -.5 * np.dot(np.dot(self.precision_, Xr[i]), Xr[i])
        log_like += fast_logdet(self.precision_) - \
                                    self.dim / 2 * np.log(2 * np.pi)
        return log_like


class RandomizedPCA(BaseEstimator):
    """Principal component analysis (PCA) using randomized SVD

    Linear dimensionality reduction using approximated Singular Value
    Decomposition of the data and keeping only the most significant
    singular vectors to project the data to a lower dimensional space.

    This implementation uses a randomized SVD implementation and can
    handle both scipy.sparse and numpy dense arrays as input.

    Parameters
    ----------
    n_components: int
        Maximum number of components to keep: default is 50.

    copy: bool
        If False, data passed to fit are overwritten

    iterated_power: int, optional
        Number of iteration for the power method. 3 by default.

    whiten: bool, optional
        When True (False by default) the components_ vectors are divided
        by the singular values to ensure uncorrelated outputs with unit
        component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making there data respect some hard-wired assumptions.

    Attributes
    ----------
    components_: array, [n_features, n_components]
        Components with maximum variance.

    explained_variance_ratio_: array, [n_components]
        Percentage of variance explained by each of the selected components.
        k is not set then all components are stored and the sum of
        explained variances is equal to 1.0

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.pca import RandomizedPCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = RandomizedPCA(n_components=2)
    >>> pca.fit(X)
    RandomizedPCA(copy=True, n_components=2, iterated_power=3, whiten=False)
    >>> print pca.explained_variance_ratio_
    [ 0.99244289  0.00755711]

    See also
    --------
    PCA
    ProbabilisticPCA

    Notes
    -------
    References:

    * Finding structure with randomness: Stochastic algorithms for
      constructing approximate matrix decompositions Halko, et al., 2009
      (arXiv:909)

    * A randomized algorithm for the decomposition of matrices
      Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert


    """

    def __init__(self, n_components, copy=True, iterated_power=3,
                 whiten=False):
        self.n_components = n_components
        self.copy = copy
        self.iterated_power = iterated_power
        self.whiten = whiten
        self.mean_ = None

    def fit(self, X, y=None, **params):
        """Fit the model to the data X.

        Parameters
        ----------
        X: array-like or scipy.sparse matrix, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._set_params(**params)
        if not hasattr(X, 'todense'):
            X = np.atleast_2d(X)

        n_samples = X.shape[0]

        if self.copy:
            X = X.copy()

        if not hasattr(X, 'todense'):
            # not a sparse matrix, ensure this is a 2D array
            X = np.atleast_2d(X)

            # Center data
            self.mean_ = np.mean(X, axis=0)
            X -= self.mean_

        U, S, V = fast_svd(X, self.n_components, q=self.iterated_power)

        self.explained_variance_ = (S ** 2) / n_samples
        self.explained_variance_ratio_ = self.explained_variance_ / \
                                        self.explained_variance_.sum()

        if self.whiten:
            n = X.shape[0]
            self.components_ = np.dot(V.T, np.diag(1.0 / S)) * np.sqrt(n)
        else:
            self.components_ = V.T

        return self

    def transform(self, X):
        """Apply the dimension reduction learned on the training data."""
        if self.mean_ is not None:
            X = X - self.mean_

        X = safe_sparse_dot(X, self.components_)
        return X

    def inverse_transform(self, X):
        """Return an reconstructed input whose transform would be X"""
        X_original = safe_sparse_dot(X, self.components_.T)
        if self.mean_ is not None:
            X_original = X_original + self.mean_
        return X_original
