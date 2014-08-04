""" Principal Component Analysis
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Mathieu Blondel <mathieu@mblondel.org>
#         Denis A. Engemann <d.engemann@fz-juelich.de>
#
# License: BSD 3 clause

from math import log, sqrt
import warnings

import numpy as np
from scipy import linalg
from scipy import sparse
from scipy.special import gammaln

from ..base import BaseEstimator, TransformerMixin
from ..utils import array2d, check_random_state, as_float_array
from ..utils import atleast2d_or_csr
from ..utils import deprecated
from ..utils.sparsefuncs import mean_variance_axis0
from ..utils.extmath import (fast_logdet, safe_sparse_dot, randomized_svd,
                             fast_dot)


def _assess_dimension_(spectrum, rank, n_samples, n_features):
    """Compute the likelihood of a rank ``rank`` dataset

    The dataset is assumed to be embedded in gaussian noise of shape(n,
    dimf) having spectrum ``spectrum``.

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
    This implements the method of `Thomas P. Minka:
    Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604`
    """
    if rank > len(spectrum):
        raise ValueError("The tested rank cannot exceed the rank of the"
                         " dataset")

    pu = -rank * log(2.)
    for i in range(rank):
        pu += (gammaln((n_features - i) / 2.)
               - log(np.pi) * (n_features - i) / 2.)

    pl = np.sum(np.log(spectrum[:rank]))
    pl = -pl * n_samples / 2.

    if rank == n_features:
        pv = 0
        v = 1
    else:
        v = np.sum(spectrum[rank:]) / (n_features - rank)
        pv = -np.log(v) * n_samples * (n_features - rank) / 2.

    m = n_features * rank - rank * (rank + 1.) / 2.
    pp = log(2. * np.pi) * (m + rank + 1.) / 2.

    pa = 0.
    spectrum_ = spectrum.copy()
    spectrum_[rank:n_features] = v
    for i in range(rank):
        for j in range(i + 1, len(spectrum)):
            pa += log((spectrum[i] - spectrum[j]) *
                      (1. / spectrum_[j] - 1. / spectrum_[i])) + log(n_samples)

    ll = pu + pl + pv + pp - pa / 2. - rank * log(n_samples) / 2.

    return ll


def _infer_dimension_(spectrum, n_samples, n_features):
    """Infers the dimension of a dataset of shape (n_samples, n_features)

    The dataset is described by its spectrum `spectrum`.
    """
    n_spectrum = len(spectrum)
    ll = np.empty(n_spectrum)
    for rank in range(n_spectrum):
        ll[rank] = _assess_dimension_(spectrum, rank, n_samples, n_features)
    return ll.argmax()


class PCA(BaseEstimator, TransformerMixin):
    """Principal component analysis (PCA)

    Linear dimensionality reduction using Singular Value Decomposition of the
    data and keeping only the most significant singular vectors to project the
    data to a lower dimensional space.

    This implementation uses the scipy.linalg implementation of the singular
    value decomposition. It only works for dense arrays and is not scalable to
    large dimensional data.

    The time complexity of this implementation is ``O(n ** 3)`` assuming
    n ~ n_samples ~ n_features.

    Parameters
    ----------
    n_components : int, None or string
        Number of components to keep.
        if n_components is not set all components are kept::

            n_components == min(n_samples, n_features)

        if n_components == 'mle', Minka\'s MLE is used to guess the dimension
        if ``0 < n_components < 1``, select the number of components such that
        the amount of variance that needs to be explained is greater than the
        percentage specified by n_components

    copy : bool
        If False, data passed to fit are overwritten and running
        fit(X).transform(X) will not yield the expected results,
        use fit_transform(X) instead.

    whiten : bool, optional
        When True (False by default) the `components_` vectors are divided
        by n_samples times singular values to ensure uncorrelated outputs
        with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making there data respect some hard-wired assumptions.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Components with maximum variance.

    `explained_variance_ratio_` : array, [n_components]
        Percentage of variance explained by each of the selected components. \
        k is not set then all components are stored and the sum of explained \
        variances is equal to 1.0

    `mean_` : array, [n_features]
        Per-feature empirical mean, estimated from the training set.

    `n_components_` : int
        The estimated number of components. Relevant when n_components is set
        to 'mle' or a number between 0 and 1 to select using explained
        variance.

    `noise_variance_` : float
        The estimated noise covariance following the Probabilistic PCA model
        from Tipping and Bishop 1999. See "Pattern Recognition and
        Machine Learning" by C. Bishop, 12.2.1 p. 574 or
        http://www.miketipping.com/papers/met-mppca.pdf. It is required to
        computed the estimated data covariance and score samples.

    Notes
    -----
    For n_components='mle', this class uses the method of `Thomas P. Minka:
    Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604`

    Implements the probabilistic PCA model from:
    M. Tipping and C. Bishop, Probabilistic Principal Component Analysis,
    Journal of the Royal Statistical Society, Series B, 61, Part 3, pp. 611-622
    via the score and score_samples methods.
    See http://www.miketipping.com/papers/met-mppca.pdf

    Due to implementation subtleties of the Singular Value Decomposition (SVD),
    which is used in this implementation, running fit twice on the same matrix
    can lead to principal components with signs flipped (change in direction).
    For this reason, it is important to always use the same estimator object to
    transform data in a consistent fashion.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.decomposition import PCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = PCA(n_components=2)
    >>> pca.fit(X)
    PCA(copy=True, n_components=2, whiten=False)
    >>> print(pca.explained_variance_ratio_) # doctest: +ELLIPSIS
    [ 0.99244...  0.00755...]

    See also
    --------
    ProbabilisticPCA
    RandomizedPCA
    KernelPCA
    SparsePCA
    TruncatedSVD
    """
    def __init__(self, n_components=None, copy=True, whiten=False):
        self.n_components = n_components
        self.copy = copy
        self.whiten = whiten

    def fit(self, X, y=None):
        """Fit the model with X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._fit(X)
        return self

    def fit_transform(self, X, y=None):
        """Fit the model with X and apply the dimensionality reduction on X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        U, S, V = self._fit(X)
        U = U[:, :self.n_components_]

        if self.whiten:
            # X_new = X * V / S * sqrt(n_samples) = U * sqrt(n_samples)
            U *= sqrt(X.shape[0])
        else:
            # X_new = X * V = U * S * V^T * V = U * S
            U *= S[:self.n_components_]

        return U

    def _fit(self, X):
        """Fit the model on X

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        Returns
        -------
        U, s, V : ndarrays
            The SVD of the input data, copied and centered when
            requested.
        """
        X = array2d(X)
        n_samples, n_features = X.shape
        X = as_float_array(X, copy=self.copy)
        # Center data
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_
        U, S, V = linalg.svd(X, full_matrices=False)
        explained_variance_ = (S ** 2) / n_samples
        explained_variance_ratio_ = (explained_variance_ /
                                     explained_variance_.sum())

        if self.whiten:
            components_ = V / (S[:, np.newaxis] / sqrt(n_samples))
        else:
            components_ = V

        n_components = self.n_components
        if n_components is None:
            n_components = n_features
        elif n_components == 'mle':
            if n_samples < n_features:
                raise ValueError("n_components='mle' is only supported "
                                 "if n_samples >= n_features")

            n_components = _infer_dimension_(explained_variance_,
                                             n_samples, n_features)
        elif not 0 <= n_components <= n_features:
            raise ValueError("n_components=%r invalid for n_features=%d"
                             % (n_components, n_features))

        if 0 < n_components < 1.0:
            # number of components for which the cumulated explained variance
            # percentage is superior to the desired threshold
            ratio_cumsum = explained_variance_ratio_.cumsum()
            n_components = np.sum(ratio_cumsum < n_components) + 1

        # Compute noise covariance using Probabilistic PCA model
        # The sigma2 maximum likelihood (cf. eq. 12.46)
        if n_components < n_features:
            self.noise_variance_ = explained_variance_[n_components:].mean()
        else:
            self.noise_variance_ = 0.

        # store n_samples to revert whitening when getting covariance
        self.n_samples_ = n_samples

        self.components_ = components_[:n_components]
        self.explained_variance_ = explained_variance_[:n_components]
        explained_variance_ratio_ = explained_variance_ratio_[:n_components]
        self.explained_variance_ratio_ = explained_variance_ratio_
        self.n_components_ = n_components

        return (U, S, V)

    def get_covariance(self):
        """Compute data covariance with the generative model.

        ``cov = components_.T * S**2 * components_ + sigma2 * eye(n_features)``
        where  S**2 contains the explained variances.

        Returns
        -------
        cov : array, shape=(n_features, n_features)
            Estimated covariance of data.
        """
        components_ = self.components_
        exp_var = self.explained_variance_
        if self.whiten:
            components_ = components_ * np.sqrt(exp_var[:, np.newaxis])
        exp_var_diff = np.maximum(exp_var - self.noise_variance_, 0.)
        cov = np.dot(components_.T * exp_var_diff, components_)
        cov.flat[::len(cov) + 1] += self.noise_variance_  # modify diag inplace
        return cov

    def get_precision(self):
        """Compute data precision matrix with the generative model.

        Equals the inverse of the covariance but computed with
        the matrix inversion lemma for efficiency.

        Returns
        -------
        precision : array, shape=(n_features, n_features)
            Estimated precision of data.
        """
        n_features = self.components_.shape[1]

        # handle corner cases first
        if self.n_components_ == 0:
            return np.eye(n_features) / self.noise_variance_
        if self.n_components_ == n_features:
            return linalg.inv(self.get_covariance())

        # Get precision using matrix inversion lemma
        components_ = self.components_
        exp_var = self.explained_variance_
        if self.whiten:
            components_ = components_ * np.sqrt(exp_var[:, np.newaxis])
        exp_var_diff = np.maximum(exp_var - self.noise_variance_, 0.)
        precision = np.dot(components_, components_.T) / self.noise_variance_
        precision.flat[::len(precision) + 1] += 1. / exp_var_diff
        precision = np.dot(components_.T,
                           np.dot(linalg.inv(precision), components_))
        precision /= -(self.noise_variance_ ** 2)
        precision.flat[::len(precision) + 1] += 1. / self.noise_variance_
        return precision

    def transform(self, X):
        """Apply the dimensionality reduction on X.

        X is projected on the first principal components previous extracted
        from a training set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            New data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        X = array2d(X)
        if self.mean_ is not None:
            X = X - self.mean_
        X_transformed = fast_dot(X, self.components_.T)
        return X_transformed

    def inverse_transform(self, X):
        """Transform data back to its original space, i.e.,
        return an input X_original whose transform would be X

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            New data, where n_samples is the number of samples
            and n_components is the number of components.

        Returns
        -------
        X_original array-like, shape (n_samples, n_features)

        Notes
        -----
        If whitening is enabled, inverse_transform does not compute the
        exact inverse operation as transform.
        """
        return fast_dot(X, self.components_) + self.mean_

    def score_samples(self, X):
        """Return the log-likelihood of each sample

        See. "Pattern Recognition and Machine Learning"
        by C. Bishop, 12.2.1 p. 574
        or http://www.miketipping.com/papers/met-mppca.pdf

        Parameters
        ----------
        X: array, shape(n_samples, n_features)
            The data.

        Returns
        -------
        ll: array, shape (n_samples,)
            Log-likelihood of each sample under the current model
        """
        X = array2d(X)
        Xr = X - self.mean_
        n_features = X.shape[1]
        log_like = np.zeros(X.shape[0])
        precision = self.get_precision()
        log_like = -.5 * (Xr * (np.dot(Xr, precision))).sum(axis=1)
        log_like -= .5 * (n_features * log(2. * np.pi)
                          - fast_logdet(precision))
        return log_like

    def score(self, X, y=None):
        """Return the average log-likelihood of all samples

        See. "Pattern Recognition and Machine Learning"
        by C. Bishop, 12.2.1 p. 574
        or http://www.miketipping.com/papers/met-mppca.pdf

        Parameters
        ----------
        X: array, shape(n_samples, n_features)
            The data.

        Returns
        -------
        ll: float
            Average log-likelihood of the samples under the current model
        """
        return np.mean(self.score_samples(X))


@deprecated("ProbabilisticPCA will be removed in 0.16. WARNING: the "
            "covariance estimation was previously incorrect, your "
            "output might be different than under the previous versions. "
            "Use PCA that implements score and score_samples. To work with "
            "homoscedastic=False, you should use FactorAnalysis.")
class ProbabilisticPCA(PCA):
    """Additional layer on top of PCA that adds a probabilistic evaluation"""
    __doc__ += PCA.__doc__

    def fit(self, X, y=None, homoscedastic=True):
        """Additionally to PCA.fit, learns a covariance model

        Parameters
        ----------
        X : array of shape(n_samples, n_features)
            The data to fit

        homoscedastic : bool, optional,
            If True, average variance across remaining dimensions
        """
        PCA.fit(self, X)

        n_samples, n_features = X.shape
        n_components = self.n_components
        if n_components is None:
            n_components = n_features

        explained_variance = self.explained_variance_.copy()
        if homoscedastic:
            explained_variance -= self.noise_variance_

        # Make the low rank part of the estimated covariance
        self.covariance_ = np.dot(self.components_[:n_components].T *
                                  explained_variance,
                                  self.components_[:n_components])

        if n_features == n_components:
            delta = 0.
        elif homoscedastic:
            delta = self.noise_variance_
        else:
            Xr = X - self.mean_
            Xr -= np.dot(np.dot(Xr, self.components_.T), self.components_)
            delta = (Xr ** 2).mean(axis=0) / (n_features - n_components)

        # Add delta to the diagonal without extra allocation
        self.covariance_.flat[::n_features + 1] += delta

        return self

    def score(self, X, y=None):
        """Return a score associated to new data

        Parameters
        ----------
        X: array of shape(n_samples, n_features)
            The data to test

        Returns
        -------
        ll: array of shape (n_samples),
            log-likelihood of each row of X under the current model
        """
        Xr = X - self.mean_
        n_features = X.shape[1]
        log_like = np.zeros(X.shape[0])
        self.precision_ = linalg.inv(self.covariance_)
        log_like = -.5 * (Xr * (np.dot(Xr, self.precision_))).sum(axis=1)
        log_like -= .5 * (fast_logdet(self.covariance_)
                          + n_features * log(2. * np.pi))
        return log_like


class RandomizedPCA(BaseEstimator, TransformerMixin):
    """Principal component analysis (PCA) using randomized SVD

    Linear dimensionality reduction using approximated Singular Value
    Decomposition of the data and keeping only the most significant
    singular vectors to project the data to a lower dimensional space.

    Parameters
    ----------
    n_components : int, optional
        Maximum number of components to keep. When not given or None, this
        is set to n_features (the second dimension of the training data).

    copy : bool
        If False, data passed to fit are overwritten and running
        fit(X).transform(X) will not yield the expected results,
        use fit_transform(X) instead.

    iterated_power : int, optional
        Number of iterations for the power method. 3 by default.

    whiten : bool, optional
        When True (False by default) the `components_` vectors are divided
        by the singular values to ensure uncorrelated outputs with unit
        component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making their data respect some hard-wired assumptions.

    random_state : int or RandomState instance or None (default)
        Pseudo Random Number generator seed control. If None, use the
        numpy.random singleton.

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Components with maximum variance.

    `explained_variance_ratio_` : array, [n_components]
        Percentage of variance explained by each of the selected components. \
        k is not set then all components are stored and the sum of explained \
        variances is equal to 1.0

    `mean_` : array, [n_features]
        Per-feature empirical mean, estimated from the training set.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.decomposition import RandomizedPCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = RandomizedPCA(n_components=2)
    >>> pca.fit(X)                 # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    RandomizedPCA(copy=True, iterated_power=3, n_components=2,
           random_state=None, whiten=False)
    >>> print(pca.explained_variance_ratio_) # doctest: +ELLIPSIS
    [ 0.99244...  0.00755...]

    See also
    --------
    PCA
    ProbabilisticPCA
    TruncatedSVD

    References
    ----------

    .. [Halko2009] `Finding structure with randomness: Stochastic algorithms
      for constructing approximate matrix decompositions Halko, et al., 2009
      (arXiv:909)`

    .. [MRT] `A randomized algorithm for the decomposition of matrices
      Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert`

    Notes
    -----
    This class supports sparse matrix input for backward compatibility, but
    actually computes a truncated SVD instead of a PCA in that case (i.e. no
    centering is performed). This support is deprecated; use the class
    TruncatedSVD for sparse matrix support.

    """

    def __init__(self, n_components=None, copy=True, iterated_power=3,
                 whiten=False, random_state=None):
        self.n_components = n_components
        self.copy = copy
        self.iterated_power = iterated_power
        self.whiten = whiten
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model with X by extracting the first principal components.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._fit(X)
        return self

    def _fit(self, X):
        """Fit the model to the data X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        Returns
        -------
        X : ndarray, shape (n_samples, n_features)
            The input data, copied, centered and whitened when requested.
        """
        random_state = check_random_state(self.random_state)
        if sparse.issparse(X):
            warnings.warn("Sparse matrix support is deprecated in 0.15"
                          " and will be dropped in 0.17. In particular"
                          " computed explained variance is incorrect on"
                          " sparse data. Use TruncatedSVD instead.",
                          DeprecationWarning)
        else:
            # not a sparse matrix, ensure this is a 2D array
            X = np.atleast_2d(as_float_array(X, copy=self.copy))

        n_samples = X.shape[0]

        if sparse.issparse(X):
            self.mean_ = None
        else:
            # Center data
            self.mean_ = np.mean(X, axis=0)
            X -= self.mean_
        if self.n_components is None:
            n_components = X.shape[1]
        else:
            n_components = self.n_components

        U, S, V = randomized_svd(X, n_components,
                                 n_iter=self.iterated_power,
                                 random_state=random_state)

        self.explained_variance_ = exp_var = (S ** 2) / n_samples
        if sparse.issparse(X):
            _, full_var = mean_variance_axis0(X)
            full_var = full_var.sum()
        else:
            full_var = np.var(X, axis=0).sum()
        self.explained_variance_ratio_ = exp_var / full_var

        if self.whiten:
            self.components_ = V / S[:, np.newaxis] * sqrt(n_samples)
        else:
            self.components_ = V

        return X

    def transform(self, X, y=None):
        """Apply dimensionality reduction on X.

        X is projected on the first principal components previous extracted
        from a training set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        # XXX remove scipy.sparse support here in 0.16
        X = atleast2d_or_csr(X)
        if self.mean_ is not None:
            X = X - self.mean_

        X = safe_sparse_dot(X, self.components_.T)
        return X

    def fit_transform(self, X, y=None):
        """Fit the model with X and apply the dimensionality reduction on X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        X = self._fit(atleast2d_or_csr(X))
        X = safe_sparse_dot(X, self.components_.T)
        return X

    def inverse_transform(self, X, y=None):
        """Transform data back to its original space.

        Returns an array X_original whose transform would be X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            New data, where n_samples in the number of samples
            and n_components is the number of components.

        Returns
        -------
        X_original array-like, shape (n_samples, n_features)

        Notes
        -----
        If whitening is enabled, inverse_transform does not compute the
        exact inverse operation of transform.
        """
        # XXX remove scipy.sparse support here in 0.16
        X_original = safe_sparse_dot(X, self.components_)
        if self.mean_ is not None:
            X_original = X_original + self.mean_
        return X_original
