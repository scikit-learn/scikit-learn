"""Incremental Principal Components Analysis."""

# Author: Kyle Kastner <kastnerkyle@gmail.com>
# License: BSD 3 clause

import numpy as np
from scipy import linalg

from .base import _BasePCA
from ..utils import check_array, gen_batches
from ..utils.extmath import svd_flip, _batch_mean_variance_update


class IncrementalPCA(_BasePCA):
    """Incremental principal components analysis (IPCA).

    Linear dimensionality reduction using Singular Value Decomposition of
    centered data, keeping only the most significant singular vectors to
    project the data to a lower dimensional space.

    Depending on the size of the input data, this algorithm can be much more
    memory efficient than a PCA.

    This algorithm has constant memory complexity, on the order
    of ``batch_size``, enabling use of np.memmap files without loading the
    entire file into memory.

    The computational overhead of each SVD is
    ``O(batch_size * n_features ** 2)``, but only 2 * batch_size samples
    remain in memory at a time. There will be ``n_samples / batch_size`` SVD
    computations to get the principal components, versus 1 large SVD of
    complexity ``O(n_samples * n_features ** 2)`` for PCA.

    Parameters
    ----------
    n_components : int or None, (default=None)
        Number of components to keep. If ``n_components `` is ``None``,
        then ``n_components`` is set to ``min(n_samples, n_features)``.

    batch_size : int or None, (default=None)
        The number of samples to use for each batch. Only used when calling
        ``fit``. If ``batch_size`` is ``None``, then ``batch_size``
        is inferred from the data and set to ``5 * n_features``, to provide a
        balance between approximation accuracy and memory consumption.

    copy : bool, (default=True)
        If False, X will be overwritten. ``copy=False`` can be used to
        save memory but is unsafe for general use.

    whiten : bool, optional
        When True (False by default) the ``components_`` vectors are divided
        by ``n_samples`` times ``components_`` to ensure uncorrelated outputs
        with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometimes
        improve the predictive accuracy of the downstream estimators by
        making data respect some hard-wired assumptions.

    Attributes
    ----------
    components_ : array, shape (n_components, n_features)
        Components with maximum variance.

    explained_variance_ : array, shape (n_components,)
        Variance explained by each of the selected components.

    explained_variance_ratio_ : array, shape (n_components,)
        Percentage of variance explained by each of the selected components.
        If all components are stored, the sum of explained variances is equal
        to 1.0

    mean_ : array, shape (n_features,)
        Per-feature empirical mean, aggregate over calls to ``partial_fit``.

    var_ : array, shape (n_features,)
        Per-feature empirical variance, aggregate over calls to ``partial_fit``.

    noise_variance_ : float
        The estimated noise covariance following the Probabilistic PCA model
        from Tipping and Bishop 1999. See "Pattern Recognition and
        Machine Learning" by C. Bishop, 12.2.1 p. 574 or
        http://www.miketipping.com/papers/met-mppca.pdf.

    n_components_ : int
        The estimated number of components. Relevant when ``n_components=None``.

    n_samples_seen_ : int
        The number of samples processed by the estimator. Will be reset on
        new calls to fit, but increments across ``partial_fit`` calls.

    Notes
    -----
    Implements the incremental PCA model from:
    `D. Ross, J. Lim, R. Lin, M. Yang, Incremental Learning for Robust Visual
    Tracking, International Journal of Computer Vision, Volume 77, Issue 1-3,
    pp. 125-141, May 2008.`
    See http://www.cs.toronto.edu/~dross/ivt/RossLimLinYang_ijcv.pdf

    This model is an extension of the Sequential Karhunen-Loeve Transform from:
    `A. Levy and M. Lindenbaum, Sequential Karhunen-Loeve Basis Extraction and
    its Application to Images, IEEE Transactions on Image Processing, Volume 9,
    Number 8, pp. 1371-1374, August 2000.`
    See http://www.cs.technion.ac.il/~mic/doc/skl-ip.pdf

    We have specifically abstained from an optimization used by authors of both
    papers, a QR decomposition used in specific situations to reduce the
    algorithmic complexity of the SVD. The source for this technique is
    `Matrix Computations, Third Edition, G. Holub and C. Van Loan, Chapter 5,
    section 5.4.4, pp 252-253.`. This technique has been omitted because it is
    advantageous only when decomposing a matrix with ``n_samples`` (rows)
    >= 5/3 * ``n_features`` (columns), and hurts the readability of the
    implemented algorithm. This would be a good opportunity for future
    optimization, if it is deemed necessary.

    References
    ----------
    D. Ross, J. Lim, R. Lin, M. Yang. Incremental Learning for Robust Visual
        Tracking, International Journal of Computer Vision, Volume 77,
        Issue 1-3, pp. 125-141, May 2008.

    G. Golub and C. Van Loan. Matrix Computations, Third Edition, Chapter 5,
        Section 5.4.4, pp. 252-253.

    See also
    --------
    PCA
    RandomizedPCA
    KernelPCA
    SparsePCA
    TruncatedSVD
    """

    def __init__(self, n_components=None, whiten=False, copy=True,
                 batch_size=None):
        self.n_components = n_components
        self.whiten = whiten
        self.copy = copy
        self.batch_size = batch_size

    def fit(self, X, y=None):
        """Fit the model with X, using minibatches of size batch_size.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples and
            n_features is the number of features.

        y: Passthrough for ``Pipeline`` compatibility.

        Returns
        -------
        self: object
            Returns the instance itself.
        """
        self.components_ = None
        self.mean_ = None
        self.singular_values_ = None
        self.explained_variance_ = None
        self.explained_variance_ratio_ = None
        self.noise_variance_ = None
        self.var_ = None
        self.n_samples_seen_ = 0
        X = check_array(X, dtype=np.float)
        n_samples, n_features = X.shape

        if self.batch_size is None:
            self.batch_size_ = 5 * n_features
        else:
            self.batch_size_ = self.batch_size

        for batch in gen_batches(n_samples, self.batch_size_):
            self.partial_fit(X[batch])
        return self

    def partial_fit(self, X):
        """Incremental fit with X. All of X is processed as a single batch.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples and
            n_features is the number of features.

        Returns
        -------
        self: object
            Returns the instance itself.
        """
        X = check_array(X, copy=self.copy, dtype=np.float)
        n_samples, n_features = X.shape
        if not hasattr(self, 'components_'):
            self.components_ = None

        if self.n_components is None:
            self.n_components_ = n_features
        elif not 1 <= self.n_components <= n_features:
            raise ValueError("n_components=%r invalid for n_features=%d, need "
                             "more rows than columns for IncrementalPCA "
                             "processing" % (self.n_components, n_features))
        else:
            self.n_components_ = self.n_components

        if (self.components_ is not None) and (self.components_.shape[0]
                                               != self.n_components_):
            raise ValueError("Number of input features has changed from %i "
                             "to %i between calls to partial_fit! Try "
                             "setting n_components to a fixed value." % (
                                 self.components_.shape[0], self.n_components_))

        if self.components_ is None:
            # This is the first pass through partial_fit
            self.n_samples_seen_ = 0
            col_var = X.var(axis=0)
            col_mean = X.mean(axis=0)
            X -= col_mean
            U, S, V = linalg.svd(X, full_matrices=False)
            U, V = svd_flip(U, V, u_based_decision=False)
            explained_variance = S ** 2 / n_samples
            explained_variance_ratio = S ** 2 / np.sum(col_var *
                                                       n_samples)
        else:
            col_batch_mean = X.mean(axis=0)
            col_mean, col_var, n_total_samples = _batch_mean_variance_update(
                X, self.mean_, self.var_, self.n_samples_seen_)
            X -= col_batch_mean
            # Build matrix of combined previous basis and new data
            mean_correction = np.sqrt((self.n_samples_seen_ * n_samples) /
                                      n_total_samples) * (self.mean_ -
                                                          col_batch_mean)
            X_combined = np.vstack((self.singular_values_.reshape((-1, 1)) *
                                    self.components_, X,
                                    mean_correction))
            U, S, V = linalg.svd(X_combined, full_matrices=False)
            U, V = svd_flip(U, V, u_based_decision=False)
            explained_variance = S ** 2 / n_total_samples
            explained_variance_ratio = S ** 2 / np.sum(col_var *
                                                       n_total_samples)
        self.n_samples_seen_ += n_samples
        self.components_ = V[:self.n_components_]
        self.singular_values_ = S[:self.n_components_]
        self.mean_ = col_mean
        self.var_ = col_var
        self.explained_variance_ = explained_variance[:self.n_components_]
        self.explained_variance_ratio_ = \
            explained_variance_ratio[:self.n_components_]
        if self.n_components_ < n_features:
            self.noise_variance_ = \
                explained_variance[self.n_components_:].mean()
        else:
            self.noise_variance_ = 0.
        return self
