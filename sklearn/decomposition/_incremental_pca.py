"""Incremental Principal Components Analysis."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Integral

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import svds

from sklearn.utils import metadata_routing

from ..base import _fit_context
from ..utils import gen_batches
from ..utils._arpack import _init_arpack_v0
from ..utils._param_validation import Interval, StrOptions
from ..utils.extmath import _incremental_mean_and_var, svd_flip
from ..utils.sparsefuncs import _implicit_column_offset, _implicit_vstack
from ..utils.sparsefuncs import (
    incr_mean_variance_axis as _sparse_incremental_mean_and_var,
)
from ..utils.validation import check_random_state, validate_data
from ._base import _BasePCA


class IncrementalPCA(_BasePCA):
    """Incremental principal components analysis (IPCA).

    Linear dimensionality reduction using Singular Value Decomposition of
    the data, keeping only the most significant singular vectors to
    project the data to a lower dimensional space. The input data is centered
    but not scaled for each feature before applying the SVD.

    Depending on the size of the input data, this algorithm can be much more
    memory efficient than a PCA.

    This algorithm has constant memory complexity, on the order of
    `batch_size * n_features`, enabling use of :class:`numpy.memmap` files without
    loading the entire file into memory. For sparse inputs, it can use the ARPACK
    implementation of the truncated SVD via :func:`scipy.sparse.linalg.svds`.

    The computational overhead of each SVD is
    ``O(batch_size * n_features ** 2)``, but only 2 * batch_size samples
    remain in memory at a time. There will be ``n_samples / batch_size`` SVD
    computations to get the principal components, versus 1 large SVD of
    complexity ``O(n_samples * n_features ** 2)`` for PCA.

    For a usage example, see
    :ref:`sphx_glr_auto_examples_decomposition_plot_incremental_pca.py`.

    Read more in the :ref:`User Guide <IncrementalPCA>`.

    .. versionadded:: 0.16

    Parameters
    ----------
    n_components : int, default=None
        Number of components to keep. If `n_components=None`, it is set to
        `min(n_samples, n_features) - 1` if the `"arpack"` solver is used and
        `min(n_samples, n_features)` otherwise.

    whiten : bool, default=False
        When True (False by default) the ``components_`` vectors are divided
        by ``n_samples`` times ``components_`` to ensure uncorrelated outputs
        with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometimes
        improve the predictive accuracy of the downstream estimators by
        making data respect some hard-wired assumptions.

    svd_solver : {"full", "arpack"}, default="full"
        The SVD solver to use.

        - `"full"`: Run exact full SVD calling the standard LAPACK solver via
          :func:`scipy.linalg.svd` and select the components by postprocessing. For
          sparse inputs, this solver will densify data in batches first.
        - `"arpack"`: Run SVD truncated to `n_components` calling ARPACK solver via
          :func:`scipy.sparse.linalg.svds`. It requires strictly
          `0 < n_components < min(X.shape)`. This solver is recommended for sparse
          inputs to avoid densifying the input data.

        .. versionadded:: 1.8

    copy : bool, default=True
        If False, X will be overwritten. ``copy=False`` can be used to
        save memory but is unsafe for general use.

    batch_size : int, default=None
        The number of samples to use for each batch. Only used when calling
        ``fit``. If ``batch_size`` is ``None``, then ``batch_size``
        is inferred from the data and set to ``5 * n_features``, to provide a
        balance between approximation accuracy and memory consumption.

    random_state : int, RandomState instance or None, default=None
        Used when the `"arpack"` solver is used. Pass an int for reproducible results
        across multiple function calls. See :term:`Glossary <random_state>`.

        .. versionadded:: 1.8

    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Principal axes in feature space, representing the directions of
        maximum variance in the data. Equivalently, the right singular
        vectors of the centered input data, parallel to its eigenvectors.
        The components are sorted by decreasing ``explained_variance_``.

    explained_variance_ : ndarray of shape (n_components,)
        Variance explained by each of the selected components.

    explained_variance_ratio_ : ndarray of shape (n_components,)
        Percentage of variance explained by each of the selected components.
        If all components are stored, the sum of explained variances is equal
        to 1.0.

    singular_values_ : ndarray of shape (n_components,)
        The singular values corresponding to each of the selected components.
        The singular values are equal to the 2-norms of the ``n_components``
        variables in the lower-dimensional space.

    mean_ : ndarray of shape (n_features,)
        Per-feature empirical mean, aggregate over calls to ``partial_fit``.

    var_ : ndarray of shape (n_features,)
        Per-feature empirical variance, aggregate over calls to
        ``partial_fit``.

    noise_variance_ : float
        The estimated noise covariance following the Probabilistic PCA model
        from Tipping and Bishop 1999. See "Pattern Recognition and
        Machine Learning" by C. Bishop, 12.2.1 p. 574 or
        http://www.miketipping.com/papers/met-mppca.pdf.

    n_components_ : int
        The estimated number of components. Relevant when
        ``n_components=None``.

    n_samples_seen_ : int
        The number of samples processed by the estimator. Will be reset on
        new calls to fit, but increments across ``partial_fit`` calls.

    batch_size_ : int
        Inferred batch size from ``batch_size``.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    PCA : Principal component analysis (PCA).
    KernelPCA : Kernel Principal component analysis (KPCA).
    SparsePCA : Sparse Principal Components Analysis (SparsePCA).
    TruncatedSVD : Dimensionality reduction using truncated SVD.

    Notes
    -----
    Implements the incremental PCA model from:
    *D. Ross, J. Lim, R. Lin, M. Yang, Incremental Learning for Robust Visual
    Tracking, International Journal of Computer Vision, Volume 77, Issue 1-3,
    pp. 125-141, May 2008.*
    See https://www.cs.toronto.edu/~dross/ivt/RossLimLinYang_ijcv.pdf

    This model is an extension of the Sequential Karhunen-Loeve Transform from:
    :doi:`A. Levy and M. Lindenbaum, Sequential Karhunen-Loeve Basis Extraction and
    its Application to Images, IEEE Transactions on Image Processing, Volume 9,
    Number 8, pp. 1371-1374, August 2000. <10.1109/83.855432>`

    We have specifically abstained from an optimization used by authors of both
    papers, a QR decomposition used in specific situations to reduce the
    algorithmic complexity of the SVD. The source for this technique is
    *Matrix Computations, Third Edition, G. Holub and C. Van Loan, Chapter 5,
    section 5.4.4, pp 252-253.*. This technique has been omitted because it is
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

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.decomposition import IncrementalPCA
    >>> from scipy import sparse
    >>> X, _ = load_digits(return_X_y=True)
    >>> transformer = IncrementalPCA(n_components=7, batch_size=200)
    >>> # either partially fit on smaller batches of data
    >>> transformer.partial_fit(X[:100, :])
    IncrementalPCA(batch_size=200, n_components=7)
    >>> # or let the fit function itself divide the data into batches
    >>> X_sparse = sparse.csr_matrix(X)
    >>> X_transformed = transformer.fit_transform(X_sparse)
    >>> X_transformed.shape
    (1797, 7)
    """

    __metadata_request__partial_fit = {"check_input": metadata_routing.UNUSED}

    _parameter_constraints: dict = {
        "n_components": [Interval(Integral, 1, None, closed="left"), None],
        "whiten": ["boolean"],
        "svd_solver": [StrOptions({"full", "arpack"})],
        "copy": ["boolean"],
        "batch_size": [Interval(Integral, 1, None, closed="left"), None],
        "random_state": ["random_state"],
    }

    def __init__(
        self,
        n_components=None,
        *,
        whiten=False,
        svd_solver="full",
        copy=True,
        batch_size=None,
        random_state=None,
    ):
        self.n_components = n_components
        self.whiten = whiten
        self.svd_solver = svd_solver
        self.copy = copy
        self.batch_size = batch_size
        self.random_state = random_state

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y=None):
        """Fit the model with X, using minibatches of size batch_size.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self.components_ = None
        self.n_samples_seen_ = 0
        self.mean_ = 0.0
        self.var_ = 0.0
        self.singular_values_ = None
        self.explained_variance_ = None
        self.explained_variance_ratio_ = None
        self.noise_variance_ = None

        X = validate_data(
            self,
            X,
            accept_sparse=["csr", "csc"],
            copy=self.copy,
            dtype=[np.float64, np.float32],
            force_writeable=True,
        )
        n_samples, n_features = X.shape

        if self.batch_size is None:
            self.batch_size_ = 5 * n_features
        else:
            self.batch_size_ = self.batch_size

        for batch in gen_batches(
            n_samples, self.batch_size_, min_batch_size=self.n_components or 0
        ):
            X_batch = X[batch]
            if sparse.issparse(X_batch) and self.svd_solver != "arpack":
                X_batch = X_batch.toarray()
            self.partial_fit(X_batch, check_input=False)

        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def partial_fit(self, X, y=None, check_input=True):
        """Incremental fit with X. All of X is processed as a single batch.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        check_input : bool, default=True
            Run check_array on X.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        first_pass = not hasattr(self, "components_")

        if check_input:
            X = validate_data(
                self,
                X,
                accept_sparse=["csr", "csc"],
                copy=self.copy,
                dtype=[np.float64, np.float32],
                force_writeable=True,
                reset=first_pass,
            )
        n_samples, n_features = X.shape
        X_is_sparse = sparse.issparse(X)

        if X_is_sparse and self.svd_solver != "arpack":
            raise TypeError(
                "IncrementalPCA.partial_fit only support sparse inputs with the "
                f'"arpack" solver, while "{self.svd_solver}" was passed. You may '
                'consider using the "arpack" solver, or converting data to dense, '
                "or using IncrementalPCA.fit to densify data in batches."
            )

        if first_pass:
            self.components_ = None

        if self.n_components is None:
            if self.components_ is None:
                if self.svd_solver == "arpack":
                    self.n_components_ = min(n_samples, n_features) - 1
                else:
                    self.n_components_ = min(n_samples, n_features)
            else:
                self.n_components_ = self.components_.shape[0]
        elif not self.n_components <= n_features:
            raise ValueError(
                "n_components=%r invalid for n_features=%d, need "
                "more rows than columns for IncrementalPCA "
                "processing" % (self.n_components, n_features)
            )
        elif self.n_components > n_samples and first_pass:
            raise ValueError(
                f"n_components={self.n_components} must be less or equal to "
                f"the batch number of samples {n_samples} for the first "
                "partial_fit call."
            )
        elif self.svd_solver == "arpack" and self.n_components == min(
            n_samples, n_features
        ):
            raise ValueError(
                f"n_components={self.n_components} must be strictly less than "
                f'{min(n_samples, n_features)=} with the "arpack" solver.'
            )
        else:
            self.n_components_ = self.n_components

        if (self.components_ is not None) and (
            self.components_.shape[0] != self.n_components_
        ):
            raise ValueError(
                "Number of input features has changed from %i "
                "to %i between calls to partial_fit! Try "
                "setting n_components to a fixed value."
                % (self.components_.shape[0], self.n_components_)
            )

        # This is the first partial_fit
        if not hasattr(self, "n_samples_seen_"):
            self.n_samples_seen_ = 0
            self.mean_ = 0.0
            self.var_ = 0.0

        # Update stats - they are 0 if this is the first step
        if X_is_sparse:
            # _sparse_incremental_mean_and_var only accepts `last_mean` and `last_var`
            # as arrays; moreover it modifies them in-place so we need to make a copy
            # since further computations still need the original values
            if np.size(self.mean_) == 1:
                last_mean = np.full(n_features, self.mean_)
            else:
                last_mean = self.mean_.copy()
            if np.size(self.var_) == 1:
                last_var = np.full(n_features, self.var_)
            else:
                last_var = self.var_.copy()

            col_mean, col_var, n_total_samples = _sparse_incremental_mean_and_var(
                X,
                axis=0,
                last_mean=last_mean,
                last_var=last_var,
                last_n=self.n_samples_seen_,
            )
        else:
            col_mean, col_var, n_total_samples = _incremental_mean_and_var(
                X,
                last_mean=self.mean_,
                last_variance=self.var_,
                last_sample_count=np.repeat(self.n_samples_seen_, n_features),
            )
        n_total_samples = n_total_samples[0]

        # Whitening
        if self.n_samples_seen_ == 0:
            # If it is the first step, simply whiten X
            if X_is_sparse:
                # Avoid densifying sparse data
                X = _implicit_column_offset(X, col_mean)
            else:
                X -= col_mean
        else:
            col_batch_mean = np.asarray(np.mean(X, axis=0))
            if X_is_sparse:
                # Avoid densifying sparse data
                X = _implicit_column_offset(X, col_batch_mean)
            else:
                X -= col_batch_mean

            # Build matrix of combined previous basis and new data
            mean_correction = np.sqrt(
                (self.n_samples_seen_ / n_total_samples) * n_samples
            ) * (self.mean_ - col_batch_mean)

            vstack_func = _implicit_vstack if X_is_sparse else np.vstack
            X = vstack_func(
                (
                    self.singular_values_.reshape((-1, 1)) * self.components_,
                    X,
                    np.atleast_2d(mean_correction),
                )
            )

        if self.svd_solver == "arpack":
            random_state = check_random_state(self.random_state)
            v0 = _init_arpack_v0(min(X.shape), random_state=random_state)
            U, S, Vt = svds(X, k=self.n_components_, v0=v0)
            S = S[::-1]
            _, Vt = svd_flip(U[:, ::-1], Vt[::-1], u_based_decision=False)
        else:  # self.svd_solver == "full"
            U, S, Vt = linalg.svd(X, full_matrices=False, check_finite=False)
            _, Vt = svd_flip(U, Vt, u_based_decision=False)

        explained_variance = S**2 / (n_total_samples - 1)
        explained_variance_ratio = S**2 / np.sum(col_var * n_total_samples)

        self.n_samples_seen_ = n_total_samples
        self.mean_ = col_mean
        self.var_ = col_var

        if self.svd_solver == "arpack":
            self.components_ = Vt
            self.singular_values_ = S
            self.explained_variance_ = explained_variance
            self.explained_variance_ratio_ = explained_variance_ratio
        else:
            self.components_ = Vt[: self.n_components_]
            self.singular_values_ = S[: self.n_components_]
            self.explained_variance_ = explained_variance[: self.n_components_]
            self.explained_variance_ratio_ = explained_variance_ratio[
                : self.n_components_
            ]

        if self.n_components_ < min(n_features, n_samples):
            if self.svd_solver == "arpack":
                total_var = self.var_.sum() * n_samples / (n_samples - 1)  # ddof=1
                self.noise_variance_ = total_var - np.sum(self.explained_variance_)
                self.noise_variance_ /= min(n_features, n_samples) - self.n_components_
            else:
                self.noise_variance_ = explained_variance[self.n_components_ :].mean()
        else:
            self.noise_variance_ = 0.0
        return self

    def transform(self, X):
        """Apply dimensionality reduction to X.

        X is projected on the first principal components previously extracted
        from a training set, using minibatches of size batch_size if X is
        sparse.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            Projection of X in the first principal components.

        Examples
        --------

        >>> import numpy as np
        >>> from sklearn.decomposition import IncrementalPCA
        >>> X = np.array([[-1, -1], [-2, -1], [-3, -2],
        ...               [1, 1], [2, 1], [3, 2]])
        >>> ipca = IncrementalPCA(n_components=2, batch_size=3)
        >>> ipca.fit(X)
        IncrementalPCA(batch_size=3, n_components=2)
        >>> ipca.transform(X) # doctest: +SKIP
        """
        if sparse.issparse(X):
            n_samples = X.shape[0]
            output = []
            for batch in gen_batches(
                n_samples, self.batch_size_, min_batch_size=self.n_components or 0
            ):
                output.append(super().transform(X[batch].toarray()))
            return np.vstack(output)
        else:
            return super().transform(X)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        # Beware that fit accepts sparse data but partial_fit doesn't
        tags.input_tags.sparse = True
        return tags
