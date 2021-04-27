"""Kernel Principal Components Analysis."""

# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: BSD 3 clause

import numpy as np
from scipy import linalg
from scipy.sparse.linalg import eigsh

from ..utils._arpack import _init_arpack_v0
from ..utils.extmath import svd_flip
from ..utils.validation import check_is_fitted, _check_psd_eigenvalues
from ..utils.deprecation import deprecated
from ..exceptions import NotFittedError
from ..base import BaseEstimator, TransformerMixin
from ..preprocessing import KernelCenterer
from ..metrics.pairwise import pairwise_kernels
from ..utils.validation import _deprecate_positional_args


class KernelPCA(TransformerMixin, BaseEstimator):
    """Kernel Principal component analysis (KPCA).

    Non-linear dimensionality reduction through the use of kernels (see
    :ref:`metrics`).

    Read more in the :ref:`User Guide <kernel_PCA>`.

    Parameters
    ----------
    n_components : int, default=None
        Number of components. If None, all non-zero components are kept.

    kernel : {'linear', 'poly', \
            'rbf', 'sigmoid', 'cosine', 'precomputed'}, default='linear'
        Kernel used for PCA.

    gamma : float, default=None
        Kernel coefficient for rbf, poly and sigmoid kernels. Ignored by other
        kernels. If ``gamma`` is ``None``, then it is set to ``1/n_features``.

    degree : int, default=3
        Degree for poly kernels. Ignored by other kernels.

    coef0 : float, default=1
        Independent term in poly and sigmoid kernels.
        Ignored by other kernels.

    kernel_params : dict, default=None
        Parameters (keyword arguments) and
        values for kernel passed as callable object.
        Ignored by other kernels.

    alpha : float, default=1.0
        Hyperparameter of the ridge regression that learns the
        inverse transform (when fit_inverse_transform=True).

    fit_inverse_transform : bool, default=False
        Learn the inverse transform for non-precomputed kernels.
        (i.e. learn to find the pre-image of a point)

    eigen_solver : {'auto', 'dense', 'arpack'}, default='auto'
        Select eigensolver to use. If n_components is much less than
        the number of training samples, arpack may be more efficient
        than the dense eigensolver.

    tol : float, default=0
        Convergence tolerance for arpack.
        If 0, optimal value will be chosen by arpack.

    max_iter : int, default=None
        Maximum number of iterations for arpack.
        If None, optimal value will be chosen by arpack.

    remove_zero_eig : bool, default=False
        If True, then all components with zero eigenvalues are removed, so
        that the number of components in the output may be < n_components
        (and sometimes even zero due to numerical instability).
        When n_components is None, this parameter is ignored and components
        with zero eigenvalues are removed regardless.

    random_state : int, RandomState instance or None, default=None
        Used when ``eigen_solver`` == 'arpack'. Pass an int for reproducible
        results across multiple function calls.
        See :term:`Glossary <random_state>`.

        .. versionadded:: 0.18

    copy_X : bool, default=True
        If True, input X is copied and stored by the model in the `X_fit_`
        attribute. If no further changes will be done to X, setting
        `copy_X=False` saves memory by storing a reference.

        .. versionadded:: 0.18

    n_jobs : int, default=None
        The number of parallel jobs to run.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

        .. versionadded:: 0.18

    Attributes
    ----------
    eigenvalues_ : ndarray of shape (n_components,)
        Eigenvalues of the centered kernel matrix in decreasing order.
        If `n_components` and `remove_zero_eig` are not set,
        then all values are stored.

    eigenvectors_ : ndarray of shape (n_samples, n_components)
        Eigenvectors of the centered kernel matrix. If `n_components` and
        `remove_zero_eig` are not set, then all components are stored.

    dual_coef_ : ndarray of shape (n_samples, n_features)
        Inverse transform matrix. Only available when
        ``fit_inverse_transform`` is True.

    X_transformed_fit_ : ndarray of shape (n_samples, n_components)
        Projection of the fitted data on the kernel principal components.
        Only available when ``fit_inverse_transform`` is True.

    X_fit_ : ndarray of shape (n_samples, n_features)
        The data used to fit the model. If `copy_X=False`, then `X_fit_` is
        a reference. This attribute is used for the calls to transform.

    .. deprecated:: 1.0
        ``lambdas_`` was renamed to ``eigenvalues_`` in version 1.0 and will be
        removed in 1.2.

    .. deprecated:: 1.0
        ``alphas_`` was renamed to ``eigenvectors_`` in version 1.0 and will be
        removed in 1.2.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.decomposition import KernelPCA
    >>> X, _ = load_digits(return_X_y=True)
    >>> transformer = KernelPCA(n_components=7, kernel='linear')
    >>> X_transformed = transformer.fit_transform(X)
    >>> X_transformed.shape
    (1797, 7)

    References
    ----------
    Kernel PCA was introduced in:
        Bernhard Schoelkopf, Alexander J. Smola,
        and Klaus-Robert Mueller. 1999. Kernel principal
        component analysis. In Advances in kernel methods,
        MIT Press, Cambridge, MA, USA 327-352.
    """
    @_deprecate_positional_args
    def __init__(self, n_components=None, *, kernel="linear",
                 gamma=None, degree=3, coef0=1, kernel_params=None,
                 alpha=1.0, fit_inverse_transform=False, eigen_solver='auto',
                 tol=0, max_iter=None, remove_zero_eig=False,
                 random_state=None, copy_X=True, n_jobs=None):
        if fit_inverse_transform and kernel == 'precomputed':
            raise ValueError(
                "Cannot fit_inverse_transform with a precomputed kernel.")
        self.n_components = n_components
        self.kernel = kernel
        self.kernel_params = kernel_params
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.alpha = alpha
        self.fit_inverse_transform = fit_inverse_transform
        self.eigen_solver = eigen_solver
        self.remove_zero_eig = remove_zero_eig
        self.tol = tol
        self.max_iter = max_iter
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.copy_X = copy_X

    # TODO: Remove in 1.1
    # mypy error: Decorated property not supported
    @deprecated("Attribute _pairwise was deprecated in "  # type: ignore
                "version 0.24 and will be removed in 1.1 (renaming of 0.26).")
    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def _get_kernel(self, X, Y=None):
        if callable(self.kernel):
            params = self.kernel_params or {}
        else:
            params = {"gamma": self.gamma,
                      "degree": self.degree,
                      "coef0": self.coef0}
        return pairwise_kernels(X, Y, metric=self.kernel,
                                filter_params=True, n_jobs=self.n_jobs,
                                **params)

    def _fit_transform(self, K):
        """ Fit's using kernel K"""
        # center kernel
        K = self._centerer.fit_transform(K)

        if self.n_components is None:
            n_components = K.shape[0]
        else:
            n_components = min(K.shape[0], self.n_components)

        # compute eigenvectors
        if self.eigen_solver == 'auto':
            if K.shape[0] > 200 and n_components < 10:
                eigen_solver = 'arpack'
            else:
                eigen_solver = 'dense'
        else:
            eigen_solver = self.eigen_solver

        if eigen_solver == 'dense':
            self.eigenvalues_, self.eigenvectors_ = linalg.eigh(
                K, eigvals=(K.shape[0] - n_components, K.shape[0] - 1))
        elif eigen_solver == 'arpack':
            v0 = _init_arpack_v0(K.shape[0], self.random_state)
            self.eigenvalues_, self.eigenvectors_ = \
                eigsh(K, n_components, which="LA", tol=self.tol,
                      maxiter=self.max_iter, v0=v0)

        # make sure that the eigenvalues are ok and fix numerical issues
        self.eigenvalues_ = _check_psd_eigenvalues(self.eigenvalues_,
                                                   enable_warnings=False)

        # flip eigenvectors' sign to enforce deterministic output
        self.eigenvectors_, _ = svd_flip(self.eigenvectors_,
                                         np.zeros_like(self.eigenvectors_).T)

        # sort eigenvectors in descending order
        indices = self.eigenvalues_.argsort()[::-1]
        self.eigenvalues_ = self.eigenvalues_[indices]
        self.eigenvectors_ = self.eigenvectors_[:, indices]

        # remove eigenvectors with a zero eigenvalue (null space) if required
        if self.remove_zero_eig or self.n_components is None:
            self.eigenvectors_ = self.eigenvectors_[:, self.eigenvalues_ > 0]
            self.eigenvalues_ = self.eigenvalues_[self.eigenvalues_ > 0]

        # Maintenance note on Eigenvectors normalization
        # ----------------------------------------------
        # there is a link between
        # the eigenvectors of K=Phi(X)'Phi(X) and the ones of Phi(X)Phi(X)'
        # if v is an eigenvector of K
        #     then Phi(X)v  is an eigenvector of Phi(X)Phi(X)'
        # if u is an eigenvector of Phi(X)Phi(X)'
        #     then Phi(X)'u is an eigenvector of Phi(X)'Phi(X)
        #
        # At this stage our self.eigenvectors_ (the v) have norm 1, we need to
        # scale them so that eigenvectors in kernel feature space (the u) have
        # norm=1 instead
        #
        # We COULD scale them here:
        #  self.eigenvectors_ = self.eigenvectors_ / np.sqrt(self.eigenvalues_)
        #
        # But choose to perform that LATER when needed, in `fit()` and in
        # `transform()`.

        return K

    def _fit_inverse_transform(self, X_transformed, X):
        if hasattr(X, "tocsr"):
            raise NotImplementedError("Inverse transform not implemented for "
                                      "sparse matrices!")

        n_samples = X_transformed.shape[0]
        K = self._get_kernel(X_transformed)
        K.flat[::n_samples + 1] += self.alpha
        self.dual_coef_ = linalg.solve(K, X, sym_pos=True, overwrite_a=True)
        self.X_transformed_fit_ = X_transformed

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        X = self._validate_data(X, accept_sparse='csr', copy=self.copy_X)
        self._centerer = KernelCenterer()
        K = self._get_kernel(X)
        self._fit_transform(K)

        if self.fit_inverse_transform:
            # no need to use the kernel to transform X, use shortcut expression
            X_transformed = self.eigenvectors_ * np.sqrt(self.eigenvalues_)

            self._fit_inverse_transform(X_transformed, X)

        self.X_fit_ = X
        return self

    def fit_transform(self, X, y=None, **params):
        """Fit the model from data in X and transform X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
        """
        self.fit(X, **params)

        # no need to use the kernel to transform X, use shortcut expression
        X_transformed = self.eigenvectors_ * np.sqrt(self.eigenvalues_)

        if self.fit_inverse_transform:
            self._fit_inverse_transform(X_transformed, X)

        return X_transformed

    def transform(self, X):
        """Transform X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
        """
        check_is_fitted(self)
        X = self._validate_data(X, accept_sparse='csr', reset=False)

        # Compute centered gram matrix between X and training data X_fit_
        K = self._centerer.transform(self._get_kernel(X, self.X_fit_))

        # scale eigenvectors (properly account for null-space for dot product)
        non_zeros = np.flatnonzero(self.eigenvalues_)
        scaled_eigenvectors = np.zeros_like(self.eigenvectors_)
        scaled_eigenvectors[:, non_zeros] = (self.eigenvectors_[:, non_zeros]
                                             / np.sqrt(self.eigenvalues_[non_zeros]))

        # Project with a scalar product between K and the scaled eigenvectors
        return np.dot(K, scaled_eigenvectors)

    def inverse_transform(self, X):
        """Transform X back to original space.

        ``inverse_transform`` approximates the inverse transformation using
        a learned pre-image. The pre-image is learned by kernel ridge
        regression of the original data on their low-dimensional representation
        vectors.

        .. note:
            :meth:`~sklearn.decomposition.fit` internally uses a centered
            kernel. As the centered kernel no longer contains the information
            of the mean of kernel features, such information is not taken into
            account in reconstruction.

        .. note::
            When users want to compute inverse transformation for 'linear'
            kernel, it is recommended that they use
            :class:`~sklearn.decomposition.PCA` instead. Unlike
            :class:`~sklearn.decomposition.PCA`,
            :class:`~sklearn.decomposition.KernelPCA`'s ``inverse_transform``
            does not reconstruct the mean of data when 'linear' kernel is used
            due to the use of centered kernel.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_components)

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_features)

        References
        ----------
        "Learning to Find Pre-Images", G BakIr et al, 2004.
        """
        if not self.fit_inverse_transform:
            raise NotFittedError("The fit_inverse_transform parameter was not"
                                 " set to True when instantiating and hence "
                                 "the inverse transform is not available.")

        K = self._get_kernel(X, self.X_transformed_fit_)
        return np.dot(K, self.dual_coef_)

    def _more_tags(self):
        return {'preserves_dtype': [np.float64, np.float32],
                'pairwise': self.kernel == 'precomputed'}

    # TODO: Remove in 1.2
    # mypy error: Decorated property not supported
    @deprecated("Attribute 'lambdas_' was deprecated in "  # type: ignore
                "version 1.0 and will be removed in 1.2. Use "
                "'eigenvalues_' instead")
    @property
    def lambdas_(self):
        return self.eigenvalues_

    # mypy error: Decorated property not supported
    @deprecated("Attribute 'alphas_' was deprecated in "  # type: ignore
                "version 1.0 and will be removed in 1.2. Use "
                "'eigenvectors_' instead")
    @property
    def alphas_(self):
        return self.eigenvectors_
