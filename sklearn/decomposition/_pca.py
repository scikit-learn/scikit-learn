"""Principal Component Analysis."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from math import log, sqrt
from numbers import Integral, Real

import numpy as np
from scipy import linalg
from scipy.sparse import issparse
from scipy.sparse.linalg import svds
from scipy.special import gammaln

from ..base import _fit_context
from ..utils import check_random_state
from ..utils._arpack import _init_arpack_v0
from ..utils._array_api import _convert_to_numpy, get_namespace
from ..utils._param_validation import Interval, RealNotInt, StrOptions
from ..utils.extmath import fast_logdet, randomized_svd, stable_cumsum, svd_flip
from ..utils.sparsefuncs import _implicit_column_offset, mean_variance_axis
from ..utils.validation import check_is_fitted, validate_data
from ._base import _BasePCA


def _assess_dimension(spectrum, rank, n_samples):
    """Compute the log-likelihood of a rank ``rank`` dataset.

    The dataset is assumed to be embedded in gaussian noise of shape(n,
    dimf) having spectrum ``spectrum``. This implements the method of
    T. P. Minka.

    Parameters
    ----------
    spectrum : ndarray of shape (n_features,)
        Data spectrum.
    rank : int
        Tested rank value. It should be strictly lower than n_features,
        otherwise the method isn't specified (division by zero in equation
        (31) from the paper).
    n_samples : int
        Number of samples.

    Returns
    -------
    ll : float
        The log-likelihood.

    References
    ----------
    This implements the method of `Thomas P. Minka:
    Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
    <https://proceedings.neurips.cc/paper/2000/file/7503cfacd12053d309b6bed5c89de212-Paper.pdf>`_
    """
    xp, _ = get_namespace(spectrum)

    n_features = spectrum.shape[0]
    if not 1 <= rank < n_features:
        raise ValueError("the tested rank should be in [1, n_features - 1]")

    eps = 1e-15

    if spectrum[rank - 1] < eps:
        # When the tested rank is associated with a small eigenvalue, there's
        # no point in computing the log-likelihood: it's going to be very
        # small and won't be the max anyway. Also, it can lead to numerical
        # issues below when computing pa, in particular in log((spectrum[i] -
        # spectrum[j]) because this will take the log of something very small.
        return -xp.inf

    pu = -rank * log(2.0)
    for i in range(1, rank + 1):
        pu += (
            gammaln((n_features - i + 1) / 2.0)
            - log(xp.pi) * (n_features - i + 1) / 2.0
        )

    pl = xp.sum(xp.log(spectrum[:rank]))
    pl = -pl * n_samples / 2.0

    v = max(eps, xp.sum(spectrum[rank:]) / (n_features - rank))
    pv = -log(v) * n_samples * (n_features - rank) / 2.0

    m = n_features * rank - rank * (rank + 1.0) / 2.0
    pp = log(2.0 * xp.pi) * (m + rank) / 2.0

    pa = 0.0
    spectrum_ = xp.asarray(spectrum, copy=True)
    spectrum_[rank:n_features] = v
    for i in range(rank):
        for j in range(i + 1, spectrum.shape[0]):
            pa += log(
                (spectrum[i] - spectrum[j]) * (1.0 / spectrum_[j] - 1.0 / spectrum_[i])
            ) + log(n_samples)

    ll = pu + pl + pv + pp - pa / 2.0 - rank * log(n_samples) / 2.0

    return ll


def _infer_dimension(spectrum, n_samples):
    """Infers the dimension of a dataset with a given spectrum.

    The returned value will be in [1, n_features - 1].
    """
    xp, _ = get_namespace(spectrum)

    ll = xp.empty_like(spectrum)
    ll[0] = -xp.inf  # we don't want to return n_components = 0
    for rank in range(1, spectrum.shape[0]):
        ll[rank] = _assess_dimension(spectrum, rank, n_samples)
    return xp.argmax(ll)


class PCA(_BasePCA):
    """Principal component analysis (PCA).

    Linear dimensionality reduction using Singular Value Decomposition of the
    data to project it to a lower dimensional space. The input data is centered
    but not scaled for each feature before applying the SVD.

    It uses the LAPACK implementation of the full SVD or a randomized truncated
    SVD by the method of Halko et al. 2009, depending on the shape of the input
    data and the number of components to extract.

    With sparse inputs, the ARPACK implementation of the truncated SVD can be
    used (i.e. through :func:`scipy.sparse.linalg.svds`). Alternatively, one
    may consider :class:`TruncatedSVD` where the data are not centered.

    Notice that this class only supports sparse inputs for some solvers such as
    "arpack" and "covariance_eigh". See :class:`TruncatedSVD` for an
    alternative with sparse data.

    For a usage example, see
    :ref:`sphx_glr_auto_examples_decomposition_plot_pca_iris.py`

    Read more in the :ref:`User Guide <PCA>`.

    Parameters
    ----------
    n_components : int, float or 'mle', default=None
        Number of components to keep.
        if n_components is not set all components are kept::

            n_components == min(n_samples, n_features)

        If ``n_components == 'mle'`` and ``svd_solver == 'full'``, Minka's
        MLE is used to guess the dimension. Use of ``n_components == 'mle'``
        will interpret ``svd_solver == 'auto'`` as ``svd_solver == 'full'``.

        If ``0 < n_components < 1`` and ``svd_solver == 'full'``, select the
        number of components such that the amount of variance that needs to be
        explained is greater than the percentage specified by n_components.

        If ``svd_solver == 'arpack'``, the number of components must be
        strictly less than the minimum of n_features and n_samples.

        Hence, the None case results in::

            n_components == min(n_samples, n_features) - 1

    copy : bool, default=True
        If False, data passed to fit are overwritten and running
        fit(X).transform(X) will not yield the expected results,
        use fit_transform(X) instead.

    whiten : bool, default=False
        When True (False by default) the `components_` vectors are multiplied
        by the square root of n_samples and then divided by the singular values
        to ensure uncorrelated outputs with unit component-wise variances.

        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making their data respect some hard-wired assumptions.

    svd_solver : {'auto', 'full', 'covariance_eigh', 'arpack', 'randomized'},\
            default='auto'
        "auto" :
            The solver is selected by a default 'auto' policy is based on `X.shape` and
            `n_components`: if the input data has fewer than 1000 features and
            more than 10 times as many samples, then the "covariance_eigh"
            solver is used. Otherwise, if the input data is larger than 500x500
            and the number of components to extract is lower than 80% of the
            smallest dimension of the data, then the more efficient
            "randomized" method is selected. Otherwise the exact "full" SVD is
            computed and optionally truncated afterwards.
        "full" :
            Run exact full SVD calling the standard LAPACK solver via
            `scipy.linalg.svd` and select the components by postprocessing
        "covariance_eigh" :
            Precompute the covariance matrix (on centered data), run a
            classical eigenvalue decomposition on the covariance matrix
            typically using LAPACK and select the components by postprocessing.
            This solver is very efficient for n_samples >> n_features and small
            n_features. It is, however, not tractable otherwise for large
            n_features (large memory footprint required to materialize the
            covariance matrix). Also note that compared to the "full" solver,
            this solver effectively doubles the condition number and is
            therefore less numerical stable (e.g. on input data with a large
            range of singular values).
        "arpack" :
            Run SVD truncated to `n_components` calling ARPACK solver via
            `scipy.sparse.linalg.svds`. It requires strictly
            `0 < n_components < min(X.shape)`
        "randomized" :
            Run randomized SVD by the method of Halko et al.

        .. versionadded:: 0.18.0

        .. versionchanged:: 1.5
            Added the 'covariance_eigh' solver.

    tol : float, default=0.0
        Tolerance for singular values computed by svd_solver == 'arpack'.
        Must be of range [0.0, infinity).

        .. versionadded:: 0.18.0

    iterated_power : int or 'auto', default='auto'
        Number of iterations for the power method computed by
        svd_solver == 'randomized'.
        Must be of range [0, infinity).

        .. versionadded:: 0.18.0

    n_oversamples : int, default=10
        This parameter is only relevant when `svd_solver="randomized"`.
        It corresponds to the additional number of random vectors to sample the
        range of `X` so as to ensure proper conditioning. See
        :func:`~sklearn.utils.extmath.randomized_svd` for more details.

        .. versionadded:: 1.1

    power_iteration_normalizer : {'auto', 'QR', 'LU', 'none'}, default='auto'
        Power iteration normalizer for randomized SVD solver.
        Not used by ARPACK. See :func:`~sklearn.utils.extmath.randomized_svd`
        for more details.

        .. versionadded:: 1.1

    random_state : int, RandomState instance or None, default=None
        Used when the 'arpack' or 'randomized' solvers are used. Pass an int
        for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

        .. versionadded:: 0.18.0

    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Principal axes in feature space, representing the directions of
        maximum variance in the data. Equivalently, the right singular
        vectors of the centered input data, parallel to its eigenvectors.
        The components are sorted by decreasing ``explained_variance_``.

    explained_variance_ : ndarray of shape (n_components,)
        The amount of variance explained by each of the selected components.
        The variance estimation uses `n_samples - 1` degrees of freedom.

        Equal to n_components largest eigenvalues
        of the covariance matrix of X.

        .. versionadded:: 0.18

    explained_variance_ratio_ : ndarray of shape (n_components,)
        Percentage of variance explained by each of the selected components.

        If ``n_components`` is not set then all components are stored and the
        sum of the ratios is equal to 1.0.

    singular_values_ : ndarray of shape (n_components,)
        The singular values corresponding to each of the selected components.
        The singular values are equal to the 2-norms of the ``n_components``
        variables in the lower-dimensional space.

        .. versionadded:: 0.19

    mean_ : ndarray of shape (n_features,)
        Per-feature empirical mean, estimated from the training set.

        Equal to `X.mean(axis=0)`.

    n_components_ : int
        The estimated number of components. When n_components is set
        to 'mle' or a number between 0 and 1 (with svd_solver == 'full') this
        number is estimated from input data. Otherwise it equals the parameter
        n_components, or the lesser value of n_features and n_samples
        if n_components is None.

    n_samples_ : int
        Number of samples in the training data.

    noise_variance_ : float
        The estimated noise covariance following the Probabilistic PCA model
        from Tipping and Bishop 1999. See "Pattern Recognition and
        Machine Learning" by C. Bishop, 12.2.1 p. 574 or
        http://www.miketipping.com/papers/met-mppca.pdf. It is required to
        compute the estimated data covariance and score samples.

        Equal to the average of (min(n_features, n_samples) - n_components)
        smallest eigenvalues of the covariance matrix of X.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    KernelPCA : Kernel Principal Component Analysis.
    SparsePCA : Sparse Principal Component Analysis.
    TruncatedSVD : Dimensionality reduction using truncated SVD.
    IncrementalPCA : Incremental Principal Component Analysis.

    References
    ----------
    For n_components == 'mle', this class uses the method from:
    `Minka, T. P.. "Automatic choice of dimensionality for PCA".
    In NIPS, pp. 598-604 <https://tminka.github.io/papers/pca/minka-pca.pdf>`_

    Implements the probabilistic PCA model from:
    `Tipping, M. E., and Bishop, C. M. (1999). "Probabilistic principal
    component analysis". Journal of the Royal Statistical Society:
    Series B (Statistical Methodology), 61(3), 611-622.
    <http://www.miketipping.com/papers/met-mppca.pdf>`_
    via the score and score_samples methods.

    For svd_solver == 'arpack', refer to `scipy.sparse.linalg.svds`.

    For svd_solver == 'randomized', see:
    :doi:`Halko, N., Martinsson, P. G., and Tropp, J. A. (2011).
    "Finding structure with randomness: Probabilistic algorithms for
    constructing approximate matrix decompositions".
    SIAM review, 53(2), 217-288.
    <10.1137/090771806>`
    and also
    :doi:`Martinsson, P. G., Rokhlin, V., and Tygert, M. (2011).
    "A randomized algorithm for the decomposition of matrices".
    Applied and Computational Harmonic Analysis, 30(1), 47-68.
    <10.1016/j.acha.2010.02.003>`

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.decomposition import PCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = PCA(n_components=2)
    >>> pca.fit(X)
    PCA(n_components=2)
    >>> print(pca.explained_variance_ratio_)
    [0.9924... 0.0075...]
    >>> print(pca.singular_values_)
    [6.30061... 0.54980...]

    >>> pca = PCA(n_components=2, svd_solver='full')
    >>> pca.fit(X)
    PCA(n_components=2, svd_solver='full')
    >>> print(pca.explained_variance_ratio_)
    [0.9924... 0.00755...]
    >>> print(pca.singular_values_)
    [6.30061... 0.54980...]

    >>> pca = PCA(n_components=1, svd_solver='arpack')
    >>> pca.fit(X)
    PCA(n_components=1, svd_solver='arpack')
    >>> print(pca.explained_variance_ratio_)
    [0.99244...]
    >>> print(pca.singular_values_)
    [6.30061...]
    """

    _parameter_constraints: dict = {
        "n_components": [
            Interval(Integral, 0, None, closed="left"),
            Interval(RealNotInt, 0, 1, closed="neither"),
            StrOptions({"mle"}),
            None,
        ],
        "copy": ["boolean"],
        "whiten": ["boolean"],
        "svd_solver": [
            StrOptions({"auto", "full", "covariance_eigh", "arpack", "randomized"})
        ],
        "tol": [Interval(Real, 0, None, closed="left")],
        "iterated_power": [
            StrOptions({"auto"}),
            Interval(Integral, 0, None, closed="left"),
        ],
        "n_oversamples": [Interval(Integral, 1, None, closed="left")],
        "power_iteration_normalizer": [StrOptions({"auto", "QR", "LU", "none"})],
        "random_state": ["random_state"],
    }

    def __init__(
        self,
        n_components=None,
        *,
        copy=True,
        whiten=False,
        svd_solver="auto",
        tol=0.0,
        iterated_power="auto",
        n_oversamples=10,
        power_iteration_normalizer="auto",
        random_state=None,
    ):
        self.n_components = n_components
        self.copy = copy
        self.whiten = whiten
        self.svd_solver = svd_solver
        self.tol = tol
        self.iterated_power = iterated_power
        self.n_oversamples = n_oversamples
        self.power_iteration_normalizer = power_iteration_normalizer
        self.random_state = random_state

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y=None):
        """Fit the model with X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Ignored.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._fit(X)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, X, y=None):
        """Fit the model with X and apply the dimensionality reduction on X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Ignored.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            Transformed values.

        Notes
        -----
        This method returns a Fortran-ordered array. To convert it to a
        C-ordered array, use 'np.ascontiguousarray'.
        """
        U, S, _, X, x_is_centered, xp = self._fit(X)
        if U is not None:
            U = U[:, : self.n_components_]

            if self.whiten:
                # X_new = X * V / S * sqrt(n_samples) = U * sqrt(n_samples)
                U *= sqrt(X.shape[0] - 1)
            else:
                # X_new = X * V = U * S * Vt * V = U * S
                U *= S[: self.n_components_]

            return U
        else:  # solver="covariance_eigh" does not compute U at fit time.
            return self._transform(X, xp, x_is_centered=x_is_centered)

    def _fit(self, X):
        """Dispatch to the right submethod depending on the chosen solver."""
        xp, is_array_api_compliant = get_namespace(X)

        # Raise an error for sparse input and unsupported svd_solver
        if issparse(X) and self.svd_solver not in ["auto", "arpack", "covariance_eigh"]:
            raise TypeError(
                'PCA only support sparse inputs with the "arpack" and'
                f' "covariance_eigh" solvers, while "{self.svd_solver}" was passed. See'
                " TruncatedSVD for a possible alternative."
            )
        if self.svd_solver == "arpack" and is_array_api_compliant:
            raise ValueError(
                "PCA with svd_solver='arpack' is not supported for Array API inputs."
            )

        # Validate the data, without ever forcing a copy as any solver that
        # supports sparse input data and the `covariance_eigh` solver are
        # written in a way to avoid the need for any inplace modification of
        # the input data contrary to the other solvers.
        # The copy will happen
        # later, only if needed, once the solver negotiation below is done.
        X = validate_data(
            self,
            X,
            dtype=[xp.float64, xp.float32],
            force_writeable=True,
            accept_sparse=("csr", "csc"),
            ensure_2d=True,
            copy=False,
        )
        self._fit_svd_solver = self.svd_solver
        if self._fit_svd_solver == "auto" and issparse(X):
            self._fit_svd_solver = "arpack"

        if self.n_components is None:
            if self._fit_svd_solver != "arpack":
                n_components = min(X.shape)
            else:
                n_components = min(X.shape) - 1
        else:
            n_components = self.n_components

        if self._fit_svd_solver == "auto":
            # Tall and skinny problems are best handled by precomputing the
            # covariance matrix.
            if X.shape[1] <= 1_000 and X.shape[0] >= 10 * X.shape[1]:
                self._fit_svd_solver = "covariance_eigh"
            # Small problem or n_components == 'mle', just call full PCA
            elif max(X.shape) <= 500 or n_components == "mle":
                self._fit_svd_solver = "full"
            elif 1 <= n_components < 0.8 * min(X.shape):
                self._fit_svd_solver = "randomized"
            # This is also the case of n_components in (0, 1)
            else:
                self._fit_svd_solver = "full"

        # Call different fits for either full or truncated SVD
        if self._fit_svd_solver in ("full", "covariance_eigh"):
            return self._fit_full(X, n_components, xp, is_array_api_compliant)
        elif self._fit_svd_solver in ["arpack", "randomized"]:
            return self._fit_truncated(X, n_components, xp)

    def _fit_full(self, X, n_components, xp, is_array_api_compliant):
        """Fit the model by computing full SVD on X."""
        n_samples, n_features = X.shape

        if n_components == "mle":
            if n_samples < n_features:
                raise ValueError(
                    "n_components='mle' is only supported if n_samples >= n_features"
                )
        elif not 0 <= n_components <= min(n_samples, n_features):
            raise ValueError(
                f"n_components={n_components} must be between 0 and "
                f"min(n_samples, n_features)={min(n_samples, n_features)} with "
                f"svd_solver={self._fit_svd_solver!r}"
            )

        self.mean_ = xp.mean(X, axis=0)
        # When X is a scipy sparse matrix, self.mean_ is a numpy matrix, so we need
        # to transform it to a 1D array. Note that this is not the case when X
        # is a scipy sparse array.
        # TODO: remove the following two lines when scikit-learn only depends
        # on scipy versions that no longer support scipy.sparse matrices.
        self.mean_ = xp.reshape(xp.asarray(self.mean_), (-1,))

        if self._fit_svd_solver == "full":
            X_centered = xp.asarray(X, copy=True) if self.copy else X
            X_centered -= self.mean_
            x_is_centered = not self.copy

            if not is_array_api_compliant:
                # Use scipy.linalg with NumPy/SciPy inputs for the sake of not
                # introducing unanticipated behavior changes. In the long run we
                # could instead decide to always use xp.linalg.svd for all inputs,
                # but that would make this code rely on numpy's SVD instead of
                # scipy's. It's not 100% clear whether they use the same LAPACK
                # solver by default though (assuming both are built against the
                # same BLAS).
                U, S, Vt = linalg.svd(X_centered, full_matrices=False)
            else:
                U, S, Vt = xp.linalg.svd(X_centered, full_matrices=False)
            explained_variance_ = (S**2) / (n_samples - 1)

        else:
            assert self._fit_svd_solver == "covariance_eigh"
            # In the following, we center the covariance matrix C afterwards
            # (without centering the data X first) to avoid an unnecessary copy
            # of X. Note that the mean_ attribute is still needed to center
            # test data in the transform method.
            #
            # Note: at the time of writing, `xp.cov` does not exist in the
            # Array API standard:
            # https://github.com/data-apis/array-api/issues/43
            #
            # Besides, using `numpy.cov`, as of numpy 1.26.0, would not be
            # memory efficient for our use case when `n_samples >> n_features`:
            # `numpy.cov` centers a copy of the data before computing the
            # matrix product instead of subtracting a small `(n_features,
            # n_features)` square matrix from the gram matrix X.T @ X, as we do
            # below.
            x_is_centered = False
            C = X.T @ X
            C -= (
                n_samples
                * xp.reshape(self.mean_, (-1, 1))
                * xp.reshape(self.mean_, (1, -1))
            )
            C /= n_samples - 1
            eigenvals, eigenvecs = xp.linalg.eigh(C)

            # When X is a scipy sparse matrix, the following two datastructures
            # are returned as instances of the soft-deprecated numpy.matrix
            # class. Note that this problem does not occur when X is a scipy
            # sparse array (or another other kind of supported array).
            # TODO: remove the following two lines when scikit-learn only
            # depends on scipy versions that no longer support scipy.sparse
            # matrices.
            eigenvals = xp.reshape(xp.asarray(eigenvals), (-1,))
            eigenvecs = xp.asarray(eigenvecs)

            eigenvals = xp.flip(eigenvals, axis=0)
            eigenvecs = xp.flip(eigenvecs, axis=1)

            # The covariance matrix C is positive semi-definite by
            # construction. However, the eigenvalues returned by xp.linalg.eigh
            # can be slightly negative due to numerical errors. This would be
            # an issue for the subsequent sqrt, hence the manual clipping.
            eigenvals[eigenvals < 0.0] = 0.0
            explained_variance_ = eigenvals

            # Re-construct SVD of centered X indirectly and make it consistent
            # with the other solvers.
            S = xp.sqrt(eigenvals * (n_samples - 1))
            Vt = eigenvecs.T
            U = None

        # flip eigenvectors' sign to enforce deterministic output
        U, Vt = svd_flip(U, Vt, u_based_decision=False)

        components_ = Vt

        # Get variance explained by singular values
        total_var = xp.sum(explained_variance_)
        explained_variance_ratio_ = explained_variance_ / total_var
        singular_values_ = xp.asarray(S, copy=True)  # Store the singular values.

        # Postprocess the number of components required
        if n_components == "mle":
            n_components = _infer_dimension(explained_variance_, n_samples)
        elif 0 < n_components < 1.0:
            # number of components for which the cumulated explained
            # variance percentage is superior to the desired threshold
            # side='right' ensures that number of features selected
            # their variance is always greater than n_components float
            # passed. More discussion in issue: #15669
            if is_array_api_compliant:
                # Convert to numpy as xp.cumsum and xp.searchsorted are not
                # part of the Array API standard yet:
                #
                # https://github.com/data-apis/array-api/issues/597
                # https://github.com/data-apis/array-api/issues/688
                #
                # Furthermore, it's not always safe to call them for namespaces
                # that already implement them: for instance as
                # cupy.searchsorted does not accept a float as second argument.
                explained_variance_ratio_np = _convert_to_numpy(
                    explained_variance_ratio_, xp=xp
                )
            else:
                explained_variance_ratio_np = explained_variance_ratio_
            ratio_cumsum = stable_cumsum(explained_variance_ratio_np)
            n_components = np.searchsorted(ratio_cumsum, n_components, side="right") + 1

        # Compute noise covariance using Probabilistic PCA model
        # The sigma2 maximum likelihood (cf. eq. 12.46)
        if n_components < min(n_features, n_samples):
            self.noise_variance_ = xp.mean(explained_variance_[n_components:])
        else:
            self.noise_variance_ = 0.0

        self.n_samples_ = n_samples
        self.n_components_ = n_components
        # Assign a copy of the result of the truncation of the components in
        # order to:
        # - release the memory used by the discarded components,
        # - ensure that the kept components are allocated contiguously in
        #   memory to make the transform method faster by leveraging cache
        #   locality.
        self.components_ = xp.asarray(components_[:n_components, :], copy=True)

        # We do the same for the other arrays for the sake of consistency.
        self.explained_variance_ = xp.asarray(
            explained_variance_[:n_components], copy=True
        )
        self.explained_variance_ratio_ = xp.asarray(
            explained_variance_ratio_[:n_components], copy=True
        )
        self.singular_values_ = xp.asarray(singular_values_[:n_components], copy=True)

        return U, S, Vt, X, x_is_centered, xp

    def _fit_truncated(self, X, n_components, xp):
        """Fit the model by computing truncated SVD (by ARPACK or randomized)
        on X.
        """
        n_samples, n_features = X.shape

        svd_solver = self._fit_svd_solver
        if isinstance(n_components, str):
            raise ValueError(
                "n_components=%r cannot be a string with svd_solver='%s'"
                % (n_components, svd_solver)
            )
        elif not 1 <= n_components <= min(n_samples, n_features):
            raise ValueError(
                "n_components=%r must be between 1 and "
                "min(n_samples, n_features)=%r with "
                "svd_solver='%s'"
                % (n_components, min(n_samples, n_features), svd_solver)
            )
        elif svd_solver == "arpack" and n_components == min(n_samples, n_features):
            raise ValueError(
                "n_components=%r must be strictly less than "
                "min(n_samples, n_features)=%r with "
                "svd_solver='%s'"
                % (n_components, min(n_samples, n_features), svd_solver)
            )

        random_state = check_random_state(self.random_state)

        # Center data
        total_var = None
        if issparse(X):
            self.mean_, var = mean_variance_axis(X, axis=0)
            total_var = var.sum() * n_samples / (n_samples - 1)  # ddof=1
            X_centered = _implicit_column_offset(X, self.mean_)
            x_is_centered = False
        else:
            self.mean_ = xp.mean(X, axis=0)
            X_centered = xp.asarray(X, copy=True) if self.copy else X
            X_centered -= self.mean_
            x_is_centered = not self.copy

        if svd_solver == "arpack":
            v0 = _init_arpack_v0(min(X.shape), random_state)
            U, S, Vt = svds(X_centered, k=n_components, tol=self.tol, v0=v0)
            # svds doesn't abide by scipy.linalg.svd/randomized_svd
            # conventions, so reverse its outputs.
            S = S[::-1]
            # flip eigenvectors' sign to enforce deterministic output
            U, Vt = svd_flip(U[:, ::-1], Vt[::-1], u_based_decision=False)

        elif svd_solver == "randomized":
            # sign flipping is done inside
            U, S, Vt = randomized_svd(
                X_centered,
                n_components=n_components,
                n_oversamples=self.n_oversamples,
                n_iter=self.iterated_power,
                power_iteration_normalizer=self.power_iteration_normalizer,
                flip_sign=False,
                random_state=random_state,
            )
            U, Vt = svd_flip(U, Vt, u_based_decision=False)

        self.n_samples_ = n_samples
        self.components_ = Vt
        self.n_components_ = n_components

        # Get variance explained by singular values
        self.explained_variance_ = (S**2) / (n_samples - 1)

        # Workaround in-place variance calculation since at the time numpy
        # did not have a way to calculate variance in-place.
        #
        # TODO: update this code to either:
        # * Use the array-api variance calculation, unless memory usage suffers
        # * Update sklearn.utils.extmath._incremental_mean_and_var to support array-api
        # See: https://github.com/scikit-learn/scikit-learn/pull/18689#discussion_r1335540991
        if total_var is None:
            N = X.shape[0] - 1
            X_centered **= 2
            total_var = xp.sum(X_centered) / N

        self.explained_variance_ratio_ = self.explained_variance_ / total_var
        self.singular_values_ = xp.asarray(S, copy=True)  # Store the singular values.

        if self.n_components_ < min(n_features, n_samples):
            self.noise_variance_ = total_var - xp.sum(self.explained_variance_)
            self.noise_variance_ /= min(n_features, n_samples) - n_components
        else:
            self.noise_variance_ = 0.0

        return U, S, Vt, X, x_is_centered, xp

    def score_samples(self, X):
        """Return the log-likelihood of each sample.

        See. "Pattern Recognition and Machine Learning"
        by C. Bishop, 12.2.1 p. 574
        or http://www.miketipping.com/papers/met-mppca.pdf

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data.

        Returns
        -------
        ll : ndarray of shape (n_samples,)
            Log-likelihood of each sample under the current model.
        """
        check_is_fitted(self)
        xp, _ = get_namespace(X)
        X = validate_data(self, X, dtype=[xp.float64, xp.float32], reset=False)
        Xr = X - self.mean_
        n_features = X.shape[1]
        precision = self.get_precision()
        log_like = -0.5 * xp.sum(Xr * (Xr @ precision), axis=1)
        log_like -= 0.5 * (n_features * log(2.0 * np.pi) - fast_logdet(precision))
        return log_like

    def score(self, X, y=None):
        """Return the average log-likelihood of all samples.

        See. "Pattern Recognition and Machine Learning"
        by C. Bishop, 12.2.1 p. 574
        or http://www.miketipping.com/papers/met-mppca.pdf

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data.

        y : Ignored
            Ignored.

        Returns
        -------
        ll : float
            Average log-likelihood of the samples under the current model.
        """
        xp, _ = get_namespace(X)
        return float(xp.mean(self.score_samples(X)))

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.transformer_tags.preserves_dtype = ["float64", "float32"]
        tags.array_api_support = True
        tags.input_tags.sparse = self.svd_solver in (
            "auto",
            "arpack",
            "covariance_eigh",
        )
        return tags
