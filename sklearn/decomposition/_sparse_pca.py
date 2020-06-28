"""Matrix factorization with Sparse PCA"""
# Author: Vlad Niculae, Gael Varoquaux, Alexandre Gramfort
# License: BSD 3 clause

import numpy as np
from scipy import linalg

from ..utils import check_random_state, check_array
from ..utils.validation import check_is_fitted
from ..utils.validation import _deprecate_positional_args
from ..linear_model import ridge_regression
from ..base import BaseEstimator, TransformerMixin
from ._dict_learning import dict_learning, dict_learning_online


def _get_explained_variance(X, components, ridge_alpha):
    """Get the explained variance.

    Get the explained variance from the principal components of the
    data. This follows the method outlined in [1] section 3.4 (Adjusted Total
    Variance). For an alternate approach (not implemented here), see [2].

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        The feature vector. n_samples and n_features are the number of
        samples and features, respectively.

    components : array of shape (n_components, n_features)
        The principal components. [1]

    ridge_alpha : float
        Amount of ridge shrinkage to apply in order to improve
        conditioning when calling the transform method.

    Notes
    -----
    Orthogonality is enforced in this case. Other variants exist that don't
    enforce this [2].

    Examples
    --------
    R implementation using elasticnet 1.1.1

    > library(elasticnet)
    > sparse_pca <- spca(x = iris[,c(1, 2, 3, 4)], K = 1, para = c(1))
    > sparse_pca$pev
    0.9090568
    > sparse_pca$loadings
                         PC1
    Sepal.Length  0.30698724
    Sepal.Width  -0.02506061
    Petal.Length  0.91173203
    Petal.Width   0.27180050

    >>> import numpy as np
    >>> from numpy.testing import assert_almost_equal
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.decomposition._sparse_pca import _get_explained_variance
    >>> data = load_iris()["data"]
    >>> data = data - data.mean(axis=0)
    >>> _, pev = _get_explained_variance(data,\
    np.array([[0.30698724, -0.02506061,  0.91173203,  0.27180050]]), 1e-8)
    >>> assert_almost_equal(pev, 0.9090568)

    References
    ----------
    .. [1] Hui Zou, Trevor Hastie, Robert Tibshirani "Sparse Principal
        Component Analysis", Journal of Computational and Graphical
        Statistics, Volume 15, Number 2, Pages 265–286.
        DOI: 10.1198/106186006X113430
    .. [2] Rodolphe Jenatton, Guillaume Obozinski, Francis Bach, "Structured
        sparse principal component analysis", Proceedings of the Thirteenth
        International Conference on Artificial Intelligence and Statistics
        (AISTATS), PMLR 9:366-373, 2010.
    """
    # Transform input
    U = ridge_regression(components.T, X.T, ridge_alpha, solver="cholesky")

    # QR decomposition of sparse PCs
    r = linalg.qr(U, overwrite_a=True, mode='r')[0]
    variance = np.square(np.diag(r)) / (X.shape[0]-1)

    # Variance in the original dataset
    cov_mat = np.cov(X.T)
    total_variance_in_x = np.trace(cov_mat) if cov_mat.size >= 2 \
        else cov_mat
    return variance,  variance / total_variance_in_x


class SparsePCA(TransformerMixin, BaseEstimator):
    """Sparse Principal Components Analysis (SparsePCA)

    Finds the set of sparse components that can optimally reconstruct
    the data.  The amount of sparseness is controllable by the coefficient
    of the L1 penalty, given by the parameter alpha.

    Read more in the :ref:`User Guide <SparsePCA>`.

    Parameters
    ----------
    n_components : int,
        Number of sparse atoms to extract.

    alpha : float,
        Sparsity controlling parameter. Higher values lead to sparser
        components.

    ridge_alpha : float,
        Amount of ridge shrinkage to apply in order to improve
        conditioning when calling the transform method.

    max_iter : int,
        Maximum number of iterations to perform.

    tol : float,
        Tolerance for the stopping condition.

    method : {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    n_jobs : int or None, optional (default=None)
        Number of parallel jobs to run.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    U_init : array of shape (n_samples, n_components),
        Initial values for the loadings for warm restart scenarios.

    V_init : array of shape (n_components, n_features),
        Initial values for the components for warm restart scenarios.

    verbose : int
        Controls the verbosity; the higher, the more messages. Defaults to 0.

    random_state : int, RandomState instance, default=None
        Used during dictionary learning. Pass an int for reproducible results
        across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    components_ : array of shape (n_components, n_features)
        Sparse components extracted from the data.

    error_ : array
        Vector of errors at each iteration.

    n_components_ : int
        Estimated number of components.

        .. versionadded:: 0.23

    n_iter_ : int
        Number of iterations run.

    mean_ : array of shape (n_features,)
        Per-feature empirical mean, estimated from the training set.
        Equal to ``X.mean(axis=0)``.

    explained_variance_ : array of shape (n_components,)
        The amount of variance explained by each of the selected components.

    explained_variance_ratio_ : array of shape (n_components,)
        Percentage of variance explained by each of the selected components.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.decomposition import SparsePCA
    >>> X, _ = make_friedman1(n_samples=200, n_features=30, random_state=0)
    >>> transformer = SparsePCA(n_components=5, random_state=0)
    >>> transformer.fit(X)
    SparsePCA(...)
    >>> X_transformed = transformer.transform(X)
    >>> X_transformed.shape
    (200, 5)
    >>> # most values in the components_ are zero (sparsity)
    >>> np.mean(transformer.components_ == 0)
    0.9666...

    Notes
    -----
    In computing the explained variance, an orthogonality-enforcing method was
    chosen [1]. However, there exist other approaches, such as in [2] (which
    are not implemented here).

    See also
    --------
    PCA
    MiniBatchSparsePCA
    DictionaryLearning

    References
    ----------
    .. [1] Hui Zou, Trevor Hastie, Robert Tibshirani "Sparse Principal
        Component Analysis", Journal of Computational and Graphical
        Statistics, Volume 15, Number 2, Pages 265–286.
        DOI: 10.1198/106186006X113430
    .. [2] Rodolphe Jenatton, Guillaume Obozinski, Francis Bach, "Structured
        sparse principal component analysis", Proceedings of the Thirteenth
        International Conference on Artificial Intelligence and Statistics
        (AISTATS), PMLR 9:366-373, 2010.
    """
    @_deprecate_positional_args
    def __init__(self, n_components=None, *, alpha=1, ridge_alpha=0.01,
                 max_iter=1000, tol=1e-8, method='lars', n_jobs=None,
                 U_init=None, V_init=None, verbose=False, random_state=None):
        self.n_components = n_components
        self.alpha = alpha
        self.ridge_alpha = ridge_alpha
        self.max_iter = max_iter
        self.tol = tol
        self.method = method
        self.n_jobs = n_jobs
        self.U_init = U_init
        self.V_init = V_init
        self.verbose = verbose
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        y : Ignored

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        random_state = check_random_state(self.random_state)
        X = self._validate_data(X)

        self.mean_ = X.mean(axis=0)
        X = X - self.mean_

        if self.n_components is None:
            n_components = X.shape[1]
        else:
            n_components = self.n_components
        code_init = self.V_init.T if self.V_init is not None else None
        dict_init = self.U_init.T if self.U_init is not None else None
        Vt, _, E, self.n_iter_ = dict_learning(X.T, n_components,
                                               alpha=self.alpha,
                                               tol=self.tol,
                                               max_iter=self.max_iter,
                                               method=self.method,
                                               n_jobs=self.n_jobs,
                                               verbose=self.verbose,
                                               random_state=random_state,
                                               code_init=code_init,
                                               dict_init=dict_init,
                                               return_n_iter=True)
        self.components_ = Vt.T
        components_norm = np.linalg.norm(
            self.components_, axis=1)[:, np.newaxis]
        components_norm[components_norm == 0] = 1
        self.components_ /= components_norm
        self.n_components_ = len(self.components_)

        self.explained_variance_, self.explained_variance_ratio_ = \
            _get_explained_variance(X, self.components_, self.ridge_alpha)
        self.error_ = E
        return self

    def transform(self, X):
        """Least Squares projection of the data onto the sparse components.

        To avoid instability issues in case the system is under-determined,
        regularization can be applied (Ridge regression) via the
        `ridge_alpha` parameter.

        Note that Sparse PCA components orthogonality is not enforced as in PCA
        hence one cannot use a simple linear projection.

        Parameters
        ----------
        X : array of shape (n_samples, n_features)
            Test data to be transformed, must have the same number of
            features as the data used to train the model.

        Returns
        -------
        X_new array of shape (n_samples, n_components)
            Transformed data.
        """
        check_is_fitted(self)

        X = check_array(X)
        X = X - self.mean_

        U = ridge_regression(self.components_.T, X.T, self.ridge_alpha,
                             solver='cholesky')

        return U


class MiniBatchSparsePCA(SparsePCA):
    """Mini-batch Sparse Principal Components Analysis

    Finds the set of sparse components that can optimally reconstruct
    the data.  The amount of sparseness is controllable by the coefficient
    of the L1 penalty, given by the parameter alpha.

    Read more in the :ref:`User Guide <SparsePCA>`.

    Parameters
    ----------
    n_components : int,
        number of sparse atoms to extract

    alpha : int,
        Sparsity controlling parameter. Higher values lead to sparser
        components.

    ridge_alpha : float,
        Amount of ridge shrinkage to apply in order to improve
        conditioning when calling the transform method.

    n_iter : int,
        number of iterations to perform for each mini batch

    callback : callable or None, optional (default: None)
        callable that gets invoked every five iterations

    batch_size : int,
        the number of features to take in each mini batch

    verbose : int
        Controls the verbosity; the higher, the more messages. Defaults to 0.

    shuffle : boolean,
        whether to shuffle the data before splitting it in batches

    n_jobs : int or None, optional (default=None)
        Number of parallel jobs to run.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    method : {'lars', 'cd'}
        lars: uses the least angle regression method to solve the lasso problem
        (linear_model.lars_path)
        cd: uses the coordinate descent method to compute the
        Lasso solution (linear_model.Lasso). Lars will be faster if
        the estimated components are sparse.

    random_state : int, RandomState instance, default=None
        Used for random shuffling when ``shuffle`` is set to ``True``,
        during online dictionary learning. Pass an int for reproducible results
        across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    components_ : array of shape (n_components, n_features)
        Sparse components extracted from the data.

    n_components_ : int
        Estimated number of components.

        .. versionadded:: 0.23

    n_iter_ : int
        Number of iterations run.

    mean_ : array of shape (n_features,)
        Per-feature empirical mean, estimated from the training set.
        Equal to ``X.mean(axis=0)``.

    explained_variance_ : array of shape (n_components,)
        The amount of variance explained by each of the selected components.

    explained_variance_ratio_ : array of shape (n_components,)
        Percentage of variance explained by each of the selected components.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.decomposition import MiniBatchSparsePCA
    >>> X, _ = make_friedman1(n_samples=200, n_features=30, random_state=0)
    >>> transformer = MiniBatchSparsePCA(n_components=5, batch_size=50,
    ...                                  random_state=0)
    >>> transformer.fit(X)
    MiniBatchSparsePCA(...)
    >>> X_transformed = transformer.transform(X)
    >>> X_transformed.shape
    (200, 5)
    >>> # most values in the components_ are zero (sparsity)
    >>> np.mean(transformer.components_ == 0)
    0.94

    See also
    --------
    PCA
    SparsePCA
    DictionaryLearning
    """
    @_deprecate_positional_args
    def __init__(self, n_components=None, *, alpha=1, ridge_alpha=0.01,
                 n_iter=100, callback=None, batch_size=3, verbose=False,
                 shuffle=True, n_jobs=None, method='lars', random_state=None):
        super().__init__(
            n_components=n_components, alpha=alpha, verbose=verbose,
            ridge_alpha=ridge_alpha, n_jobs=n_jobs, method=method,
            random_state=random_state)
        self.n_iter = n_iter
        self.callback = callback
        self.batch_size = batch_size
        self.shuffle = shuffle

    def fit(self, X, y=None):
        """Fit the model from data in X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        y : Ignored

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        random_state = check_random_state(self.random_state)
        X = self._validate_data(X)

        self.mean_ = X.mean(axis=0)
        X = X - self.mean_

        if self.n_components is None:
            n_components = X.shape[1]
        else:
            n_components = self.n_components
        Vt, _, self.n_iter_ = dict_learning_online(
            X.T, n_components, alpha=self.alpha,
            n_iter=self.n_iter, return_code=True,
            dict_init=None, verbose=self.verbose,
            callback=self.callback,
            batch_size=self.batch_size,
            shuffle=self.shuffle,
            n_jobs=self.n_jobs, method=self.method,
            random_state=random_state,
            return_n_iter=True)
        self.components_ = Vt.T
        components_norm = np.linalg.norm(
            self.components_, axis=1)[:, np.newaxis]
        components_norm[components_norm == 0] = 1
        self.components_ /= components_norm
        self.n_components_ = len(self.components_)

        self.explained_variance_, self.explained_variance_ratio_ = \
            _get_explained_variance(X, self.components_, self.ridge_alpha)

        return self
