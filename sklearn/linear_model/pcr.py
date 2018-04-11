import copy

import numpy as np
from ..linear_model import LinearRegression
from ..decomposition import PCA
from ..base import BaseEstimator
from ..base import RegressorMixin
from ..utils.validation import check_is_fitted
from ..utils.validation import check_X_y


class PCR(BaseEstimator, RegressorMixin):
    """
    Principal Components Regression algorithm as presented in
    Introduction to Statistical Learning by Tibshirani et al (2013).
    The underlying model is a linear regression model operating on
    reduced number of components general by PCA.

    Parameters
    ----------
    n_components : int, float, None or string
        Number of components to keep.
        if n_components is not set all components are kept::
            n_components == min(n_samples, n_features)

    copy : bool (default True)
        If False, data passed to fit are overwritten and running
        fit(X).transform(X) will not yield the expected results,
        use fit_transform(X) instead.

    whiten : bool, optional (default False)
        When True (False by default) the `components_` vectors are multiplied
        by the square root of n_samples and then divided by the singular values
        to ensure uncorrelated outputs with unit component-wise variances.
        Whitening will remove some information from the transformed signal
        (the relative variance scales of the components) but can sometime
        improve the predictive accuracy of the downstream estimators by
        making their data respect some hard-wired assumptions.

    svd_solver : string {'auto', 'full', 'arpack', 'randomized'}
        auto :
            the solver is selected by a default policy based on `X.shape` and
            `n_components`: if the input data is larger than 500x500 and the
            number of components to extract is lower than 80% of the smallest
            dimension of the data, then the more efficient 'randomized'
            method is enabled. Otherwise the exact full SVD is computed and
            optionally truncated afterwards.
        full :
            run exact full SVD calling the standard LAPACK solver via
            `scipy.linalg.svd` and select the components by postprocessing
        arpack :
            run SVD truncated to n_components calling ARPACK solver via
            `scipy.sparse.linalg.svds`. It requires strictly
            0 < n_components < min(X.shape)
        randomized :
            run randomized SVD by the method of Halko et al.

    tol : float >= 0, optional (default .0)
        Tolerance for singular values computed by svd_solver == 'arpack'.

    iterated_power : int >= 0, or 'auto', (default 'auto')
        Number of iterations for the power method computed by
        svd_solver == 'randomized'.

    random_state : int, RandomState instance or None, optional (default None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used when ``svd_solver`` == 'arpack' or 'randomized'.

    fit_intercept : boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    n_jobs : int, optional, default 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. This will only provide speedup for
        n_targets > 1 and sufficient large problems.

    Attributes
    ----------
    components_ : array, shape (n_components, n_features)
        Principal axes in feature space, representing the directions of
        maximum variance in the data. The components are sorted by
        ``explained_variance_``.
    explained_variance_ : array, shape (n_components,)
        The amount of variance explained by each of the selected components.
        Equal to n_components largest eigenvalues
        of the covariance matrix of X.
        .. versionadded:: 0.18
    explained_variance_ratio_ : array, shape (n_components,)
        Percentage of variance explained by each of the selected components.
        If ``n_components`` is not set then all components are stored and the
        sum of the ratios is equal to 1.0.
    singular_values_ : array, shape (n_components,)
        The singular values corresponding to each of the selected components.
        The singular values are equal to the 2-norms of the ``n_components``
        variables in the lower-dimensional space.
    mean_ : array, shape (n_features,)
        Per-feature empirical mean, estimated from the training set.
        Equal to `X.mean(axis=0)`.
    n_components_ : int
        The estimated number of components. When n_components is set
        to 'mle' or a number between 0 and 1
        (with svd_solver == 'full') this number is estimated from
        input data. Otherwise it equals the parameter n_components,
        or the lesser value of n_features and n_samples if n_components
        is None.
    noise_variance_ : float
        The estimated noise covariance following the Probabilistic PCA model
        from Tipping and Bishop 1999. See "Pattern Recognition and
        Machine Learning" by C. Bishop, 12.2.1 p. 574 or
        http://www.miketipping.com/papers/met-mppca.pdf. It is required to
        computed the estimated data covariance and score samples.
        Equal to the average of (min(n_features, n_samples) - n_components)
        smallest eigenvalues of the covariance matrix of X.
    coef_ : array, shape (n_features, ) or (n_targets, n_features)
        Estimated coefficients for the linear regression problem.
        If multiple targets are passed during the fit (y 2D), this
        is a 2D array of shape (n_targets, n_features), while if only
        one target is passed, this is a 1D array of length n_features.
    intercept_ : array
        Independent term in the linear model.
    """

    def __init__(self, fit_intercept=True, copy_X=True, normalize=False,
                 n_jobs=1, n_components=None, copy=True, whiten=False,
                 svd_solver='auto', tol=0.0, iterated_power='auto',
                 random_state=None):

        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.normalize = normalize
        self.n_jobs = n_jobs
        self.n_components = n_components
        self.copy = copy
        self.whiten = whiten
        self.svd_solver = svd_solver
        self.tol = tol
        self.iterated_power = iterated_power
        self.random_state = random_state

    def _get_dict(self, obj):
        """
         Get the attributes associated with the PCA or Linear Regression model.

        Parameters
        ----------
        obj : The PCA or regression model whose attributes needs to looked up.

        Returns
        -------

        """

        for attr, value in obj.__dict__.items():
            self.__dict__[attr] = value

    def _set_attributes(self, model):
        """
         Set the attributes associated with the PCA or Linear Regression
         model to passed estimator.

        Parameters
        ----------
        model : Specifies the PCA or linear regression model whose attributes
        to be set

        Returns
        -------

        """

        if model:
            self._get_dict(model)

    def get_regression_model(self):
        """
        The associated linear regression model for PCR.

        Parameters
        ----------

        Returns
        -------
        _lr: Returns the unfitted linear regression model of PCR.

        """

        return self._lr

    def get_transformed_data(self):
        """
        Generates the reduced principal components from the training data
        using PCA.

        Parameters
        ----------

        Returns
        -------
        transformed: Returns a transformed numpy array or sparse matrix. The
        data has been reduced to the specified number of components using the
        `n_components` parameter of the PCR model.

        Notes
        -------
        The algorithm should have first been executed on a training dataset.

        """

        transformed = self._X_reduced
        return transformed

    def _reduce(self, X):
        """
        Fit the Principal Components Analysis model.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            Training data

        Returns
        -------
        """
        if not self.random_state:
            random_state = 0
        else:
            random_state = self.random_state

        self._pca = PCA(n_components=self.n_components, copy=self.copy,
                        whiten=self.whiten, svd_solver=self.svd_solver,
                        tol=self.tol, iterated_power=self.iterated_power,
                        random_state=random_state)

        self._pca.fit(X)
        self._X_reduced = self._pca.transform(X)
        self._set_attributes(self._pca)

    def _regress(self, X, y):
        """
        Fit the Linear Regression model.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples, n_targets]
            Target values

        Returns
        -------
        """
        self._model.fit(X, np.ravel(y))
        self._set_attributes(self._model)

    def fit(self, X, y):
        """
        Fit the Principal Components Regression model.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples, n_targets]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_X_y(X, y)

        self._model = LinearRegression(fit_intercept=self.fit_intercept,
                                       copy_X=self.copy_X,
                                       normalize=self.normalize,
                                       n_jobs=self.n_jobs)

        self._lr = copy.deepcopy(self._model)

        self._reduce(X)
        self._regress(self._X_reduced, y)

        return self

    def predict(self, X):
        """
        Predict using the PCR model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.

        Returns
        -------
        predictions : array, shape = (n_samples,)
            Returns predicted values.
        """
        check_is_fitted(self, "coef_")
        X_reduced_test = self._pca.transform(X)
        predictions = self._model.predict(X_reduced_test)
        return predictions
