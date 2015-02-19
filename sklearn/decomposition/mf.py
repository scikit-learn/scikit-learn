""" Matrix factorization with missing values
"""

import numpy as np
import scipy.sparse as sparse

from ..utils import as_float_array, get_mask
from ..utils.validation import check_array, check_random_state, check_is_fitted
from ..base import BaseEstimator, TransformerMixin

import _mf

_rmse = _mf._rmse


class MatrixFactorization(BaseEstimator, TransformerMixin):
    """Matrix Factorization

    Parameters
    ----------
    n_components : int or None
        Estimated rank of the matrix that will be factorized, if it's not
        set, the matrix will be considered to have full rank

    algorithm : {'sgd', 'sgd_adagrad', 'als', 'als1'}, default: 'sgd'
        Algorithm to perform optimization.
        sgd and sgd_adagrad are variants of stochastic gradient descent,
        the later one uses AdaGrad learning rate adjustment.
        ALS and ALS1 are variants of Alternating Least Squares, the later
        one avoids inverting matrices, thus should be faster for large n_components

    n_iter : int, default: 100
        Number of iterations to compute.

    missing_values : integer or string, optional (default="NaN")
        The placeholder for the missing values. If the matrix to factorize
        is sparse, this parameter will be ignored (the placeholder will be 0)

    learning_rate : double, default: 1e-3
        The learning_rate for SGD algorithm

    regularization : double, default: 1e-4
        The regularization coefficient. Note: bias (if turned on)
        is not regularized.

    fit_intercept : boolean, default: True
        Specifies if a bias (a.k.a. intercept) should be added to the factorization.

    init_L : array, [n_samples, n_components]
        An initial guess for the first factor

    init_R : array, [n_components, n_features]
        An initial guess for the second

    init_sample_biases : array, [n_samples, n_components]
        An initial guess for biases for the first (sample) factor

    init_feature_biases : array, [n_components, n_features]
        An initial guess for the biases for the second (feature) factor

    verbose : int, default: 0
        Level of verbosity of output. 0 means no output, 1 - output
        loss each k iterations where k is chosen so that no more than
        100 messages are outputted.

    Attributes
    ----------
    `components_` : array, [n_components + 2, n_features]
        The second factor of the data. The bias of the features is located in
            the second row. The first row will always be 1

    Examples
    --------

    TODO

    References
    ----------
    This implements

    Istvan Pilaszy, David Zibriczky: ALS and ALS1: Fast als-based matrix factorization
    for explicit and implicit feedback datasets
    http://dl.acm.org/citation.cfm?id=1864726
    """

    def __init__(self, n_components=None, n_iter=100, missing_values="NaN",
                 learning_rate=1e-3, regularization=1e-4,
                 fit_intercept=True,
                 init_L=None, init_R=None,
                 init_sample_biases=None,
                 init_feature_biases=None,
                 random_state=None,
                 algorithm='sgd',
                 verbose=0):

        self.n_components = n_components
        self.n_iter = n_iter

        self.learning_rate = learning_rate
        self.regularization = regularization
        self.missing_values = missing_values
        self.fit_intercept = fit_intercept

        self.init_L = init_L
        self.init_R= init_R

        self.init_sample_biases = init_sample_biases
        self.init_feature_biases = init_feature_biases

        self.random_state = random_state
        self.verbose = verbose
        self.algorithm = algorithm

    def _extract_data(self, X):
        """
        Converts data matrix into 3 lists of present values.
        The first list is a list of values, and the last two correspond
        to lists of row and column indices.
        """
        if sparse.issparse(X):
            X_data = X.data
            X_rows = X.row
            X_cols = X.col

            if self.missing_values != 0:
                mask = np.logical_not(get_mask(X_data, self.missing_values))
                X_data = X_data[mask]
                indices = np.where(mask)
                X_rows = X_rows[indices]
                X_cols = X_cols[indices]
        else:
            mask = np.logical_not(get_mask(X, self.missing_values))
            X_data = X[mask]
            X_rows, X_cols = list(np.where(mask))

        # Standardizes data order
        # TODO: do we really need it?
        standardized = sorted(zip(X_data, X_rows, X_cols))
        X_data = np.array([x[0] for x in standardized])
        X_rows = np.array([x[1] for x in standardized])
        X_cols = np.array([x[2] for x in standardized])
        return X_data, X_rows, X_cols

    def fit_transform(self, X, y=None, update=True):
        """Factorize the matrix and return the first factor

        This is more efficient than calling fit followed by transform.

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be factorized

        Returns
        -------
        data: array, [n_samples, n_components + 2]
            The first factor of the data. The bias of the samples is located in
            the first column. The second column will always be 1

            If fit_intercept is False, the bias column is zero.
        """
        X = check_array(X, force_all_finite=False, accept_sparse='coo')
        n_samples, n_features = X.shape
        random_state = check_random_state(self.random_state)

        n_components = self.n_components
        if n_components is None:
            n_components = min(n_samples, n_features)

        # Extract the non-missing values and their indices
        X_data, X_rows, X_cols = self._extract_data(X)

        if len(set(X_rows)) < n_samples:
            raise ValueError("X should not have rows of missing values only")

        if len(set(X_cols)) < n_features:
            raise ValueError("X should not have columns of missing values only")

        mean = X_data.mean()
        var  = X_data.var()

        # Initialize L and R
        # Add some noise with var = var(X)
        if self.init_L is None:
            L = random_state.normal(0, np.sqrt(var), (n_samples, n_components))
        else:
            L = self.init_L.copy()

        if self.init_R is None:
            R = random_state.normal(0, np.sqrt(var), (n_components, n_features))
        else:
            R = self.init_R.copy()

        # Initialize the biases
        # TODO: should we add another (global) bias into the model?
        if self.init_sample_biases is None:
            bias_samples = np.zeros(n_samples)
        else:
            bias_samples = self.init_sample_biases

        if self.init_feature_biases is None:
            bias_features = np.zeros(n_features)
        else:
            bias_features = self.init_feature_biases

        # Correct biases' means
        if self.fit_intercept:
            correction = mean - bias_features.mean() - bias_samples.mean()
            bias_features += correction / 2
            bias_samples += correction / 2
        else:
            v = mean - L.mean() * R.mean()
            correction = np.sqrt(abs(v))
            sign = np.sign(v)
            R += correction * sign
            L += correction

        # Factorize the matrix
        if self.algorithm in ["als", "als1"]:
            als1 = self.algorithm == "als1"
            _mf.factorize_matrix_als(
                X_data,
                X_rows.astype(np.int32),
                X_cols.astype(np.int32),
                L, R, bias_samples, bias_features,
                n_samples, n_features, n_components,
                self.n_iter, int(self.fit_intercept),
                self.regularization, self.verbose, als1
            )
        elif self.algorithm in ["sgd", "sgd_adagrad"]:
            adagrad = self.algorithm == "sgd_adagrad"
            _mf.factorize_matrix_sgd(
                X_data,
                X_rows.astype(np.int32),
                X_cols.astype(np.int32),
                L, R, bias_samples, bias_features,
                n_samples, n_features, n_components,
                self.n_iter, int(self.fit_intercept),
                self.regularization, self.learning_rate,
                random_state, self.verbose, adagrad
            )
        else:
            raise ValueError("Unknown algorithm: %s" % self.algorithm)

        if update:
            self.feature_vectors_ = R
            self.components_ = np.vstack([
                bias_features.reshape(1, n_features),
                np.ones((1, n_features)),
                R])

        return np.hstack([
            np.ones((n_samples, 1)),
            bias_samples.reshape(n_samples, 1),
            L])

    def fit(self, X, y=None):
        """Factorize the matrix

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be factorized
        """
        self.fit_transform(X)
        return self

    def transform(self, X, y=None):
        """Factorize the matrix and return the first factor

        Due to details of the factorization problem, the factorization
        needs to be performed on each transform call, so it's not a cheap
        operation.

        Parameters
        ----------

        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data matrix to be factorized

        Returns
        -------
        data: array, [n_samples, n_components + 2]
            The first factor of the data. The bias of the samples is located in
            the first column. The second column will always be 1

            If fit_intercept is False, the bias column is zero.
        """
        X = as_float_array(X, copy=False, force_all_finite=False)
        check_is_fitted(self, 'feature_vectors_')

        if X.shape[1] != self.feature_vectors_.shape[1]:
            raise ValueError("X needs to have the same number of features")

        old_init_R = self.init_R

        self.init_R = self.feature_vectors_
        res = self.fit_transform(X, update=False)

        self.init_R = old_init_R

        return res