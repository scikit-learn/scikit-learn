""" Matrix factorization with missing values
"""

import numpy as np
import scipy.sparse as sparse

from ..utils.validation import check_array, check_random_state
from ..base import BaseEstimator, TransformerMixin

import _mf


_rmse = _mf._rmse


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    # TODO: if at the end this function still does the same job that the one in
    # sklearn.preprocessing.imputation, call it directly
    if value_to_mask == "NaN" or np.isnan(value_to_mask):
        return np.isnan(X)
    else:
        return X == value_to_mask

class MatrixFactorization(BaseEstimator, TransformerMixin):
    """Matrix Factorization

    Parameters
    ----------
    n_components : int or None
        Estimated rank of the matrix that will be factorized, if it's not
        set, the matrix will be considered to have full rank

    n_iter : int, default: 100
        Number of iterations to compute.

    missing_values : integer or string, optional (default="NaN")
        The placeholder for the missing values. If the matrix to factorize
        is sparse, this parameter will be ignored (the placeholder will be 0)

    learning_rate : double, default: TODO
        The learning_rate

    regularization : double, default: TODO
        The regularization coefficient

    bias_regularization : double, default: TODO
        The regularization coefficient of the sample/features bias

    init_L : array, [n_samples, n_components]
        An initial guess for the first factor

    init_R : array, [n_components, n_features]
        An initial guess for the second


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

    TODO
    """

    def __init__(self, n_components=None, n_iter=100, missing_values="NaN",
                 learning_rate=1e-3, regularization=1e-4,
                 bias_regularization=0,
                 init_L=None, init_R=None,
                 random_state=None, verbose=0):

        self.n_components = n_components
        self.n_iter = n_iter

        self.learning_rate = learning_rate
        self.regularization = regularization
        self.missing_values = missing_values
        self.bias_regularization = bias_regularization

        self.init_L = init_L
        self.init_R= init_R

        self.random_state = random_state
        self.verbose = verbose

    def fit_transform(self, X, y=None):
        # TODO: how come this is more efficient?
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
        """
        n_samples, n_features = X.shape
        random_state = check_random_state(self.random_state)

        n_components = self.n_components
        if n_components is None:
            n_components = min(n_samples, n_features)

        # Extract the non-missing values and their indices
        if sparse.issparse(X):
            X = X.tocoo()
            X_data = X.data
            X_rows = X.row
            X_cols = X.col
        else:
            X = check_array(X, force_all_finite=False)
            mask = np.logical_not(_get_mask(X, self.missing_values))
            X_data = X[mask]
            X_rows, X_cols = list(np.where(mask))

        mean = X_data.mean()
        var  = X_data.var()

        # Initialize the factors
        if self.init_L is None or self.init_R is None:

            # Initialize L and R
            # such that for all i, j: L[i, :] * R[:, j] = X.mean()
            initial_value = np.sqrt(mean / n_components)
            L = np.zeros((n_samples, n_components))
            R = np.zeros((n_components, n_features))
            L.fill(initial_value)
            R.fill(initial_value)

            # Add some noise with var = var(X)
            # TODO: shouldn't we add gaussian noise here?
            r = np.sqrt(var) * np.sqrt(12) / 2.0
            L += random_state.uniform(-r, r, L.shape)
            R += random_state.uniform(-r, r, R.shape)

        else:
            L = self.init_L.copy()
            R = self.init_R.copy()

        # Initialize the bias
        bias_samples = np.zeros(n_samples)
        bias_features = np.zeros(n_features)

        # Factorize the matrix
        # TODO: how to handle possible overflows? Should they really happen?
        _mf.factorize_matrix(
             X_data,
             X_rows.astype(np.int32),
             X_cols.astype(np.int32),
             L, R, bias_samples, bias_features,
             n_samples, n_features, n_components, self.n_iter,
             self.regularization, self.bias_regularization,
             self.learning_rate,
             random_state, self.verbose,
        )

        self.components_ = np.vstack([
            bias_features.reshape(1, n_features),
            np.ones((1, n_features)),
            R])

        return np.hstack([
            np.ones((n_samples, 1)),
            bias_samples.reshape(n_samples, 1),
            L])

    def fit(self, X, y=None):

        self.fit_transform(X)
        return self