""" Matrix factorization with missing values
"""

import numpy as np
import scipy.sparse as sparse

from ..utils.validation import array2d, check_random_state
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

    bias_learning_rate : double, default: TODO
        The learning_rate of the sample/features bias

    bias_regularization : double, default: TODO
        The regularization coefficient of the sample/features bias

    init_L : array, [n_samples, n_components]
        An initial guess for the first factor

    init_R : array, [n_components, n_features]
        An initial guess for the second

    init_bias_samples : array, [n_samples, 1]
        An initial guess for the bias fo the samples

    init_bias_features : array, [n_features, 1]
        An initial guess for the bias fo the features

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
                 bias_learning_rate=0, bias_regularization=0,
                 init_L=None, init_R=None,
                 init_bias_samples=None, init_bias_features=None,
                 random_state=None, verbose=0):

        self.n_components = n_components
        self.n_iter = n_iter

        self.learning_rate = learning_rate
        self.regularization = regularization
        self.missing_values = missing_values
        self.bias_learning_rate = bias_learning_rate
        self.bias_regularization = bias_regularization

        self.start_L = init_L
        self.start_R= init_R
        self.start_bias_samples = init_bias_samples
        self.start_bias_features = init_bias_features

        self.random_state = random_state
        self.verbose = verbose

    def fit_transform(self, X, y=None):
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
            X_content = np.vstack((X.data, X.row, X.col)).T
        else:
            X = array2d(X, force_all_finite=False)
            mask = np.logical_not(_get_mask(X, self.missing_values))
            X_content = np.vstack([X[mask]] + list(np.where(mask))).T
        mean = X_content[:, 0].mean()
        var  = X_content[:, 0].var()

        print mean

        # Initialize the factors
        if self.start_L is None or self.start_R is None:

            # Initialize L and R
            # such that for all i, j: L[i, :] * R[:, j] = X.mean()
            initial_value = np.sqrt(mean / n_components)
            L = np.zeros((n_samples, n_components))
            R = np.zeros((n_components, n_features))
            L.fill(initial_value)
            R.fill(initial_value)

            # Add some noise
            r = np.sqrt(var)
            L = L + random_state.uniform(
                -r, r,
                n_samples * n_components).reshape(L.shape)
            R = R + random_state.uniform(
                -r, r,
                n_components * n_features).reshape(R.shape)

        else:
            L = self.start_L.copy()
            R = self.start_R.copy()

        # Initialize the bias
        if self.start_bias_samples is None or self.start_bias_features is None:
            bias_samples = np.zeros(n_samples)
            bias_features = np.zeros(n_features)
        else:
            bias_samples = self.start_bias_samples.copy()
            bias_features = self.start_bias_features.copy()

        # Factorize the matrix
        # TODO: how to handle possible overflows? Should they really happen?
        _mf.factorize_matrix(
             X_content[:, 0],
             X_content[:, 1].astype(np.int32),
             X_content[:, 2].astype(np.int32),
             L, R, bias_samples, bias_features,
             n_samples, n_features, self.n_iter, n_components,
             self.learning_rate, self.regularization,
             self.bias_learning_rate, self.bias_regularization,
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

    def fit(self, X, y=None, **params):


        self.fit_transform(X, **params)
        return self