import warnings
from numbers import Real, Integral
import numpy as np
from .base import BaseEstimator
from .utils import check_array
from .utils._param_validation import Interval


class MarginalSumsRegression(BaseEstimator):
    """
    Marginal Sums Regression for aggregated and non aggreagted data.

    This is a multiplicative model based on a product of factors and a base value.
    f1 * f2 * ... * fn * b

    The base value is the mean of the target variable.
    There is a factor initialzed with 1 for each feature. The features are expected to
    be onehot encoded. The first column of X can be a weight vector (see add_weights).

    The features get updated sequentially. First the marginal sum for each feature is
    calculated by X.T dot y. These sums are fixed for the rest of the algorithm. For
    each feature X gets masked, such that all rows with the current feature = 1 are
    selected and all columns except for the current feature. This masked X gets
    multiplied elementwise with the corresponding weights and current factors of the
    other features.
    SUM(weights * PROD(factors) * y_mean)
    The original marginal sum for the current feature gets divided by this estimated
    marginal sum. The result is the updated factor
    marginal sum / estimated marginal sum

    This update is done iterativley until either a number of maximum iterations is
    reached or the algorithm converges.

    Parameters
    ----------

    add_weights : bool, default=False
        If True, a numpy array of shape (X.shape[0],) gets initialized with ones.
        If False, the first column of X is considered to contain weights.

    nax_iter : int, default=100
        Number of maximum iterations, in case the algorithm does not converge. One
        iteration consists of at least one factor update of each feature.

    min_factor_change : float, default=0.001
        Criteria for early stopping. Minimal change of at least one factor in the last
        iteration.
    """

    _parameter_constraints: dict = {
        "add_weights": ["boolean"],
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "min_factor_change": [Interval(Real, 0, None, closed="left")],
    }

    def __init__(self, add_weights=False, max_iter=100, min_factor_change=0.001):
        self.add_weights = add_weights
        self.max_iter = max_iter
        self.min_factor_change = min_factor_change

    def _fit(self, X, y):
        """
        Fit the estimator by iterativly optimizing the factors for each feature.

        Parameters
        ----------

        X : array of shape (n,m)
            Input array with n observations and m features. All features need to be
            onehot encoded.

        y : array of shape (n,1)
            Target variable.
        """

        for i in range(self.max_iter):
            for feature in range(X.shape[1]):
                # Create a mask to select all rows with the current feature = 1 and
                # all columns except for the current feature.
                col_mask = [True if i != feature else False for i in range(X.shape[1])]
                row_mask = X[:, feature] > 0

                # Mask X and multiply elementwise with the current factors
                X_factor = np.multiply(self.factors[col_mask], X[row_mask][:, col_mask])

                # Calculate the marginal sum with the current factors
                # SUM(weights * PROD(factors) * y_mean)
                calc_marginal_sum = (
                    self.weights[row_mask]
                    * np.prod(X_factor, axis=1, where=X_factor > 0)
                    * self.y_mean
                ).sum()

                # Update the factor
                updated_factor = self.marginal_sums[feature] / calc_marginal_sum
                self.factors_change[feature] = np.absolute(
                    self.factors[feature] - updated_factor
                )
                self.factors[feature] = updated_factor

            # Check early stopping criteria after each iteration
            if np.max(self.factors_change) < self.min_factor_change:
                print(f"Converged after {i+1} iterations.")
                break

            if i == self.max_iter - 1:
                warnings.warn(
                    f"Did not converge after {self.max_iter} iterations.", UserWarning
                )

    def fit(self, X, y):
        """
        Wrapper for the fit_ method to calculated weights, marginal sums, the mean
        target and initialize factors.

        Parameters
        ----------

        X : array of shape (n,m)
            Input array with n observations and either m features (no weight vector) or
            m - 1 features if the first row is a weight vector. All features, except
            for the weight vector, need to be onehot encoded. The weights must not
            contain zeros.

        y : array of shape (n,1)
            Target variable.
        """

        self._validate_params()

        # ensure ndarray
        if hasattr(X, "toarray"):
            X = X.toarray()

        X = check_array(X)
        y = check_array(y)

        # init weight vector
        if self.add_weights:
            self.weights = np.ones(X.shape[0])
        else:
            self.weights = X[:, 0]
            if 0 in self.weights:
                raise ValueError("0 detected in first column. Expected weights > 0.")
            X = X[:, 1:]

        # check if array is onehot encoded
        if not ((X == 0) | (X == 1)).all():
            raise ValueError(
                "Value different from 1 or 0 detected. Only onehot encoded values"
                " expected."
            )

        # init factors
        self.factors = np.ones(X.shape[1])
        self.factors_change = np.zeros(X.shape[1])

        # calculate marginal sums of original data
        self.marginal_sums = np.dot(X.T, y)

        # calculate mean y
        self.y_mean = np.sum(y) / np.sum(self.weights)

        self._fit(X, y)

    def predict(self, X):
        """
        Predict based on the fitted model. This method expects no weight vector, only
        the onehot encoded features.

        Parameters
        -----------

        X : array of shape (n,m)
            Input array with n observations and m features. All features need to be
            onehot encoded.
        """
        if hasattr(X, "toarray"):
            X = X.toarray()
        X_factor = np.multiply(self.factors, X)
        return np.prod(X_factor, axis=1, where=X_factor > 0) * self.y_mean

    def fit_predict(self, X, y):
        """
        Fit & predict. See the fit and predcit methods for details.

        Parameters
        ----------

        X : array of shape (n,m)
            Input array with n observations and either m features (no weight vector) or
            m - 1 features if the first row is a weight vector. All features, except
            for the weight vector, need to be onehot encoded.

        y : array of shape (n,1)
            Target variable.
        """

        self.fit(X, y)
        # remove weight vector for prediction
        if not self.add_weights:
            X = X[:, 1:]
        return self.predict(X)
