import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted

class LinearRegressor(BaseEstimator, RegressorMixin):
    """
    Linear regression model with optional L1/L2 regularization.

    Parameters
    ----------
    epochs : int, default=1000
        Number of iterations for gradient descent.

    learning_rate : float, default=0.01
        Step size for gradient updates.

    lamd : float, default=0.0
        Regularization strength.

    reg : {'L1', 'L2', None}, default=None
        Type of regularization to apply.

    gd : {'batch', 'stochastic', 'mini-batch'}, default='batch'
        Gradient descent method.

    batch_size : int, default=1
        Size of batches for mini-batch gradient descent.

    seed : int, default=None
        Random seed for shuffling.

    Attributes
    ----------
    coef_ : ndarray of shape (n_features,)
        Estimated coefficients for the linear model.

    intercept_ : float
        Intercept term.

    loss_history_ : list
        List of loss values at each epoch.
    """

    def __init__(self, epochs=1000, learning_rate=0.01, lamd=0.0, reg=None,
                 gd='batch', batch_size=1, seed=None):
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.lamd = lamd
        self.reg = reg
        self.gd = gd
        self.batch_size = batch_size
        self.seed = seed

    def _apply_regularization(self, grad, m):
        if self.reg is None:
            return grad
        elif self.reg == 'L1':
            return grad + self.lamd * np.sign(m)
        elif self.reg == 'L2':
            return grad + 2 * self.lamd * m
        else:
            raise ValueError(f"Unsupported regularization type: {self.reg}")

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        y = y.reshape(-1, 1)
        m = np.zeros(X.shape[1])
        b = 0.0
        self.loss_history_ = []

        rng = np.random.RandomState(self.seed)

        for _ in range(self.epochs):
            if self.gd == 'batch':
                y_pred = X @ m + b
                error = y - y_pred
                grad_m = (-2 / len(y)) * (X.T @ error).flatten()
                grad_m = self._apply_regularization(grad_m, m)
                m -= self.learning_rate * grad_m
                b -= self.learning_rate * (-2 / len(y)) * np.sum(error)

            elif self.gd == 'stochastic':
                for j in range(len(y)):
                    y_pred = np.dot(X[j], m) + b
                    error = y[j] - y_pred
                    grad_m = -2 * X[j] * error
                    grad_m = self._apply_regularization(grad_m, m)
                    m -= self.learning_rate * grad_m
                    b -= self.learning_rate * (-2) * error

            elif self.gd == 'mini-batch':
                indices = rng.permutation(len(y))
                X_shuffled = X[indices]
                y_shuffled = y[indices]
                for i in range(0, len(y), self.batch_size):
                    xb = X_shuffled[i:i+self.batch_size]
                    yb = y_shuffled[i:i+self.batch_size]
                    y_pred = xb @ m + b
                    error = yb - y_pred
                    grad_m = (-2 / len(yb)) * (xb.T @ error).flatten()
                    grad_m = self._apply_regularization(grad_m, m)
                    m -= self.learning_rate * grad_m
                    b -= self.learning_rate * (-2 / len(yb)) * np.sum(error)

            else:
                raise ValueError(f"Unsupported gradient descent type: {self.gd}")

            loss = np.mean((y - (X @ m + b).reshape(-1,1)) ** 2)
            self.loss_history_.append(loss)

        self.coef_ = m
        self.intercept_ = b
        return self

    def predict(self, X):
        check_is_fitted(self, ['coef_', 'intercept_'])
        X = check_array(X)
        return X @ self.coef_ + self.intercept_

   

        
    