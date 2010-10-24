# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with sparse data."""

import numpy as np
from scipy import sparse

from ...base import ClassifierMixin
from ..base import LinearModel
from . import sgd_fast_sparse


class SGD(LinearModel, ClassifierMixin):
    """Linear Model trained by minimizing a regularized training
    error using SGD.

    This implementation works on scipy.sparse X and dense coef_.

    Parameters
    ----------
    loss : str, ('hinge'|'log'|'modifiedhuber')
        The loss function to be used.
    penalty : str, ('l2'|'l1'|'elasticnet')
        The penalty (aka regularization term) to be used.
    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001
    rho : float
        The Elastic Net mixing parameter, with 0 < rho <= 1.
    coef_ : ndarray of shape n_features
        The initial coeffients to warm-start the optimization
    intercept_ : float
        The initial intercept to warm-start the optimization
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.
    n_iter: int
        The number of passes over the training data (aka epochs).
    shuffle: bool
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        Weights asigned to the features.

    `intercept_` : float
        Constants in decision function.

    """

    def _get_loss_function(self):
        loss_functions = {"hinge" : sgd_fast_sparse.Hinge(),
                          "log" : sgd_fast_sparse.Log(),
                          "modifiedhuber" : sgd_fast_sparse.ModifiedHuber(),
                          }
        try:
            self.loss_function = loss_functions[self.loss]
        except KeyError:
            raise ValueError("The loss %s is not supported. " % self.loss)

    def _set_coef(self, coef_):
        self.coef_ = coef_
        if coef_ is None:
            self.sparse_coef_ = None
        else:
            # sparse representation of the fitted coef for the predict method
            self.sparse_coef_ = sparse.csr_matrix(coef_)

    def fit(self, X, Y, **params):
        """Fit current model with SGD

        X is expected to be a sparse matrix. For maximum efficiency, use a
        sparse matrix in CSR format (scipy.sparse.csr_matrix)
        """
        self._set_params(**params)
        X = sparse.csr_matrix(X)
        Y = np.asanyarray(Y, dtype=np.float64)
        classes = np.unique(Y)
        if len(classes) != 2:
            raise ValueError("SGD supports binary classification only.")
        self.classes = classes

        # encode original class labels as 1 (classes[0]) or -1 (classes[1]).
        Y_new = np.ones(Y.shape, dtype=np.float64)
        Y_new[Y == classes[1]] = - 1.0
        Y = Y_new

        n_samples, n_features = X.shape[0], X.shape[1]
        if self.coef_ is None:
            self.coef_ = np.zeros(n_features, dtype=np.float64, order="c")

        X_data = np.array(X.data, dtype=np.float64, order="c")
        X_indices = X.indices
        X_indptr = X.indptr
        verbose = 0#2 # XXX : shouldn't verbose be a instance param
        coef_, intercept_ = sgd_fast_sparse.plain_sgd(self.coef_,
                                                      self.intercept_,
                                                      self.loss_function,
                                                      self.penalty_type,
                                                      self.alpha, self.rho,
                                                      X_data,
                                                      X_indices, X_indptr, Y,
                                                      self.n_iter,
                                                      int(self.fit_intercept),
                                                      verbose,
                                                      int(self.shuffle))

        # update self.coef_ and self.sparse_coef_ consistently
        self._set_coef(coef_)
        self.intercept_ = intercept_

        # return self for chaining fit and predict calls
        return self

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples]
           Array containing the predicted class labels (either -1 or 1).
        """
        indices = np.array(self.predict_margin(X) < 0, dtype=np.int)
        return self.classes[indices]

    def predict_margin(self, X):
        """Predict signed 'distance' to the hyperplane (aka confidence score).

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] with signed 'distances' to the hyperplane.
        """
        # np.dot only works correctly if both arguments are sparse matrices
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        return np.ravel(np.dot(self.sparse_coef_, X.T).todense()
                        + self.intercept_)

    def predict_proba(self, X):
        """Predict class membership probability.

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples]
            Contains the membership probabilities of the positive class.
        """
        # how can this be, logisitic *does* implement this
        raise NotImplementedError(
                'sgd does not provide this functionality')
