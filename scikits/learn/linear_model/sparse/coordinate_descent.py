# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.
"""Implementation of coordinate descent for the Elastic Net with sparse data.
"""

import warnings
import numpy as np
from scipy import sparse

from ..base import LinearModel
from . import cd_fast_sparse


class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer

    This implementation works on scipy.sparse X and dense coef_.

    rho=1 is the lasso penalty. Currently, rho <= 0.01 is not
    reliable, unless you supply your own sequence of alpha.

    Parameters
    ----------
    alpha : float
        Constant that multiplies the L1 term. Defaults to 1.0
    rho : float
        The ElasticNet mixing parameter, with 0 < rho <= 1.
    coef_ : ndarray of shape n_features
        The initial coeffients to warm-start the optimization
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

        TODO: fit_intercept=True is not yet implemented
    """

    def __init__(self, alpha=1.0, rho=0.5, fit_intercept=False):
        if fit_intercept:
            raise NotImplementedError("fit_intercept=True is not implemented")
        self.alpha = alpha
        self.rho = rho
        self.fit_intercept = fit_intercept
        self.intercept_ = 0.0
        self._set_coef(None)

    def _set_coef(self, coef_):
        self.coef_ = coef_
        if coef_ is None:
            self.sparse_coef_ = None
        else:
            # sparse representation of the fitted coef for the predict method
            self.sparse_coef_ = sparse.csr_matrix(coef_)

    def fit(self, X, y, max_iter=1000, tol=1e-4, **params):
        """Fit current model with coordinate descent

        X is expected to be a sparse matrix. For maximum efficiency, use a
        sparse matrix in CSC format (scipy.sparse.csc_matrix)
        """
        self._set_params(**params)
        X = sparse.csc_matrix(X)
        y = np.asanyarray(y, dtype=np.float64)

        # NOTE: we are explicitly not centering the data the naive way to
        # avoid breaking the sparsity of X

        n_samples, n_features = X.shape[0], X.shape[1]
        if self.coef_ is None:
            self.coef_ = np.zeros(n_features, dtype=np.float64)

        alpha = self.alpha * self.rho * n_samples
        beta = self.alpha * (1.0 - self.rho) * n_samples
        X_data = np.array(X.data, np.float64)

        # TODO: add support for non centered data
        coef_, self.dual_gap_, self.eps_ = \
                cd_fast_sparse.enet_coordinate_descent(
                    self.coef_, alpha, beta, X_data, X.indices, X.indptr, y,
                    max_iter, tol)

        # update self.coef_ and self.sparse_coef_ consistently
        self._set_coef(coef_)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want'
                                'to increase the number of interations')

        # XXX TODO: implement intercept_ fitting

        # return self for chaining fit and predict calls
        return self

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] with the predicted real values
        """
        # np.dot only works correctly if both arguments are sparse matrices
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        return np.ravel(np.dot(self.sparse_coef_, X.T).todense()
                        + self.intercept_)


class Lasso(ElasticNet):
    """Linear Model trained with L1 prior as regularizer

    This implementation works on scipy.sparse X and dense coef_. Technically
    this is the same as Elastic Net with the L2 penalty set to zero.

    Parameters
    ----------
    alpha : float
        Constant that multiplies the L1 term. Defaults to 1.0
    coef_ : ndarray of shape n_features
        The initial coeffients to warm-start the optimization
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    """

    def __init__(self, alpha=1.0, fit_intercept=False):
        super(Lasso, self).__init__(alpha=alpha, rho=1.0,
                                    fit_intercept=fit_intercept)
