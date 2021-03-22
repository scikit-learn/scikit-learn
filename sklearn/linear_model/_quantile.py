# Authors: David Dale dale.david@mail.ru
# License: BSD 3 clause
import warnings

import numpy as np
from scipy.optimize import linprog

from ..base import BaseEstimator, RegressorMixin
from ._base import LinearModel
from ..exceptions import ConvergenceWarning
from ..utils.validation import _check_sample_weight
from ..utils.fixes import sp_version, parse_version


class QuantileRegressor(LinearModel, RegressorMixin, BaseEstimator):
    """Linear regression model that predicts conditional quantiles.

    The linear Quantile Regressor optimizes the pinball loss for a
    desired `quantile` and is robust to outliers.

    Read more in the :ref:`User Guide <quantile_regression>`

    .. versionadded:: 1.0

    Parameters
    ----------
    quantile : float, default=0.5
        The quantile that the model tries to predicts. It must be strictly
        between 0 and 1.

    alpha : float, default=1e-4
        Constant that multiplies L1 penalty term.

    fit_intercept : bool, default=True
        Whether or not to fit the intercept. This can be set to False
        if the data is already centered around the origin.

    copy_X : boolean, default=True
        If True, X will be copied; else, it may be overwritten.

    solver : str, optional, default 'auto'
        Name of the solver used by scipy.optimize.linprog.
        If it is 'auto', will use 'highs' with scipy>=1.6.0
        and 'interior-point' with older versions.

    solver_options : dict, default=None
        Additional parameters passed to scipy.optimize.linprog as options.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the features.

    intercept_ : float
        The intercept of the model, aka bias term.

    n_iter_ : int
        The actual number of iterations performed by the solver.

    References
    ----------
    .. [1] Koenker, R., & Bassett Jr, G. (1978). Regression quantiles.
            Econometrica: journal of the Econometric Society, 33-50.

    .. [2] Chen, C., & Wei, Y. (2005).
           Computational issues for quantile regression.
           Sankhya: The Indian Journal of Statistics, 399-417.
    """

    def __init__(
            self,
            quantile=0.5,
            alpha=0.0001,
            fit_intercept=True,
            copy_X=True,
            solver='auto',
            solver_options=None,
    ):
        self.quantile = quantile
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.solver = solver
        self.solver_options = solver_options

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        y : array-like, shape (n_samples,)
            Target vector relative to X.

        sample_weight : array-like, shape (n_samples,)
            Weight given to each sample.

        Returns
        -------
        self : object
            Returns self.
        """

        X, y = self._validate_data(X, y, accept_sparse=False,
                                   y_numeric=True, multi_output=False)

        sample_weight = _check_sample_weight(sample_weight, X)

        X, y, X_offset, y_offset, X_scale = self._preprocess_data(
            X, y, fit_intercept=self.fit_intercept, copy=self.copy_X,
            sample_weight=sample_weight)

        if self.quantile >= 1.0 or self.quantile <= 0.0:
            raise ValueError(
                "Quantile should be strictly between 0.0 and 1.0, got %f"
                % self.quantile)

        n_obs, n_slopes = X.shape
        n_params = n_slopes

        X_full = X
        if self.fit_intercept:
            n_params += 1
            X_full = np.concatenate([np.ones([n_obs, 1]), X], axis=1)

        # the linear programming formulation of quantile regression
        # follows https://stats.stackexchange.com/questions/384909/
        c_vector = np.concatenate([
            np.ones(n_params * 2) * self.alpha,
            sample_weight * self.quantile,
            sample_weight * (1 - self.quantile),
        ])
        # do not penalize the intercept
        if self.fit_intercept:
            c_vector[0] = 0
            c_vector[n_params] = 0

        a_eq_matrix = np.concatenate([
            X_full,
            -X_full,
            np.eye(n_obs),
            -np.eye(n_obs),
        ], axis=1)
        b_eq_vector = y

        method = self.solver
        if method == 'auto':
            if sp_version < parse_version('1.0.0'):
                method = 'simplex'
            elif sp_version < parse_version('1.6.0'):
                method = 'interior-point'
            else:
                method = 'highs'

        result = linprog(
            c=c_vector,
            A_eq=a_eq_matrix,
            b_eq=b_eq_vector,
            method=method,
            options=self.solver_options,
        )
        solution = result.x
        if not result.success:
            warnings.warn(
                'Linear programming for Quantile regression did not converge. '
                'Status is {}'.format(result.status), ConvergenceWarning
            )
            if solution is np.nan:
                solution = np.zeros(a_eq_matrix.shape[1])

        params_pos = solution[:n_params]
        params_neg = solution[n_params:2 * n_params]
        params = params_pos - params_neg

        self.n_iter_ = result.nit

        self.coef_ = params[self.fit_intercept:]
        # do not use self.set_intercept_, because it assumes intercept is zero
        # if the data is normalized, which is false in this case
        if self.fit_intercept:
            self.coef_ = self.coef_ / X_scale
            self.intercept_ = params[0] + y_offset \
                - np.dot(X_offset, self.coef_.T)
        else:
            self.intercept_ = 0.0
        return self

    def _more_tags(self):
        if sp_version >= parse_version('1.0.0'):
            return {}
        return {
            '_xfail_checks': {
                'check_regressors_train':
                'scipy.optimize.linprog is unstable in versions before 1.0.0',
                'check_regressor_data_not_an_array':
                'scipy.optimize.linprog is unstable in versions before 1.0.0',
            }
        }
