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

    solver : {'highs-ds', 'highs-ipm', 'highs', 'interior-point', \
            'revised simplex', 'simplex'}, default='interior-point'
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
            solver="interior-point",
            solver_options=None,
    ):
        self.quantile = quantile
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.solver = solver
        self.solver_options = solver_options

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = self._validate_data(X, y, accept_sparse=False,
                                   y_numeric=True, multi_output=False)
        sample_weight = _check_sample_weight(sample_weight, X)

        n_samples, n_features = X.shape
        n_params = n_features

        if self.fit_intercept:
            n_params += 1
            # Note that centering y and X with _preprocess_data does not work
            # for quantile regression.

        # The objective is defined as 1/n * sum(pinball loss) + alpha * L1.
        # So we rescale the penalty term, which is equivalent.
        if self.alpha >= 0:
            alpha = np.sum(sample_weight) * self.alpha
        else:
            raise ValueError(f"Penalty alpha must be a non-negative number, "
                             f"got {self.alpha}".format(self.alpha))

        if self.quantile >= 1.0 or self.quantile <= 0.0:
            raise ValueError(
                f"Quantile should be strictly between 0.0 and 1.0, got "
                f"{self.quantile}")

        if not isinstance(self.fit_intercept, bool):
            raise ValueError(f"The argument fit_intercept must be bool, "
                             f"got {self.fit_intercept}")

        if self.solver not in (
            "highs-ds", "highs-ipm", "highs", "interior-point",
            "revised simplex", "simplex"
        ):
            raise ValueError(f"Invalid value for argument solver, "
                             f"got {self.solver}")
        elif (
            self.solver == "revised simplex"
            and sp_version < parse_version('1.3.0')
        ):
            raise ValueError(f"Solver 'revised simplex' is only available "
                             f"with scipy>=1.3.0, got {sp_version}")
        elif (
            self.solver in ("highs-ds", "highs-ipm", "highs")
            and sp_version < parse_version('1.6.0')
        ):
            raise ValueError(f"Solver {self.solver} is only available "
                             f"with scipy>=1.6.0, got {sp_version}")

        # Use linear programming formulation of quantile regression
        #     min_x c x
        #           A_eq x = b_eq
        #                0 <= x
        # x = (s0, s, t0, t, u, v) = slack variables
        # intercept = s0 + t0
        # coef = s + t
        # c = (alpha * 1_p, alpha * 1_p, quantile * 1_n, (1-quantile) * 1_n)
        # residual = y - X@coef - intercept = u - v
        # A_eq = (1_n, X, -1_n, -X, diag(1_n), -diag(1_n))
        # b_eq = y
        # p = n_features + fit_intercept
        # n = n_samples
        # 1_n = vector of length n with entries equal one
        # see https://stats.stackexchange.com/questions/384909/
        c = np.concatenate([
            np.ones(n_params * 2) * alpha,
            sample_weight * self.quantile,
            sample_weight * (1 - self.quantile),
        ])
        if self.fit_intercept:
            # do not penalize the intercept
            c[0] = 0
            c[n_params] = 0

            A_eq = np.concatenate([
                np.ones((n_samples, 1)),
                X,
                -np.ones((n_samples, 1)),
                -X,
                np.eye(n_samples),
                -np.eye(n_samples),
            ], axis=1)
        else:
            A_eq = np.concatenate([
                X,
                -X,
                np.eye(n_samples),
                -np.eye(n_samples),
            ], axis=1)

        b_eq = y

        result = linprog(
            c=c,
            A_eq=A_eq,
            b_eq=b_eq,
            method=self.solver,
            options=self.solver_options,
        )
        solution = result.x
        if not result.success:
            warnings.warn(
                'Linear programming for Quantile regression did not converge. '
                'Status is {}'.format(result.status), ConvergenceWarning
            )
            if solution is np.nan:
                solution = np.zeros(A_eq.shape[1])

        # positive - negative
        params = solution[:n_params] - solution[n_params:2 * n_params]

        self.n_iter_ = result.nit

        if self.fit_intercept:
            self.coef_ = params[1:]
            self.intercept_ = params[0]
        else:
            self.coef_ = params
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
