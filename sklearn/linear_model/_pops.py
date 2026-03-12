"""
POPS (Pointwise Optimal Parameter Sets) Regression.

Bayesian regression for low-noise data accounting for model misspecification.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Real

import numpy as np
from scipy.linalg import eigh
from scipy.stats import qmc

from sklearn.base import _fit_context
from sklearn.linear_model._base import _preprocess_data
from sklearn.linear_model._bayes import BayesianRidge
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.validation import (
    _check_sample_weight,
    check_is_fitted,
    validate_data,
)


class POPSRegression(BayesianRidge):
    _parameter_constraints: dict = {
        **BayesianRidge._parameter_constraints,
        "mode_threshold": [Interval(Real, 0, None, closed="neither")],
        "resample_density": [Interval(Real, 0, None, closed="neither")],
        "resampling_method": [StrOptions({"uniform", "sobol", "latin", "halton"})],
        "percentile_clipping": [Interval(Real, 0, 50.0, closed="both")],
        "leverage_percentile": [Interval(Real, 0.0, 100.0, closed="left")],
        "posterior": [StrOptions({"hypercube", "ensemble"})],
    }

    def __init__(
        self,
        *,
        max_iter=300,
        tol=1.0e-3,
        alpha_1=1.0e-6,
        alpha_2=1.0e-6,
        lambda_1=1.0e-6,
        lambda_2=1.0e-6,
        alpha_init=None,
        lambda_init=None,
        compute_score=False,
        fit_intercept=False,
        copy_X=True,
        verbose=False,
        mode_threshold=1.0e-8,
        resample_density=1.0,
        resampling_method="uniform",
        percentile_clipping=0.0,
        leverage_percentile=50.0,
        posterior="hypercube",
    ):
        super().__init__(
            max_iter=max_iter,
            tol=tol,
            alpha_1=alpha_1,
            alpha_2=alpha_2,
            lambda_1=lambda_1,
            lambda_2=lambda_2,
            alpha_init=alpha_init,
            lambda_init=lambda_init,
            compute_score=compute_score,
            fit_intercept=fit_intercept,
            copy_X=copy_X,
            verbose=verbose,
        )
        self.mode_threshold = mode_threshold
        self.resample_density = resample_density
        self.resampling_method = resampling_method
        self.percentile_clipping = percentile_clipping
        self.leverage_percentile = leverage_percentile
        self.posterior = posterior

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        pops_fit_intercept = self.fit_intercept
        if self.fit_intercept:
            X = np.asarray(X)
            X = np.hstack([X, np.ones((X.shape[0], 1))])
            self.fit_intercept = False
        try:
            super().fit(X, y, sample_weight=sample_weight)
            X_pp, y_pp = validate_data(
                self, X, y, dtype=[np.float64, np.float32], reset=False
            )
            if sample_weight is not None:
                sw = _check_sample_weight(sample_weight, X_pp, dtype=X_pp.dtype)
            else:
                sw = None
            X_pp, y_pp, _, _, _, _ = _preprocess_data(
                X_pp,
                y_pp,
                fit_intercept=False,
                copy=True,
                sample_weight=sw,
                rescale_with_sw=True,
            )
            scaled_sigma_ = self.alpha_ * self.sigma_
            errors = y_pp - X_pp @ self.coef_
            pointwise_correction = np.dot(X_pp, scaled_sigma_)
            leverage_scores = np.sum(pointwise_correction * X_pp, axis=1)
            safe_leverage = np.where(leverage_scores > 0, leverage_scores, np.inf)
            pointwise_correction *= (errors / safe_leverage)[:, None]
            leverage_mask = leverage_scores >= np.percentile(
                leverage_scores, self.leverage_percentile
            )
            if not np.any(leverage_mask):
                leverage_mask = np.ones(X_pp.shape[0], dtype=bool)
            self._pointwise_correction = pointwise_correction
            self._leverage_scores = leverage_scores
            self._leverage_mask = leverage_mask
            self.posterior_samples_, self.misspecification_sigma_ = (
                self._build_posterior()
            )
            self._fitted_with_intercept = pops_fit_intercept
        finally:
            self.fit_intercept = pops_fit_intercept
        return self

    def _build_posterior(self):
        pc = self._pointwise_correction[self._leverage_mask]
        if self.posterior == "ensemble":
            sigma = pc.T @ pc / pc.shape[0]
            return pc.T, sigma
        elif self.posterior == "hypercube":
            self._hypercube_support, self._hypercube_bounds = self._fit_hypercube(pc)
            return self._sample_hypercube()

    def _fit_hypercube(self, pc):
        e_values, e_vectors = eigh(pc.T @ pc)
        mask = e_values > self.mode_threshold * e_values.max()
        e_vectors = e_vectors[:, mask]
        projections = e_vectors.copy()
        projected = pc @ projections
        bounds = [
            np.percentile(projected, self.percentile_clipping, axis=0),
            np.percentile(projected, 100.0 - self.percentile_clipping, axis=0),
        ]
        return projections, bounds

    def _sample_hypercube(self, size=None, resampling_method=None):
        if resampling_method is None:
            resampling_method = self.resampling_method
        low, high = self._hypercube_bounds
        n_resample = (
            size
            if size is not None
            else max(100, int(self.resample_density * self._leverage_scores.size))
        )
        if resampling_method == "latin":
            samples = qmc.LatinHypercube(d=low.size).random(n_resample).T
        elif resampling_method == "sobol":
            n_resample = 2 ** int(np.log2(n_resample))
            samples = qmc.Sobol(d=low.size).random(n_resample).T
        elif resampling_method == "halton":
            samples = qmc.Halton(d=low.size).random(n_resample).T
        else:
            samples = np.random.uniform(size=(low.size, n_resample))
        samples = low[:, None] + (high - low)[:, None] * samples
        hypercube_samples = self._hypercube_support @ samples
        hypercube_sigma = (
            hypercube_samples @ hypercube_samples.T / hypercube_samples.shape[1]
        )
        return hypercube_samples, hypercube_sigma

    def predict(
        self, X, return_std=False, return_bounds=False, return_epistemic_std=False
    ):
        check_is_fitted(self)
        if getattr(self, "_fitted_with_intercept", False):
            X = np.asarray(X)
            X = np.hstack([X, np.ones((X.shape[0], 1))])
        y_mean = self._decision_function(X)
        result = [y_mean]
        if return_std or return_bounds or return_epistemic_std:
            y_epistemic_var = (np.dot(X, self.sigma_) * X).sum(axis=1)
            if return_std:
                y_misspecification_var = (
                    np.dot(X, self.misspecification_sigma_) * X
                ).sum(axis=1)
                result.append(
                    np.sqrt(
                        y_misspecification_var + y_epistemic_var + 1.0 / self.alpha_
                    )
                )
            if return_bounds:
                y_posterior = X @ self.posterior_samples_
                result.extend(
                    [y_posterior.max(axis=1) + y_mean, y_posterior.min(axis=1) + y_mean]
                )
            if return_epistemic_std:
                result.append(np.sqrt(y_epistemic_var))
        return result[0] if len(result) == 1 else tuple(result)
