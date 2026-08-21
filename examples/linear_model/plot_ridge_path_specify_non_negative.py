"""
================================================================================================
Plot Ridge coefficients as a function of the regularization
and positive coefficient constraints
================================================================================================

Shows the effect of collinearity in the coefficients of an estimator.
As well as a comparison between options to force positive coefficients.

.. currentmodule:: sklearn.linear_model

:class:`Ridge` Regression is the estimator used in this example.
Each color represents a different feature of the
coefficient vector, and this is displayed as a function of the
regularization parameter.

This example is an extension of the `plot_ridge_path` example. We
add the option to force some coefficients to be positive. This is
implemented by either setting the `positive` parameter of the Ridge estimator
to `True` or by setting the `positive` parameter to an array containing
`True` for the features that should be positive and `False` for the others.

The original example also shows the usefulness of applying Ridge regression
to highly ill-conditioned matrices. For such matrices, a slight
change in the target variable can cause huge variances in the
calculated weights. In such cases, it is useful to set a certain
regularization (alpha) to reduce this variation (noise).

This example expands on the alpha sensitivity exploration. We show the
effect of different levels of allowed minimum alpha values on the
coefficients.

We show at relatively high alpha values, the coefficients are similar
between the two Ridge solvers, and the coefficients are relatively
stable in relation to the alpha value. This is true for both positive
and non-positive coefficient constraints.

However, at very low alpha values, the magnitude of the coefficients
is very different between the two Ridge solvers. Additionally, in the
case of positive coefficient constraints, the coefficient weights
have a period of instability as alpha decreases - but eventually
converge to more stable values.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from sklearn import linear_model

# Set random seed for reproducibility
np.random.seed(123)

# X is the 10x10 Hilbert matrix
X = 1.0 / (np.arange(1, 11) + np.arange(0, 10)[:, np.newaxis])
y = np.ones(10)

# %%
# Positive coefficient constraints
# Assign a True boolean value to any feature index to enforce a positive
# coefficient. In this case, the first, third, fifth, seventh, and ninth
# features will be forced to have positive coefficients.
# Coefficient array must have same number of values and X.shape[1]
# -------------------------------
X_positive_coefficient = [
    True,
    False,
    True,
    False,
    True,
    False,
    True,
    False,
    True,
    False,
]


# %%
# Construct function using lbfgs solver without positive coefficient constraints
# Ridge implementation of lbfgs solver requires positive=True or positive=[bool,...]
# This is because other solvers are far more efficient for this problem.
# Here we have an example, both to show how positive constraints are implemented
# and to show the effect of regularization on the coefficients.
# -------------------------------


def _solve_lbfgs_unconstrained(
    X,
    y,
    alpha,
    tol=1e-4,
):
    # replicate transformations in _BaseRidge
    alpha = np.array([alpha] * y.shape[0])
    y = np.reshape(y, (-1, 1))

    # get n_features for coefs shape
    n_samples, n_features = X.shape

    coefs = np.empty((y.shape[1], n_features), dtype=X.dtype)

    # scipy.optimize config options
    config = {
        "method": "L-BFGS-B",
        "tol": tol,
        "jac": True,
    }

    # Optimize objective function
    for i in range(y.shape[1]):
        x0 = np.zeros((n_features,))
        y_column = y[:, i]

        def func(w):
            residual = X.dot(w) - y_column
            f = 0.5 * residual.dot(residual) + 0.5 * alpha[i] * w.dot(w)
            grad = X.T @ residual + alpha[i] * w
            return f, grad

        result = optimize.minimize(func, x0, **config)
        coefs[i] = result["x"]

    return coefs


# %%
# Generate coefficients for Ridge regression with different alphas
# for different scenarios
# 1. Default Ridge regression, no positive constraint
# 2. Ridge regression with lbfgs solver, no positive constraint
# 3. Ridge regression with all positive constraints
# 4. Ridge regression with some positive constraints
# -------------


def generate_path(X, y, alphas, positive=False, scenario="default"):
    coefs = []
    for a in alphas:
        if scenario == "default":
            ridge = linear_model.Ridge(alpha=a, fit_intercept=False, positive=positive)
            ridge.fit(X, y)
            coefs.append(ridge.coef_)
        elif scenario == "lbfgs_unconstrained":
            ridge = _solve_lbfgs_unconstrained(X, y, alpha=a)
            coefs.append(ridge.flatten())
    return coefs


def plot_results(coefs, alphas, title):
    """Plot the coefficients as a function of alpha"""
    ax = plt.gca()
    ax.plot(alphas, coefs)
    ax.set_xscale("log")
    ax.set_xlim(ax.get_xlim()[::-1])  # reverse axis
    plt.xlabel("alpha")
    plt.ylabel("weights")
    plt.title(title)
    plt.axis("tight")
    plt.show()


# %%
# Generate coefficients for Ridge regression with different alphas
# for different levels of allowed minimum alpha
# Default Ridge regression, with no positive constraint and lbfgs solver
# Generalize to similar coefficeint weights at high alpha values
# However, at very low alpha values, the coefficients are very different
# both between the two Ridge solvers and between positive and non-positive
# coefficient constraints.
# -------------

# Min alpha 1e-3, max alpha 1
n_alphas = 200
alphas = np.logspace(-3, 0, n_alphas)

coefs = generate_path(X, y, alphas, positive=False, scenario="default")
coefs_lbfgs_unconstrained = generate_path(
    X, y, alphas, positive=False, scenario="lbfgs_unconstrained"
)
coefs_all_positive = generate_path(X, y, alphas, positive=True, scenario="default")
coefs_some_positive = generate_path(
    X, y, alphas, positive=X_positive_coefficient, scenario="default"
)


min_alpha = alphas.min()
plot_results(
    coefs,
    alphas,
    f"Ridge coefs without positive coefficient constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_lbfgs_unconstrained,
    alphas,
    f"Approx. of lbfgs Ridge without positive constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_some_positive,
    alphas,
    f"Ridge coefs with conditional positive constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_all_positive,
    alphas,
    f"Ridge coefs with all positive constraints - alpha_min={min_alpha}",
)


# %%
# Min alpha 1e-5, max alpha 1
n_alphas = 200
alphas = np.logspace(-5, 0, n_alphas)

coefs = generate_path(X, y, alphas, positive=False, scenario="default")
coefs_lbfgs_unconstrained = generate_path(
    X, y, alphas, positive=False, scenario="lbfgs_unconstrained"
)
coefs_all_positive = generate_path(X, y, alphas, positive=True, scenario="default")
coefs_some_positive = generate_path(
    X, y, alphas, positive=X_positive_coefficient, scenario="default"
)

min_alpha = alphas.min()
plot_results(
    coefs,
    alphas,
    f"Ridge coefs without positive coefficient constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_lbfgs_unconstrained,
    alphas,
    f"Approx. of lbfgs Ridge without positive constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_some_positive,
    alphas,
    f"Ridge coefs with conditional positive constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_all_positive,
    alphas,
    f"Ridge coefs with all positive constraints - alpha_min={min_alpha}",
)


# %%
# Min alpha 1e-10, max alpha 1
n_alphas = 200
alphas = np.logspace(-10, 0, n_alphas)

coefs = generate_path(X, y, alphas, positive=False, scenario="default")
coefs_lbfgs_unconstrained = generate_path(
    X, y, alphas, positive=False, scenario="lbfgs_unconstrained"
)
coefs_all_positive = generate_path(X, y, alphas, positive=True, scenario="default")
coefs_some_positive = generate_path(
    X, y, alphas, positive=X_positive_coefficient, scenario="default"
)

min_alpha = alphas.min()
plot_results(
    coefs,
    alphas,
    f"Ridge coefs without positive coefficient constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_lbfgs_unconstrained,
    alphas,
    f"Approx. of lbfgs Ridge without positive constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_some_positive,
    alphas,
    f"Ridge coefs with conditional positive constraints - alpha_min={min_alpha}",
)
plot_results(
    coefs_all_positive,
    alphas,
    f"Ridge coefs with all positive constraints - alpha_min={min_alpha}",
)
# %%
