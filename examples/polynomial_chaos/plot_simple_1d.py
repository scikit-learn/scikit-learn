r"""
===============================================================================
Polynomial Chaos expansion of a one-dimensional model
===============================================================================

This example illustrates how Polynomial Chaos expansions work using a simple
one-dimensional example. Consider

.. math::
    y = x \sin(x)

where the input feature :math:`x \in [0, 10]`.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Let's fix the random seed for reproducibility.
import numpy as np

np.random.seed(2023)


# %%
# First, define the true model :math:`y = x \sin(x)`.
def model(X):
    return (X * np.sin(X)).ravel()


# %%
# Suppose :math:`x` is an input feature that is uniformly distributed between 0
# and 10, i.e., :math:`x\sim\mathcal{U}[0, 10]`. Given a set of training
# samples :math:`X`, we can evaluate the true model to compute the associated
# model outputs :math:`y`.
from scipy.stats import uniform

distribution = uniform(scale=10)
X = distribution.rvs((50, 1))
y = model(X)

# %%
# In reality, only *noisy* model outputs are available, so we add noise
# :math:`\eta\sim\mathcal{N}(0, \sigma)` to the model outputs, with
# :math:`\sigma = 0.5`.
sigma = 0.5
eta = sigma * np.random.randn(*y.shape)
y = y + eta

# %%
# Next, we generate a Polynomial Chaos surrogate model to approximate the model
# output. Since the input follows a uniform distribution, the associated
# orthogonal polynomials are the Legendre polynomials. We construct the
# Polynomial Chaos surrogate for an increasing maximum degree of the orthogonal
# polynomials.
from sklearn.polynomial_chaos import PolynomialChaosExpansion

degrees = range(1, 21)
pces = [PolynomialChaosExpansion(distribution, degree=k) for k in degrees]

# %%
# Next, we fit the Polynomial Chaos expansions to the available data :math:`(X,
# y)`.
for pce in pces:
    pce.fit(X, y)

# %%
# Finally, we plot the various polynomial approximations, together with the
# true model and the data.
import matplotlib.pyplot as plt

plt.scatter(X, y)
X_test = np.linspace(0, 10, 100).reshape(-1, 1)
plt.plot(X_test, model(X_test), label="exact")
for pce in pces[:9:2]:
    y_test = pce.predict(X_test)
    plt.plot(X_test, y_test, label=f"degree {pce.degree}")
plt.legend(bbox_to_anchor=(1, 1), frameon=False, loc="upper left")
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
plt.title("Polynomial Chaos approximations")
plt.tight_layout(pad=1.2)
plt.show()

# %%
# As expected, the approximations improve as the maximum degree of the
# orthogonal polynomials increases. However, if the degree of the orthogonal
# polynomial is chosen too high, we risk overfitting the data.

# %%
# The code below evaluate the training and test error for all Polynomial Chaos
# expansions.
from sklearn.metrics import mean_absolute_error

X_test = np.linspace(0, 10, 100).reshape(-1, 1)
y_test = model(X_test)
training_errors = list()
test_errors = list()
for pce in pces:
    training_errors.append(mean_absolute_error(y, pce.predict(X)))
    test_errors.append(mean_absolute_error(y_test, pce.predict(X_test)))

# %%
# Plotting the test error as a function of the polynomial degree clearly
# illustrates the effect of overfitting.
from matplotlib.ticker import MaxNLocator

plt.plot(degrees, training_errors, marker="o", label="training")
plt.plot(degrees, test_errors, marker="o", label="test")
plt.xlabel("polynomial degree")
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("mean absolute error")
plt.legend(frameon=False)
plt.title("using LinearRegression")
plt.show()

# %%
# As the polynomial degree increases, the training error initially decreases,
# but when the polynomial degree is too large, the test error starts
# increasing, although the training error remains low.

# %%
# Overfitting can be avoided by adding `regularization
# <https://en.wikipedia.org/wiki/Regularization_(mathematics)>`_ when solving
# for the coefficients. In a Polynomial Chaos expansion, we can use, for
# example, `LASSO (least absolute shrinkage and selection operator)
# <https://en.wikipedia.org/wiki/Lasso_(statistics)>`_. This is implemented in
# `scikit-learn` as :class:`~sklearn.linear_model.LassoCV`.
from sklearn.linear_model import LassoCV

estimator = LassoCV(fit_intercept=False, max_iter=100000, tol=1e-1)

# %%
# Let's refit the Polynomimal Chaos expansions using this new estimator.
for pce in pces:
    pce.set_params(estimator=estimator)
    pce.fit(X, y)

# %%
# Next, let's recompute the training and test errors ...
X_test = np.linspace(0, 10, 100).reshape(-1, 1)
y_test = model(X_test)
training_errors = list()
test_errors = list()
for pce in pces:
    training_errors.append(mean_absolute_error(y, pce.predict(X)))
    test_errors.append(mean_absolute_error(y_test, pce.predict(X_test)))

# ... and plot them on a graph.
plt.plot(degrees, training_errors, marker="o", label="training")
plt.plot(degrees, test_errors, marker="o", label="test")
plt.xlabel("polynomial degree")
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("mean absolute error")
plt.legend(frameon=False)
plt.title("using LassoCV")
plt.show()

# %%
# This time, the test error remains low as the degree of the polynomial
# grows.

# %%
#  The
# :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_index_sets.py` example
# illustrates how we can use the :class:`~sklearn.model_selection.GridSearchCV`
# method to help us select the best polynomial degree for a given set of
# measurements.

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_index_sets.py` for an
#     example with more than one feature.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_noisy_data.py`
#     for an example with more details on how to avoid overfitting.
