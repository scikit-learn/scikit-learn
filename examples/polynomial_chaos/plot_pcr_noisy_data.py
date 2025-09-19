"""
===============================================================================
Polynomial Chaos expansions with noisy data
===============================================================================

This example illustrates how to use Polynomial Chaos regression with noisy
data. The model is taken from [Storlie2009]_.

.. topic: References

    .. [Storlie2009] `C. B. Storlie, L. P. Swiler, J. C. Helton, and C. J.
      Sallaberry. "Implementation and evaluation of nonparametric regression
      procedures for sensitivity analysis of computationally demanding models."
      Reliability Engineering and System Safety, 94(11), 2009.`

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Let's fix the random seed for reproducibility.
import numpy as np

np.random.seed(2023)


# %%
# The Sobol rational function is defined as
def Sobol(x_1, x_2):
    return (x_2 + 0.5) ** 4 / (x_1 + 0.5) ** 2


# %%
# This function is monotonic in each of the inputs. However, there is
# substantial interaction between the inputs :math:`x_1` and :math:`x_2`
# which makes sensitivity analysis difficult. Below is plot of the
# isocontourlines of the function.
import matplotlib.pyplot as plt

_, ax = plt.subplots()
x = np.linspace(0, 1, 128)
X, Y = np.meshgrid(x, x)
Z = Sobol(X, Y)
levels = [0, 0.25, 1, 3, 6, 9, 12, 15, 18, 20]
c = ax.contourf(X, Y, Z, levels=levels, alpha=0.4, antialiased=True)
c = ax.contour(X, Y, Z, levels=levels)
ax.clabel(c, c.levels, inline=True, fontsize=10)
ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$")
plt.show()

# %%
# The inputs :math:`x_1` and :math:`x_2` are assumed to vary according to a
# standard uniform distribution.
from scipy.stats import uniform

distribution = uniform()

# %%
# Now, let's generate some samples from the input distribution, and
# evaluate the polynomial model.
X = distribution.rvs((56, 2))
y = Sobol(X[:, 0], X[:, 1])

# %%
# Let's approximate the polynomial using an 8th order Polynomial Chaos
# expansion.
from sklearn.polynomial_chaos import PolynomialChaosRegressor

distribution = uniform()
degree = 8
pce = PolynomialChaosRegressor(distribution, degree=degree)

# %%
# Next, we perform the Polynomial Chaos regression by fitting the model to the
# data :math:`(X, y)`
pce.fit(X, y)

# %%
# Let's also generate some test data.
X_test = distribution.rvs((1000, 2))
y_test = Sobol(X_test[:, 0], X_test[:, 1])

# %%
# We can visualize the goodness-of-fit for the test data as follows.
from sklearn.metrics import PredictionErrorDisplay

PredictionErrorDisplay.from_predictions(
    y_true=y_test, y_pred=pce.predict(X_test), kind="actual_vs_predicted"
)
plt.show()

# %%
# Let's compute the mean absolute error for training and test data
from sklearn.metrics import mean_absolute_error

print("training error:", mean_absolute_error(y, pce.predict(X)))
print("test error:", mean_absolute_error(y_test, pce.predict(X_test)))

# %%
# We conclude that our Polynomial Chaos surrogate is reasonably accurate, and
# we can trust the values for the main sensitivity indices extracted from the
# surrogate.
pce.main_sens()

# %%
# Now, supposed that the exact model evaluations are not accessible, and that
# we can only access noisy approximations of the model output.
noise = 20 * 0.01  # 1% noise
y_noisy = y + noise * np.random.randn(*y.shape)

# %%
# Let's fit our :class:`~sklearn.polynomial_chaos.PolynomialChaosRegressor` to
# the noisy data instead.
pce.fit(X, y_noisy)

# %%
# We can again visualize the goodness-of-fit for the test data, shown below.
PredictionErrorDisplay.from_predictions(
    y_true=y_test, y_pred=pce.predict(X_test), kind="actual_vs_predicted"
)
plt.show()

# %%
# When only noisy outputs are available, we are unable to construct a good
# surrogate model. Let's compute the mean absolute error for training and test
# data.
print("training error:", mean_absolute_error(y, pce.predict(X)))
print("test error:", mean_absolute_error(y_test, pce.predict(X_test)))

# %%
# We can somewhat alleviate this issue by using a sparsity-promoting estimator.
from sklearn.linear_model import LassoCV

estimator = LassoCV(
    fit_intercept=False, alphas=np.logspace(-12, 2, 25), max_iter=500000
)
pce.set_params(estimator=estimator)
pce.fit(X, y_noisy)
PredictionErrorDisplay.from_predictions(
    y_true=y_test, y_pred=pce.predict(X_test), kind="actual_vs_predicted"
)
plt.show()

# %%
# Let's again compute the mean absolute error for training and test
# data.
print("training error:", mean_absolute_error(y, pce.predict(X)))
print("test error:", mean_absolute_error(y_test, pce.predict(X_test)))

# %%
# Using the sparse estimator, training and test data are much smaller. They are
# also of similar magnitude, indicating that our surrogate model is not
# overfiting. Hence, despite the noisy outputs, we are able to construct a
# reasonably accurate surrogate.

# %%
# We also find values for the sensivitiy indices that are reasonably close to
# the ones found using the noiseless outputs.
pce.main_sens()

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_ishigami.py` for
#     a global sensitivity analysis of the Ishigami function, a
#     well-known test problem.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_sobol_g.py` for
#     an example of how to adaptively construct the multiindex set in the
#     Polynomial Chaos expansion to compute sensitivity indices.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_feature_selection_g.py`
#     for an example of how to use pruning to remove small basis terms from
#     the expansion.
