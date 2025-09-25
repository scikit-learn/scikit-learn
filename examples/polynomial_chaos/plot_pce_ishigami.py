r"""
===============================================================================
Global sensitivity analysis of the Ishigami function
===============================================================================

This example illustrates how to use a Polynomial Chaos Expansion to perform a
Global Sensitivity Analysis of the Ishigami function. The Ishigami function is
a well-known test problem in uncertainty quantification and sensitivity
analysis [Saltelli2000]_. It is defined as

.. math::
    y = \sin x_1 + a \sin^2 x_2 + b x_3^4 \sin x_1

The inputs are usually assumed to be independent and identically distributed,
following a uniform distribution on :math:`[-\pi, \pi]`. The parameters are
usually chosen as :math:`a = 7` and :math:`b = 0.1`.

.. topic: References

    .. [Saltelli2000] `Saltelli, Andrea, et al. "Global sensitivity analysis:
       the primer." John Wiley & Sons, 2008.`

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Let's fix the random seed for reproducibility.
import numpy as np

np.random.seed(2023)

# %%
# First, let's define the model.
a = 7
b = 0.1


def ishigami(X):
    y = np.sin(X[:, 0]) + a * np.sin(X[:, 1]) ** 2 + b * X[:, 2] ** 4 * np.sin(X[:, 0])
    return y


# %%
# Next, let's generate some input/output data. If we store the input features
# in a `DataFrame`, we can easily access them by name later.
from pandas import DataFrame
from scipy.stats import uniform

distribution = uniform(loc=-np.pi, scale=2 * np.pi)
X = distribution.rvs((660, 3))
y = ishigami(X)
X = DataFrame(data=X, columns=("$x_1$", "$x_2$", "$x_3$"))

# %%
# Define a Polynomial Chaos expansion and fit the data. We use a 9th order
# total degree polynomial basis. The number of basis terms :math:`P` can be
# computed as
#
# .. math::
#   P = \frac{(d + k)!}{d!k!}
#
# where :math:`d` is the dimension (i.e., the number of input features), and
# :math:`k` is the degree of the polynomial. In our case, we have :math:`d = 3`
# and :math:`k = 9`, so there are :math:`P = 220` basis terms.
from sklearn.polynomial_chaos import PolynomialChaosExpansion

pce = PolynomialChaosExpansion(distribution, degree=9)
pce.fit(X, y)

# %%
# Next, let's generate some test data.
X_test = distribution.rvs((660, 3))
y_test = ishigami(X_test)
X_test = DataFrame(data=X_test, columns=("$x_1$", "$x_2$", "$x_3$"))

# %%
# We can visualize the goodness-of-fit as follows.
import matplotlib.pyplot as plt

from sklearn.metrics import PredictionErrorDisplay

PredictionErrorDisplay.from_predictions(
    y_true=y_test, y_pred=pce.predict(X_test), kind="actual_vs_predicted"
)
plt.show()

# %%
# It seems like a polynomial basis with polynomials of degree 9 is sufficient
# to obtain a reasonably accurate surrogate model. Hence, we expect the
# sensitivity information, extracted from the Polynomial Chaos surrogate, to
# be reasonably accurate as well.

# %%
# After the expansion is fitted, we can easily extract estiamtes for statistics
# (such as mean and variance) of the model output. The exact values are
#
# .. math::
#       \mathbb{E}[y] = \frac{a}{2} = 3.5
#
# and
#
# .. math::
#   \mathbb{V}[y] = \frac{a^2}{8} + \frac{b\pi^4}{5} + \frac{b^2\pi^8}{18} +
#       \frac{1}{2} \approx 13.844588
#
mean = a / 2
var = a**2 / 8 + b * np.pi**4 / 5 + b**2 * np.pi**8 / 18 + 1 / 2
DataFrame(
    {"predicted": [pce.mean(), pce.var()], "exact": [mean, var]},
    index=["mean", "variance"],
)

# %%
# We can also extract the main-effect sensitivity indices from the Polynomial
# Chaos expansion. The :math:`j`\ th main-effect index indicates the effect on
# the model output of varying the :math:`j`\ th parameter alone. The exact
# values are
#
# .. math::
#   S_1 &= \left(\frac{b \pi^4}{5} + \frac{b^2 \pi^8}{50} + 1/2 \right)
#          \bigg/ \mathbb{V}[y] \approx 0.313905 \\
#   S_2 &= \frac{a^2}{8} \bigg/ \mathbb{V}[y] \approx 0.442411 \\
#   S_3 &= 0
#
S1 = (b * np.pi**4 / 5 + b**2 * np.pi**8 / 50 + 1 / 2) / var
S2 = (a**2 / 8) / var
S3 = 0
df = DataFrame(
    {"predicted": pce.main_sens().flatten(), "exact": [S1, S2, S3]},
    index=["S1", "S2", "S3"],
)
df["predicted"].plot.bar()
plt.ylabel("main sensitivity indices")
df

# %%
# Similarly, we can extract the total-effect sensitivity indices from the
# Polynomial Chaos expansion. The :math:`j`\ th total-effect index indicates
# the effect on the model output of varying the :math:`j`\ th parameter,
# including interaction effects. The exact values are
#
# .. math::
#    S_1^T &= \left( \frac{1}{2} \left(1 + {b\pi^4}{5}\right)^2 +
#          \frac{8b^2\pi^8}{225} \right) \bigg/ \mathbb{V}[y] \approx
#          0.557589 \\
#    S_2^T &= \frac{a^2}{8} \bigg/ \mathbb{V}[y] \approx 0.442411 \\
#    S_3^T &= \left(\frac{8b^2\pi^8}{225}\right) \bigg/ \mathbb{V}[y]
#          \approx 0.442411 \\
#
S1_t = (1 / 2 * (1 + b * np.pi**4 / 5) ** 2 + 8 * b**2 * np.pi**8 / 225) / var
S2_t = (a**2 / 8) / var
S3_t = (8 * b**2 * np.pi**8 / 225) / var
df = DataFrame(
    {"predicted": pce.total_sens().flatten(), "exact": [S1_t, S2_t, S3_t]},
    index=["S1_t", "S2_t", "S3_t"],
)
df["predicted"].plot.bar()
plt.ylabel("total sensitivity indices")
df

# %%
# Only one of the joint sensitivity indices, :math:`S_{13}`, is nonzero. It's
# exact value is
#
# .. math::
#   S_{13} = \left(\frac{8b^2\pi^8}{225}\right) \bigg/ \mathbb{V}[y]
#
S12 = 0
S13 = 8 * b**2 * np.pi**8 / 225 / var
S23 = 0
joint = list()
joint.append(pce.joint_sens(0, 1))
joint.append(pce.joint_sens(0, 2))
joint.append(pce.joint_sens(1, 2))
DataFrame({"predicted": joint, "exact": [S12, S13, S23]}, index=["S12", "S13", "S23"])

# %%
# Convergence of the sensitivity indices
# --------------------------------------
# Now, we investigate how the values of the (total) sensitivity indices
# change as a function of the degree of the Polynomial Chaos expansion. For a
# given degree, we set the number of learning points as
#
# .. math::
#   N = d \times \frac{(d + k)!}{d!k!}
#
# where :math:`d` is the dimension (i.e., the number of input features), and
# :math:`k` is the degree of the polynomial. This is to avoid an
# underdetermined system when computing the values of the coefficients in the
# expansion.
from math import factorial

total_sens = list()
n_samples = list()
for degree in range(3, 12):
    n = 3 * factorial(3 + degree) // (6 * factorial(degree))
    X = distribution.rvs((n, 3))
    y = ishigami(X)
    pce = PolynomialChaosExpansion(distribution, degree=degree)
    pce.fit(X, y)
    total_sens.append(pce.total_sens())
    n_samples.append(n)
total_sens = np.vstack(total_sens).T

# %%
# Next, we plot the values of the sensitivity indices as function of the number
# of learning points. We observe that the sensitivity indices converge rapidly
# to the exact values.
for j, sens in enumerate([S1_t, S2_t, S3_t]):
    plt.plot(
        n_samples,
        total_sens[j],
        label=f"$S_{j + 1}$",
        color=f"C{j}",
        marker="o",
    )
    plt.plot(
        n_samples,
        np.full(len(n_samples), sens),
        color=f"C{j}",
        linestyle="dashed",
    )
plt.xlabel("number of training points")
plt.ylabel("sensitivity indices")
plt.title("Global sensitivity analysis (Polynomial Chaos)")
plt.ylim(0, 1)
plt.legend(frameon=False)
plt.show()

# %%
# Thus, for a polynomial basis with polynomials of degree > 6, using
# 504 training samples, the computed values of the sensitivity indices are
# almost indistinguishable from the exact values.

# %%
# We can compare our results to a sampling-based approach to obtain the Sobol
# sensitivity indices, implemented in `scipy`. This is approach also known as
# the "pick-and-freeze" approach. We evaluate the total-effect sensitivity
# indices for an increasing number of samples :math:`n` (but note that the
# total number of model evaluations is :math:`n \times (d + 2)`, where, for the
# Ishigami function, the dimension :math:`d = 3`. The method only works if
# :math:`n` is a power of 2.
try:
    from scipy.stats import sobol_indices

    scipy_has_sobol_indices = True
except ImportError:
    scipy_has_sobol_indices = False

if scipy_has_sobol_indices:
    from scipy.stats import sobol_indices

    def ishigami(x):
        return np.sin(x[0]) + a * np.sin(x[1]) ** 2 + b * (x[2] ** 4) * np.sin(x[0])

    total_sens_sampling = list()
    n_samples_sampling = list()
    for log_n in range(2, 12):
        n = 2**log_n
        indices = sobol_indices(func=ishigami, n=n, dists=[distribution] * 3)
        total_sens_sampling.append(indices.total_order)
        n_samples_sampling.append(n * 5)  # 5 is dimension + 2
    total_sens_sampling = np.vstack(total_sens_sampling).T

    # %%
    # Now we recreate the same figure as before.
    for j, sens in enumerate([S1_t, S2_t, S3_t]):
        plt.plot(
            n_samples_sampling,
            total_sens_sampling[j],
            label=f"$S_{j + 1}$",
            color=f"C{j}",
            marker="o",
        )
        plt.plot(
            n_samples_sampling,
            np.full(len(n_samples_sampling), sens),
            color=f"C{j}",
            linestyle="dashed",
        )
    plt.xlabel("number of training points")
    plt.ylabel("sensitivity indices")
    plt.title("Global sensitivity analysis (sampling approach)")
    plt.ylim(0, 1)
    plt.legend(frameon=False)
    plt.show()

# %%
# Notice how the sensitivity indices converge much slower this time. (Note the
# scale of the :math:`x`-axis!) We need a lot more samples to reach similar
# accuracy compared to the Polynomial Chaos approach.

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_sobol_g.py` for
#     an example of how to adaptively construct the multiindex set in the
#     Polynomial Chaos expansion to compute sensitivity indices.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_noisy_data.py` for
#     an example of how to use sparse estimators to deal with noisy measurements.
