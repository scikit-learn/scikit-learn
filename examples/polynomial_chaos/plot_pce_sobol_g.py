r"""
===============================================================================
Global sensitivity analysis of the Sobol function
===============================================================================

This example illustrates how to use Polynomial Chaos regression with adaptive
basis growth to perform a Global Sensitivity Analysis of the Sobol function.
The Sobol function (also known as the :math:`G` function, because it was
denoted as such in the original paper by Sobol), is a well-known test problem
in uncertainty quantification and sensitivity analysis [Saltelli2000]_. It is
defined as

.. math::
    y = \prod_{j=1}^d \frac{|4x_j - 2| + a_j}{1 + a_j}

where the inputs :math:`x_j, j = 1, 2, \ldots, d` follow a standard uniform
distribution, and the coefficients :math:`a_j` are nonnegative.

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
# First, let's define the parameters in the model.
a = np.array([1, 2, 5, 20])
dimension = len(a)

# %%
# The inputs :math:`x_j, j=1, 2, \ldots, d` follow a standard uniform
# distribution. Let's generate some training data by sampling from the input
# distributions.
from scipy.stats import uniform

distribution = uniform()
X = distribution.rvs((85, dimension))

# %%
# Next, evaluate the model in each of the training samples.
from math import prod

y = prod((abs(4 * X_j - 2) + a_j) / (1 + a_j) for a_j, X_j in zip(a, X.T))

# %%
# Then, define a Polynomial Chaos expansion and fit the data. We use a 2nd
# order total degree polynomial basis (the default value). The number of basis
# terms :math:`P` can be computed as
#
# .. math::
#   P = \frac{(d + k)!}{d!k!}
#
# where :math:`d` is the dimension (i.e., the number of input features), and
# :math:`k` is the degree of the polynomial. In our case, we have :math:`d = 4`
# and :math:`k = 2`, so there are :math:`P = 15` basis terms.
from sklearn.polynomial_chaos import PolynomialChaosExpansion

pce = PolynomialChaosExpansion(distribution)
pce.fit(X, y)

# %%
# Ater the expansion is fitted, we can extract the main-effect sensitivity
# indices from the Polynomial Chaos expansion. The :math:`j`\ th main-effect
# index indicates the effect on the model output of varying the :math:`j`\ th
# parameter alone. In this example, the variance of the output and the Sobol
# sensitivity indices can be computed analytically as
#
# .. math::
#    S_j = \frac{D_j}{\mathbb{V}[y]}
#
# for :math:`j = 1, 2, \ldots, d`, where
#
# .. math::
#    \mathbb{V}[y] = \prod_{j=1}^d (D_j + 1) - 1
#
# with
#
# .. math::
#    D_j = \frac{1}{3(1 + a_j)^2}.
#
import matplotlib.pyplot as plt
from pandas import DataFrame

D = 1 / (3 * (1 + a) ** 2)
V = prod(D + 1) - 1
S = D / V
df = DataFrame(
    {"predicted": pce.main_sens().flatten(), "exact": S},
    index=[f"S{j}" for j in range(dimension)],
)
df["predicted"].plot.bar()
plt.ylabel("main sensitivity indices")
df

# %%
# Even with a low-order Polynomial Chaos expansion, we obtain a reasonably good
# approximation for the sensitivity indices. We can improve the approximation
# by using a higher-order PC expansion. However, even in this
# moderate-dimensional case (dimension 4), the number of basis terms grows
# rapidly with the maximum degree of the polynomials.
from math import factorial

N = [
    factorial(dimension + degree) // factorial(dimension) // factorial(degree)
    for degree in range(10)
]
DataFrame(N, columns=["# basis terms"]).rename_axis(index="degree")

# %%
# For an order-4 polynomial basis, the number of basis terms (and hence the
# number of unknown coefficients) would be larger than the number of data
# points. As a consequence, the resulting polynomial system would become
# underdetermined, and regularization is required to obtain a solution.
# We may also try to change from `total_degree` multiindex sets to a more
# sparse multiindex set shape, such as a `hyperbolic_cross`. However, we would
# inevitably run into the same problem.

# %%
# An other approach to avoid this problem altogether is to use an adaptive
# construction of the multiindex set. In the code below, we track how the
# relative error in the main sensitivity indices changes as a function of the
# number of iterations in the adaptive algorithm. We also use an :math:`L1`
# prior for regularization, as implemented in the Lasso method.
from sklearn.linear_model import LassoCV

pce = PolynomialChaosExpansion(distribution, degree=1)
pce.fit(X, y)
errors = [np.linalg.norm(pce.main_sens() - S) / np.linalg.norm(S)]
iters = [1]
estimator = LassoCV(
    fit_intercept=False, alphas=np.logspace(-12, 2, 25), max_iter=100000
)
for i in range(10, 35, 5):
    pce = PolynomialChaosExpansion(distribution, degree=1, estimator=estimator)
    pce.fit(X, y, n_iter=i)
    errors.append(np.linalg.norm(pce.main_sens() - S) / np.linalg.norm(S))
    iters.append(i)

# %%
# Finally, we plot the decay of the relative error as a function of the number
# of iterations in the adaptive algorithm, indicating an initial rapid decay in
# the relative error.
plt.plot(iters, errors)
plt.yscale("log")
plt.xlabel("number of iterations")
plt.ylabel("relative error in main sensitivity indices")
_ = plt.title("Convergence of adaptive basis construction")

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_ishigami.py` for
#     a global sensitivity analysis of the Ishigami function, another
#     well-known test problem.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_noisy_data.py` for
#     an example of how to use sparse estimators to deal with noisy measurements.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_feature_selection_g.py`
#     for an example of how to use pruning to remove small basis terms from
#     the expansion.
