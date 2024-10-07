"""
========================================
Lasso, Lasso-LARS, and Elastic Net paths
========================================

This example shows how to compute the "paths" of coefficients along the Lasso,
Lasso-LARS, and Elastic Net regularization paths. It illustrates the
relationship between the regularization parameter :math:`\\alpha`
and the coefficients :math:`w`.

When performing linear regression on a given dataset
:math:`(X, y)`, regularization terms can be added to
control the model's complexity.
Scikit-learn provides the following regularization techniques:

- :func:`~sklearn.linear_model.Lasso`
- :func:`~sklearn.linear_model.LassoLars`
- :func:`~sklearn.linear_model.ElasticNet` with default `l1_ratio=0.5`

Mathematically, these are formulated by minimising the constrained
least-squares penalty:

.. math::

    \\min_{w} \\frac{1}{2n_{\\operatorname{sample}}}
    \\vert \\vert Xw - y \\vert \\vert^2_2 +
    \\left\\{
    \\begin{array}{cl}
    \\alpha \\vert \\vert w \\vert \\vert_1 & \\mbox{Lasso(-LARS)} \\\\
    \\frac{\\alpha}{2} \\vert \\vert w \\vert \\vert_1  +
    \\frac{\\alpha}{4} \\vert \\vert w \\vert \\vert^2_2& \\mbox{Elastic Net} \\\\
    \\end{array}
    \\right.

Thus, the Lasso model includes the :math:`\\ell^1`-norm of the regression
coefficients in the penalty, while the Elastic Net model
incorporates both :math:`\\ell^1`- and :math:`\\ell^2`-norms.

Any solution to this optimisation problem depends on :math:`\\alpha`.
For example, in Lasso, a large :math:`\\alpha` forces the least-squares
penalty to stay small, which in turn keeps the norm
:math:`\\vert \\vert w \\vert \\vert_1`
small. Conversely, a smaller :math:`\\alpha` allows the norm
:math:`\\vert \\vert w \\vert \\vert_1`
to grow larger.

This suggests that the regression coefficients :math:`w` evolve as
:math:`\\alpha` increases, and we are interested in knowing
:math:`w` across a range of :math:`\\alpha` values. This is known
as the **regularization path**: a list of :math:`w` values corresponding to
different :math:`\\alpha` values, ranging from small to large.
In this example, we plot the regularization paths to show how the sizes of the
coefficients change as the regularization parameter increases.

We use the diabetes dataset :func:`~sklearn.datasets.load_diabetes` to plot
the regression coefficients for Lasso and Elastic Net. Scikit-learn provides the
following functions to compute multiple :math:`w` values for various
:math:`\\alpha` values efficiently:

- :func:`~sklearn.linear_model.lasso_path`
- :func:`~sklearn.linear_model.lars_path`
- :func:`~sklearn.linear_model.enet_path`

The :func:`~sklearn.linear_model.lasso_path` and
:func:`~sklearn.linear_model.enet_path` functions compute
:math:`w` with coordinate decent: for each entry of :math:`w`,
the function solves for it optimal value while keeping the others
fixed. Since the algorithm iterates until convergence,
Lasso doesn't operate in a fixed number of steps based solely
on the dataset's size, which can make it take longer to run.
In contrast, the Lasso-LARS model computes the Lasso solution in fewer steps.

The Lasso-LARS model uses the Least Angle Regression (LARS) algorithm
(see [1]_) to compute the Lasso solution in
:math:`\\min \\left\\{
n_{\\operatorname{sample}}-1,n_{\\operatorname{feature}}
\\right\\}`
steps. This provides an efficient algorithm for computing the entire Lasso path, and
is implemented as :func:`~sklearn.linear_model.LassoLars`
and :func:`~sklearn.linear_model.lars_path`.

We now present the visualisation of the regularization path for the diabetes dataset.
Each model is represented by 10 curves, corresponding to the number of features in the
dataset. Each curve shows how a particular coefficient :math:`w_i` changes as
:math:`\\alpha` increases.

- In the "Lasso vs LARS Paths" visual, the Lasso and LARS paths appear identical because
  both models solve the same constrained problem. However, LARS reaches the solution
  faster than Lasso.

- The "Lasso vs Elastic-Net Paths" visual is more notable. Elastic Net's coefficients
  tend to have smaller absolute values than those of Lasso. Additionally, Elastic Net
  maintains more non-zero coefficients than Lasso towards the end. This demonstrates
  how the :math:`\\ell^1`-norm constraint encourages sparsity in the solution, while
  combining it with the :math:`\\ell^2`-norm provides a balanced compromise.



.. rubric:: References

.. [1] "The Elements of Statistical Learning", Hastie T., Tibshirani R.,
       Friedman J., Algorithm 3.2, p. 74, 2008.


"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt

from sklearn.datasets import load_diabetes
from sklearn.linear_model import enet_path, lars_path, lasso_path

X, y = load_diabetes(return_X_y=True)
X /= X.std(axis=0)  # Standardize data (easier to set the l1_ratio parameter)

# Compute paths

eps = 5e-3  # the smaller it is the longer is the path

print("Computing regularization path using the lasso...")
alphas_lasso, coefs_lasso, _ = lasso_path(X, y, eps=eps)

print("Computing regularization path using the positive lasso...")
alphas_positive_lasso, coefs_positive_lasso, _ = lasso_path(
    X, y, eps=eps, positive=True
)

print("Computing regularization path using the LARS...")
alphas_lars, _, coefs_lars = lars_path(X, y, method="lasso")

print("Computing regularization path using the positive LARS...")
alphas_positive_lars, _, coefs_positive_lars = lars_path(
    X, y, method="lasso", positive=True
)

print("Computing regularization path using the elastic net...")
alphas_enet, coefs_enet, _ = enet_path(X, y, eps=eps, l1_ratio=0.8)

print("Computing regularization path using the positive elastic net...")
alphas_positive_enet, coefs_positive_enet, _ = enet_path(
    X, y, eps=eps, l1_ratio=0.8, positive=True
)

# Display results

plt.figure(1)
for coef_lasso, coef_lars in zip(coefs_lasso, coefs_lars):
    l1 = plt.semilogx(alphas_lasso, coef_lasso, c="#0072B2")
    l2 = plt.semilogx(alphas_lars, coef_lars, linestyle="--", c="#D55E00")

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("Lasso vs LARS Paths")
plt.legend((l1[-1], l2[-1]), ("Lasso", "LARS"), loc="lower right")
plt.axis("tight")

plt.figure(2)
for coef_l, coef_e in zip(coefs_lasso, coefs_enet):
    l1 = plt.semilogx(alphas_lasso, coef_l, c="#0072B2")
    l2 = plt.semilogx(alphas_enet, coef_e, linestyle="--", c="#D55E00")

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("Lasso vs Elastic-Net Paths")
plt.legend((l1[-1], l2[-1]), ("Lasso", "Elastic-Net"), loc="lower right")
plt.axis("tight")

plt.figure(3)
for coef_l, coef_pl in zip(coefs_lasso, coefs_positive_lasso):
    l1 = plt.semilogy(alphas_lasso, coef_l, c="#0072B2")
    l2 = plt.semilogy(alphas_positive_lasso, coef_pl, linestyle="--", c="#D55E00")

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("Lasso vs Positive Lasso Paths")
plt.legend((l1[-1], l2[-1]), ("Lasso", "positive Lasso"), loc="lower right")
plt.axis("tight")

plt.figure(4)
for coef_lars, coef_positive_lars in zip(coefs_lars, coefs_positive_lars):
    l1 = plt.semilogx(alphas_lars, coef_lars, c="#0072B2")
    l2 = plt.semilogx(
        alphas_positive_lars, coef_positive_lars, linestyle="--", c="#D55E00"
    )

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("LARS vs Positive LARS Paths")
plt.legend((l1[-1], l2[-1]), ("LARS", "Positive LARS"), loc="lower right")
plt.axis("tight")

plt.figure(5)
for coef_e, coef_pe in zip(coefs_enet, coefs_positive_enet):
    l1 = plt.semilogx(alphas_enet, coef_e, c="#0072B2")
    l2 = plt.semilogx(alphas_positive_enet, coef_pe, linestyle="--", c="#D55E00")

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("Elastic-Net vs Positive Elastic-Net Paths")
plt.legend((l1[-1], l2[-1]), ("Elastic-Net", "positive Elastic-Net"), loc="lower right")
plt.axis("tight")
plt.show()
