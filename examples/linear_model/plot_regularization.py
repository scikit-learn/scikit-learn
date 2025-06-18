# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""
=================================================
Regularization in Linear and Logistic Regressions
=================================================

This example explores regularization techniques for linear and logistic regression
in both regression and classification tasks. It demonstrates how the
regularization parameter :math:`\\alpha` can be adjusted to control the complexity
of the trained coefficients :math:`w` and reduce overfitting.
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_classification, make_regression
from sklearn.linear_model import (
    LogisticRegression,
    Ridge,
    enet_path,
    lars_path,
    lasso_path,
)
from sklearn.metrics import mean_squared_error
from sklearn.svm import l1_min_c

# %%
#
# Regularization in Linear Regression
# -----------------------------------
#
# When performing linear regression on a given dataset
# :math:`(X, y)`, regularization terms can be added to
# control the model's complexity and mitigate overfitting.
# Scikit-learn provides the following regularization techniques:
#
# - :func:`~sklearn.linear_model.Lasso`
# - :func:`~sklearn.linear_model.LassoLars`
# - :func:`~sklearn.linear_model.Ridge`
# - :func:`~sklearn.linear_model.ElasticNet` (with `l1_ratio=0.5`
#   in this example)
#
# Mathematically, these are formulated by minimizing the constrained
# least-squares penalty:
#
# .. math::
#
#     \min_{w} \frac{1}{2n_{\operatorname{sample}}}
#     \Vert Xw - y \Vert^2_2 +
#     \left\{
#     \begin{array}{cl}
#     \alpha \Vert w \Vert_1 & \mbox{Lasso(-LARS)} \\
#     \alpha \Vert w \Vert_2^2 & \mbox{Ridge} \\
#     \frac{\alpha}{2} \Vert w \Vert_1  +
#     \frac{\alpha}{4} \Vert w \Vert^2_2& \mbox{Elastic Net} \\
#     \end{array}
#     \right.
#
# Thus, the Lasso model (resp. Ridge model) includes the :math:`\ell^1`-norm
# (resp. :math:`\ell^2`-norm) of the regression
# coefficients in the penalty, while the Elastic Net model
# incorporates both :math:`\ell^1`- and :math:`\ell^2`-norms.
#
# We can interpret the :math:`\ell^p`-norms as minimising the least squares penalty
# under different geometries. This is illustrated by plotting the unit circles
#
# .. math::
#
#     \left\{
#     \begin{array}{cl}
#     \Vert w \Vert_1 &=1 \\
#     \Vert w \Vert_2^2 &=1 \\
#     0.5 \Vert w \Vert_1 + 0.25 \Vert w \Vert_2^2 &= 1 \\
#     \end{array}
#     \right.
#
# in :math:`\mathbb{R}^2`:
line = np.linspace(-1.2, 1.2, 1001)
xx, yy = np.meshgrid(line, line)

l1 = np.abs(xx) + np.abs(yy)
l2 = xx**2 + yy**2
elastic_net = 0.5 * l1 + 0.25 * l2

plt.figure()
ax = plt.gca()

l1_contour = plt.contour(xx, yy, l1, levels=[1], colors="#0072B2")
l2_contour = plt.contour(xx, yy, l2, levels=[1], colors="#D55E00")
elastic_net_contour = plt.contour(xx, yy, elastic_net, levels=[1], colors="#009E73")

ax.set_aspect("equal")
ax.spines["left"].set_position("center")
ax.spines["right"].set_color("none")
ax.spines["bottom"].set_position("center")
ax.spines["top"].set_color("none")

plt.clabel(l1_contour, inline=1, fmt={1.0: r"$\ell^1$"}, manual=[(-1, -1)])
plt.clabel(l2_contour, inline=1, fmt={1.0: r"$\ell^2$"}, manual=[(-1, -1)])
plt.clabel(
    elastic_net_contour,
    inline=1,
    fmt={1.0: "Elastic Net"},
    manual=[(1, -2)],
)

plt.title(r"Unit Circles in $\mathbb{R}^2$")

plt.tight_layout()
_ = plt.show()

# %%
# Algebraically, any solution to this optimization problem depends
# on :math:`\alpha`. For example, in Lasso, a large :math:`\alpha` forces
# the least-squares penalty to stay small, which in turn keeps the norm
# :math:`\Vert w \Vert_1`
# small. Conversely, a smaller :math:`\alpha` allows the norm
# :math:`\Vert w \Vert_1`
# to grow larger.
#
# This suggests that the regression coefficients :math:`w` evolve as
# :math:`\alpha` increases, and we are interested in knowing
# :math:`w` across a range of :math:`\alpha` values. This is known
# as the **regularization path**: a list of :math:`w` values corresponding to
# different :math:`\alpha` values, ranging from small to large.
#
# In this example, we plot the regularization paths to show how the magnitudes of
# the coefficients change as the regularization parameter :math:`\alpha` increases.
# This demonstrates how model complexity varies with :math:`\alpha`. We then compare
# the trained coefficients with the true coefficients used to generate the training set,
# illustrating how regularization helps mitigate overfitting.
#
# Creating a Noise-free Regression Dataset
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# We generate a toy dataset with 400 samples and 10 features, suitable for
# regression analysis. Since the data is noise-free in this example,
# we can expect our regression model to recover the true coefficients `w` exactly.

X, y, w = make_regression(n_samples=400, n_features=10, coef=True, random_state=42)

# %%
#
# Impact of Regularization Parameter on Model Complexity
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Lasso(-LARS) and Elastic Net Models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Scikit-learn provides the following functions to compute multiple
# :math:`w` values for various :math:`\alpha` values efficiently:
#
# - :func:`~sklearn.linear_model.lasso_path`
# - :func:`~sklearn.linear_model.lars_path`
# - :func:`~sklearn.linear_model.enet_path` with `l1_ratio=0.5`
#

eps = 3e-4  # the smaller it is the longer is the path

alphas_lasso, coefs_lasso, _ = lasso_path(X, y, eps=eps)

alphas_enet, coefs_enet, _ = enet_path(X, y, eps=eps, l1_ratio=0.5)

# %%
# The :func:`~sklearn.linear_model.lasso_path` and
# :func:`~sklearn.linear_model.enet_path` functions compute
# :math:`w` with **coordinate decent**: for each entry of :math:`w`,
# the function solves for it optimal value while keeping the others
# fixed. Since the algorithm iterates until convergence,
# Lasso doesn't operate in a fixed number of steps based solely
# on the dataset's size, which can make it take longer to run.
# In contrast, the Lasso-LARS model computes the Lasso solution in fewer steps.

alphas_lars, _, coefs_lars = lars_path(X, y, method="lasso")

# %%
# The Lasso-LARS model uses the **Least Angle Regression (LARS)** algorithm
# (see [1]_ Algorithm 3.2 on page 74) to compute the Lasso solution in
# :math:`\min \left\{
# n_{\operatorname{sample}}-1,n_{\operatorname{feature}}
# \right\}`
# steps. This provides an efficient algorithm for computing the entire Lasso path, and
# is implemented as :func:`~sklearn.linear_model.LassoLars`
# and :func:`~sklearn.linear_model.lars_path`.
#
# Ridge Model
# ~~~~~~~~~~~
#
# Next, we compute the coefficients for the Ridge model using the :math:`\alpha`
# from Elastic Net:

coefs_ridge = []
for a in alphas_enet:
    ridge = Ridge(alpha=a)
    ridge.fit(X, y)
    coefs_ridge.append(ridge.coef_)

coefs_ridge = np.asarray(coefs_ridge)

# %%
# Plotting the Regularization Paths
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# We now present the visualization of the regularization paths for the dataset.
# Each model is represented by 10 curves, corresponding to the number of features in the
# dataset. Each curve shows how a particular coefficient :math:`w_i` changes as
# :math:`\alpha` increases.

model_names = ["Lasso", "Lasso-LARS", "Elastic Net", "Ridge"]

model_coefficients = [coefs_lasso.T, coefs_lars.T, coefs_enet.T, coefs_ridge]

model_alphas = [alphas_lasso, alphas_lars, alphas_enet, alphas_enet]

fig, axes = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(25, 10))

for i in range(len(model_names)):
    for j in range(len(model_names)):
        if i == j:
            axes[i, i].semilogx(model_alphas[i], model_coefficients[i], c="#0072B2")

            axes[i, i].set_title(f"{model_names[i]} Paths", fontsize=14)

        elif j < i:
            l1 = axes[i, j].semilogx(
                model_alphas[i], model_coefficients[i], c="#0072B2"
            )

            l2 = axes[i, j].semilogx(
                model_alphas[j], model_coefficients[j], linestyle="--", c="#D55E00"
            )

            axes[i, j].set_title(
                f"{model_names[j]} vs {model_names[i]} Paths", fontsize=14
            )

            axes[i, j].legend(
                (l1[-1], l2[-1]),
                (f"{model_names[i]}", f"{model_names[j]}"),
                loc="upper right",
            )

        else:
            fig.delaxes(axes[i, j])

fig.text(0.5, 0.02, r"$\alpha$", fontsize=18, ha="center")
fig.text(0, 0.5, "Coefficients", fontsize=18, va="center", rotation=90)

fig.suptitle(
    "Comparing Regularization Paths: Lasso(-LARS), Ridge, and Elastic Net", fontsize=20
)

fig.tight_layout(pad=3.0)
_ = plt.show()

# %%
#
# * In the "Lasso vs Lasso-LARS Paths" visual,
#   the Lasso and Lasso-LARS paths appear identical towards the end
#   because both models solve the same constrained problem.
#   However, Lasso-LARS reaches the solution faster than Lasso.
#
# * The "Lasso vs Elastic-Net Paths" visual is more notable.
#   Elastic Net's coefficients tend to have smaller absolute values than those of Lasso.
#   Additionally, Elastic Net maintains more non-zero coefficients than Lasso towards
#   the end.
#
# * In the "Lasso(-LARS) vs Ridge Paths" and "Elastic Net vs Ridge Paths" visuals, the
#   Ridge model focuses on shrinking all coefficients uniformly, rather than setting
#   some to exactly zero. As a result, the Ridge model retains all features after
#   training, unlike the Lasso(-LARS) or Elastic Net models.
#
# This demonstrates how different regularization techniques govern
# the model's complexity:
#
# 1. the :math:`\ell^1`-norm constraint encourages sparsity in the solution.
#
# 2. the :math:`\ell^2`-norm constraint focuses on shrinkage of the magnitude
#    of the solution.
#
# 3. the Elastic Net constraint provides a balanced compromise.
#
# Mitigating Overfitting with Regularization
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Recall that the true coefficient `w` refers to the coefficients of the linear model
# used to generate the training dataset. In this section, we compare the trained
# coefficients of Lasso(-LARS), Ridge, and Elastic Net with `w` to demonstrate how
# regularization can mitigate overfitting. This is achieved by computing the
# :func:`~sklearn.metrics.mean_squared_error` (MSE) between the true and trained
# coefficients.

lasso_mse, lars_mse, enet_mse, ridge_mse = [], [], [], []

for coef_lasso, coef_enet, coef_ridge in zip(coefs_lasso.T, coefs_enet.T, coefs_ridge):
    lasso_mse.append(mean_squared_error(coef_lasso, w))
    enet_mse.append(mean_squared_error(coef_enet, w))
    ridge_mse.append(mean_squared_error(coef_ridge, w))

for coef_lars in coefs_lars.T:
    lars_mse.append(mean_squared_error(coef_lars, w))

lasso_mse = np.asarray(lasso_mse)
lars_mse = np.asarray(lars_mse)
enet_mse = np.asarray(enet_mse)
ridge_mse = np.asarray(ridge_mse)

# %%
#
# The idea is that a smaller MSE between the true and trained coefficients implies
# greater similarity between the coefficients. Thus, if the MSE is small, the
# trained model captures the underlying pattern of the training data well.
# However, this can also indicate that the trained model may not perform well on
# generalised data, as the pattern may not hold for unseen data.
# This is essentially the overfitting problem.
#
# The following visualization demonstrates how the MSE changes for different trained
# models as the regularization parameter :math:`\alpha` increases.

plt.figure()
l1 = plt.semilogx(alphas_lasso, lasso_mse.T, c="#0072B2")
l2 = plt.semilogx(alphas_lars, lars_mse.T, c="#D55E00")
l3 = plt.semilogx(alphas_enet, enet_mse.T, c="#009E73")
l4 = plt.semilogx(alphas_enet, ridge_mse, c="#F0E442")

plt.xlabel(r"$\alpha$")
plt.ylabel("Mean Squared Error")
plt.title("Coefficient Error Across Regularization Strengths")
plt.legend(
    (l1[-1], l2[-1], l3[-1], l4[-1]),
    ("Lasso", "LARS", "Elastic Net", "Ridge"),
    loc="upper left",
)

plt.axis("tight")
_ = plt.show()

# %%
#
# In the visualization, for small values of :math:`\alpha`, since our synthetic data is
# noise-free, the trained coefficients of Lasso(-LARS), Ridge, and Elastic Net are
# similar to the true coefficient `w` (with MSE close to 0). This indicates that the
# models capture the intricate details of the training data well.
#
# As :math:`\alpha` increases, the MSE also increases. This improves the models' ability
# to generalise to unseen data (e.g., if the data were noisy), but it also risks
# degrading model performance if the regularization becomes too strong.
#
# Regularization in Logistic Regression
# -------------------------------------
#
# Regularization can also be applied to Logistic Regression when working on
# classification tasks. scikit-learn's :func:`~sklearn.linear_model.LogisticRegression`
# enables users to apply regularization using the `penalty` parameter:
#
# * `l1`: :math:`\ell^1`-regularization, similar to the Lasso model
# * `l2`: :math:`\ell^2`-regularization, similar to the Ridge model
# * `elasticnet`: Combined with the `l1_ratio` parameter for a mix of :math:`\ell^1`
#   and :math:`\ell^2`
#
# Additionally, the `C` parameter controls the inverse of the regularization strength.
# Smaller values of `C` apply stronger regularization.
#
# We demonstrate the effect of regularization by creating a synthetic classification
# dataset.
#

X, y = make_classification(
    n_samples=400,
    n_features=64,
    n_informative=64,
    n_redundant=0,
    n_classes=2,
    n_clusters_per_class=1,
    random_state=42,
)

# %%
#
# In this synthetic binary classification dataset, there are 400 samples,
# each with 64 features. This toy dataset is noise-free to maintain consistency with
# our earlier regression example.
#
# As noted in the regression example, :math:`\ell^1`-regularization may set some
# coefficients exactly to zero. For extreme values of `C`, the trained coefficients
# may even become the zero vector. To address this, scikit-learn provides the
# :func:`~sklearn.svm.l1_min_c` function, which computes the minimum value of the
# regularization strength `C` at which the model begins to learn meaningful patterns
# (i.e., some coefficients become non-zero).
#

cs = l1_min_c(X, y, loss="log") * np.logspace(0, 10, 16)

# %%
#
# We now plot heatmaps to represent the sparsity for each `penalty` and each value
# of `C`.
#

l1_ratio = 0.5  # l1 weight in the Elastic-Net regularization

fig, axes = plt.subplots(3, 3)

# Set regularization parameter
for i, (C, axes_row) in enumerate(zip((1, 0.1, 0.01), axes)):
    # Increase tolerance for short training time
    clf_l1_LR = LogisticRegression(C=C, penalty="l1", tol=0.01, solver="saga")
    clf_l2_LR = LogisticRegression(C=C, penalty="l2", tol=0.01, solver="saga")
    clf_en_LR = LogisticRegression(
        C=C, penalty="elasticnet", solver="saga", l1_ratio=l1_ratio, tol=0.01
    )
    clf_l1_LR.fit(X, y)
    clf_l2_LR.fit(X, y)
    clf_en_LR.fit(X, y)

    coef_l1_LR = clf_l1_LR.coef_.ravel()
    coef_l2_LR = clf_l2_LR.coef_.ravel()
    coef_en_LR = clf_en_LR.coef_.ravel()

    sparsity_l1_LR = np.mean(coef_l1_LR == 0) * 100
    sparsity_l2_LR = np.mean(coef_l2_LR == 0) * 100
    sparsity_en_LR = np.mean(coef_en_LR == 0) * 100

    if i == 0:
        axes_row[0].set_title(r"$\ell^1$ penalty")
        axes_row[1].set_title(f"Elastic-Net\n {l1_ratio = }")
        axes_row[2].set_title(r"$\ell^2$ penalty")

    for ax, coefs in zip(axes_row, [coef_l1_LR, coef_en_LR, coef_l2_LR]):
        ax.imshow(
            np.abs(coefs.reshape(8, 8)),
            interpolation="nearest",
            cmap="binary",
            vmax=1,
            vmin=0,
        )
        ax.set_xticks(())
        ax.set_yticks(())

    axes_row[0].set_ylabel(f"C = {C:.2f}")

_ = plt.show()
# %%
#
# Each heatmap organizes the 64 coefficients (the number of features in our synthetic
# classification dataset) into an 8×8 grid. It is constructed by taking the absolute
# values of the coefficients and displaying them in a black-and-white scale, where
# lower values appear white and higher values appear black.
#
# We can see that larger values of `C` (i.e., weaker regularization) give the model
# more freedom, while smaller values of `C` impose stronger constraints, leading to
# increased sparsity. As expected, the Elastic-Net penalty results in a level of
# sparsity between that of :math:`\ell^1` and :math:`\ell^2`.
#
# .. rubric:: References
#
# .. [1] Hastie, T., Tibshirani, R., & Friedman, J. (2009). The Elements of Statistical
#        Learning: Data Mining, Inference, and Prediction. New York,
#        NY: Springer New York.
