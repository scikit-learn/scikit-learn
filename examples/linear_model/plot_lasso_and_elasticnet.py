"""
========================================
Lasso and Elastic Net for Sparse Signals
========================================

The present example compares three l1-based regression models on a manually
generated sparse signal corrupted with an additive gaussian noise:

 - a :ref:`lasso`;
 - an :ref:`automatic_relevance_determination`;
 - an :ref:`elastic_net`.

We compute the :math:`R^2` score and the fitting time of the models. The
estimated coefficients of each model are compared with the ground-truth.
"""

# %%
# Generate synthetic dataset
# --------------------------
#
# We generate a dataset where the number of samples is lower than the total
# number of features. This leads to a singular matrix that cannot be dealt with
# an :ref:`ordinary_least_squares` as it requires some regularization.
#
# Here `X` and `y` are linearly linked. Only 100 out of the 1_000 features in
# `X` will be used to generate `y` while the rest are not useful at predicting
# `y`, resulting in a sparse feature space. Some correlations between the
# informative features are added by means of the `effective_rank` parameter.
# Finally, gaussian noise is added.
#
# Such sparse, noisy and correlated features can be obtained, for instance, from
# sensor nodes monitoring some environmental variables, as they typically register
# similar values depending on their positions (spatial correlations).
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

X, y, true_weights = make_regression(
    n_samples=500,
    n_features=1000,
    n_informative=100,
    effective_rank=15,
    noise=1,
    coef=True,
    random_state=42,
)

# %%
# We split data in train set and test set for simplicity. In practice one should
# use cross-validation to estimate the variance of the test score.
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.5, random_state=42
)

# %%
# Lasso
# -----
#
# In this example we demo a :class:`~sklearn.linear_model.Lasso` with a fix
# value of the regularization parameter `alpha`. To select the optimal value for
# `alpha` we used a :class:`~sklearn.linear_model.LassoCV`. This is not shown
# here to keep the example simple.
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score
from time import time

common_params = {"alpha": 0.00002, "max_iter": 10_000}
t0 = time()
lasso = Lasso(**common_params).fit(X_train, y_train)
print(f"Lasso fit done in {(time() - t0):.3f}s")

y_pred_lasso = lasso.predict(X_test)
r2_score_lasso = r2_score(y_test, y_pred_lasso)
print(f"Lasso r^2 on test data : {r2_score_lasso:.3f}")

# %%
# ARD
# ---
#
# An ARD regression is the bayesian version of the Lasso. It can produce
# interval estimates for all of the parameters, including the error variance, if
# required.

from sklearn.linear_model import ARDRegression

t0 = time()
ard = ARDRegression().fit(X_train, y_train)
print(f"ARD fit done in {(time() - t0):.3f}s")

y_pred_ard = ard.predict(X_test)
r2_score_ard = r2_score(y_test, y_pred_ard)
print(f"ARD r^2 on test data : {r2_score_ard:.3f}")

# %%
# ElasticNet
# ----------
#
# Additionally to the l1-penalty, the class
# :class:`~sklearn.linear_model.ElasticNet` introduces a l2-penalty by means of
# the parameter `l1_ratio`.  For `l1_ratio = 0` the penalty is l2 and the model
# is equivalent to a :class:`~sklearn.linear_model.RidgeRegression`. For
# `l1_ratio = 1` it is an l1 penalty and the model is equivalent to a
# :class:`~sklearn.linear_model.Lasso`.  For `0 < l1_ratio < 1`, the penalty is
# a combination of l1 and l2.
#
# Similarly to the Lasso model, here we train the model with fix values for
# `alpha` and `l1_ratio`. To select their optimal value we used an
# :class:`~sklearn.linear_model.ElasticNetCV`, not shown here to keep the
# example simple.

from sklearn.linear_model import ElasticNet

t0 = time()
enet = ElasticNet(l1_ratio=0.8, **common_params).fit(X_train, y_train)
print(f"ElasticNet fit done in {(time() - t0):.3f}s")

y_pred_enet = enet.predict(X_test)
r2_score_enet = r2_score(y_test, y_pred_enet)
print(f"ElasticNet r^2 on test data : {r2_score_enet:.3f}")

# %%
# Plot results
# ------------
#
# In this section we use a heatmap to visualize the sparsity of the true
# and estimated coefficients of the respective linear models.

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import SymLogNorm

df = pd.DataFrame(
    {
        "True weights": true_weights,
        "Lasso": lasso.coef_,
        "ARDRegression": ard.coef_,
        "ElasticNet": enet.coef_,
    }
)

plt.figure(figsize=(10, 6))
ax = sns.heatmap(
    df.T,
    norm=SymLogNorm(linthresh=10e-4, vmin=-80, vmax=80),
    cbar_kws={"label": "coefficients' values"},
    cmap="seismic_r",
)
plt.ylabel("linear model")
plt.xlabel("coefficients")
plt.title(
    f"Models' coefficients\nLasso $R^2$: {r2_score_lasso:.3f}, "
    f"ARD $R^2$: {r2_score_ard:.3f}, "
    f"ElasticNet $R^2$: {r2_score_enet:.3f}"
)
plt.tight_layout()

# %%
# Conclusions
# -----------
#
# :class:`~sklearn.linear_model.Lasso` is known to recover sparse data very well
# but does not perform well with highly correlated features. Indeed, if several
# correlated features contribute to the target, :class:`~sklearn.linear_model.Lasso`
# will end-up selecting a single one.
#
# :class:`~sklearn.linear_model.ARDRegression` is better when handling gaussian
# noise, which translates in a slightly better score than Lasso, but is still
# unable to handle correlated features and requires a large amount of time due
# to fitting a prior.
#
# :class:`~sklearn.linear_model.ElasticNet` does not have a problem with
# correlated features but is unable to recover the true generating process of
# sparse data as good as Lasso. In this case, this is the model with the best
# score and lowest fitting time. In the absence of correlations, an
# :class:`~sklearn.linear_model.ElasticNetCV` would find an optimal `l1_ratio =
# 1`, i.e. the optimal model would be a :class:`~sklearn.linear_model.Lasso`.
