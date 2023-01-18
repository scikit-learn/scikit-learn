"""
========================================
Lasso and Elastic Net for Sparse Signals
========================================

Estimates Lasso and Elastic-Net regression models on a manually generated
sparse signal corrupted with an additive noise. Estimated coefficients are
compared with the ground-truth.

"""

# %%
# Data Generation
# ---------------------------------------------------

from sklearn.datasets import make_regression

X, y, true_weights = make_regression(
    n_samples=500,
    n_features=1000,
    n_informative=100,
    effective_rank=15,
    noise=1,
    coef=True,
    random_state=42,
)

# Split data in train set and test set
n_samples = X.shape[0]
X_train, y_train = X[: n_samples // 2], y[: n_samples // 2]
X_test, y_test = X[n_samples // 2 :], y[n_samples // 2 :]

# %%
# Lasso
# ---------------------------------------------------

from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score

common_params = {
    "alpha": 0.00002,
    "max_iter": 10_000,
}
lasso = Lasso(**common_params).fit(X_train, y_train)

y_pred_lasso = lasso.predict(X_test)
r2_score_lasso = r2_score(y_test, y_pred_lasso)

print("r^2 on test data : %f" % r2_score_lasso)

# %%
# ElasticNet
# ---------------------------------------------------

from sklearn.linear_model import ElasticNet

enet = ElasticNet(l1_ratio=0.8, **common_params).fit(X_train, y_train)

y_pred_enet = enet.predict(X_test)
r2_score_enet = r2_score(y_test, y_pred_enet)

print("r^2 on test data : %f" % r2_score_enet)


# %%
# Plot
# ---------------------------------------------------

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import SymLogNorm

df = pd.DataFrame(
    {
        "True weights": true_weights,
        "Lasso": lasso.coef_,
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
    f"ElasticNet $R^2$: {r2_score_enet:.3f}"
)
plt.tight_layout()
