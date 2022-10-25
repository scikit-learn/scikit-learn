"""
=====================
Monotonic Constraints
=====================

This example illustrates the effect of monotonic constraints on a gradient
boosting estimator.

We build an artificial dataset where the target value is in general
positively correlated with the first feature (with some random and
non-random variations), and in general negatively correlated with the second
feature.

By imposing a monotonic increase or a monotonic decrease constraint, respectively,
on the features during the learning process, the estimator is able to properly follow
the general trend instead of being subject to the variations.

This example was inspired by the `XGBoost documentation
<https://xgboost.readthedocs.io/en/latest/tutorials/monotonic.html>`_.

"""

from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import PartialDependenceDisplay
import numpy as np
import matplotlib.pyplot as plt


rng = np.random.RandomState(0)

n_samples = 5000
f_0 = rng.rand(n_samples)
f_1 = rng.rand(n_samples)
X = np.c_[f_0, f_1]
noise = rng.normal(loc=0.0, scale=0.01, size=n_samples)

# y is positively correlated with f_0, and negatively correlated with f_1
y = 5 * f_0 + np.sin(10 * np.pi * f_0) - 5 * f_1 - np.cos(10 * np.pi * f_1) + noise

fig, ax = plt.subplots()


# Without any constraint
gbdt = HistGradientBoostingRegressor()
gbdt.fit(X, y)
disp = PartialDependenceDisplay.from_estimator(
    gbdt,
    X,
    features=[0, 1],
    feature_names=(
        "First feature",
        "Second feature",
    ),
    line_kw={"linewidth": 4, "label": "unconstrained", "color": "tab:blue"},
    ax=ax,
)

# With monotonic increase (1) and a monotonic decrease (-1) constraints, respectively.
gbdt = HistGradientBoostingRegressor(monotonic_cst=[1, -1])
gbdt.fit(X, y)

PartialDependenceDisplay.from_estimator(
    gbdt,
    X,
    features=[0, 1],
    line_kw={"linewidth": 4, "label": "constrained", "color": "tab:orange"},
    ax=disp.axes_,
)

for f_idx in (0, 1):
    disp.axes_[0, f_idx].plot(
        X[:, f_idx], y, "o", alpha=0.3, zorder=-1, color="tab:green"
    )
    disp.axes_[0, f_idx].set_ylim(-6, 6)

plt.legend()
fig.suptitle("Monotonic constraints effect on partial dependences")

plt.show()
