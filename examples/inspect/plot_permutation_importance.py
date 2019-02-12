"""
==================================================
Permutation Importance vs Random Forest Importance
==================================================

The random forest `feature_importances_`, are computed from train set
statistics and are subject to bias with the cardinality of the feature. The
permutation importance of a feature is calculated by measuring how much the
model performance decreases when the feature is permutated.

In this example, we add a column of random numbers to the diabetes dataset.
Then we fit a :class:`sklearn.ensemble.RandomForestRegressor` to this modified
dataset. The feature importance from the random forest is plotted. In this
case, the ``RANDOM`` feature is considerd more important than the ``age`` or
``sex`` feature.

Next, we use :func:`sklearn.inspect.permutation_importance` to calcuate the
permutation importance for each feature.
The `sklearn.inspect.permutation_importance` returns a numpy array where
values in each row are the cross-validated scores for a feature. The
permutation importance for the random forest is plotted. In this case,
The ``RANDOM`` feature is less important than ``sex`` and ``age``.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_diabetes
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspect import permutation_importance


def plot_importances(importances, features, highlight=None, ax=None):
    N = features.shape[0]

    if ax is None:
        _, ax = plt.subplots()
    y_ticks = range(1, N + 1)
    arg_sorted = np.argsort(importances)

    color = ["blue" for _ in range(N)]
    labels = features[arg_sorted]

    if highlight is not None:
        for idx, label in enumerate(labels):
            if label == highlight:
                color[idx] = "red"

    ax.barh(y_ticks, importances[arg_sorted], color=color)
    ax.set_yticks(y_ticks)
    ax.set_xlim(0, np.max(importances)*1.05)
    ax.set_ylim(0, N + 1)
    ax.set_yticklabels(features[arg_sorted])


ds = load_diabetes()
X, y = ds.data, ds.target
features = np.array(ds.feature_names + ["RAND"])
rng = np.random.RandomState(42)
X = np.hstack([X, rng.normal(scale=1, size=(X.shape[0], 1))])

rf = RandomForestRegressor(n_estimators=50, random_state=rng)
rf.fit(X, y)

fig, (ax1, ax2) = plt.subplots(1, 2)
plot_importances(rf.feature_importances_, features, highlight="RAND", ax=ax1)
ax1.set_title("Feature importance from random forest")

perm_importances = permutation_importance(rf, X, y, random_state=rng,
                                          scoring="explained_variance")
perm_importances_mean = perm_importances.mean(axis=1)
plot_importances(perm_importances_mean, features, highlight="RAND", ax=ax2)
ax2.set_title("Permutation importance")
fig.tight_layout()
plt.show()
