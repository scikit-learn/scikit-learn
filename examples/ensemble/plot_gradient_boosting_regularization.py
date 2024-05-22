"""
================================
Gradient Boosting regularization
================================

Illustration of the effect of different regularization strategies
for Gradient Boosting. The example is taken from Chapter 10.12 of
Hastie et al 2009 [1]_.

The loss function used is log loss, aka binomial deviance. Regularization via
shrinkage (``learning_rate < 1.0``) improves performance considerably.
In combination with shrinkage, stochastic gradient boosting
(``subsample < 1.0``) can produce more accurate models by reducing the
variance via bagging.
Subsampling without shrinkage usually does poorly.
Another strategy to reduce the variance is by subsampling the features
analogous to the random splits in Random Forests
(via the ``max_features`` parameter).

.. [1] T. Hastie, R. Tibshirani and J. Friedman, "Elements of Statistical
    Learning Ed. 2", Springer, 2009.

"""

# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn import datasets
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split

X, y = datasets.make_hastie_10_2(n_samples=4000, random_state=1)

# map labels from {-1, 1} to {0, 1}
labels, y = np.unique(y, return_inverse=True)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.8, random_state=0)

original_params = {
    "n_estimators": 400,
    "max_leaf_nodes": 4,
    "max_depth": None,
    "random_state": 2,
    "min_samples_leaf": 2,
}

fig, axes = plt.subplots(ncols=2, figsize=(10, 5), sharex=True, sharey=True)

for j, model_class in enumerate(
    [GradientBoostingClassifier, HistGradientBoostingClassifier]
):
    for label, color, setting in [
        ("No shrinkage", "orange", {"learning_rate": 1.0, "subsample": 1.0}),
        ("learning_rate=0.2", "turquoise", {"learning_rate": 0.2, "subsample": 1.0}),
        ("subsample=0.5", "blue", {"learning_rate": 1.0, "subsample": 0.5}),
        (
            "learning_rate=0.2, subsample=0.5",
            "gray",
            {"learning_rate": 0.2, "subsample": 0.5},
        ),
        (
            "learning_rate=0.2, max_features=2",
            "magenta",
            {"learning_rate": 0.2, "max_features": 2},
        ),
    ]:
        params = dict(original_params)
        params.update(setting)
        n_iter = params["n_estimators"]
        if model_class == HistGradientBoostingClassifier:
            params["max_iter"] = params.pop("n_estimators")
            if "max_features" in params:
                params["max_features"] = float(
                    params["max_features"] / X_train.shape[1]
                )

        clf = model_class(**params)
        clf.fit(X_train, y_train)

        # compute test set deviance
        test_loss = np.zeros((n_iter,), dtype=np.float64)

        for i, y_proba in enumerate(clf.staged_predict_proba(X_test)):
            test_loss[i] = 2 * log_loss(y_test, y_proba[:, 1])

        axes[j].plot(
            (np.arange(test_loss.shape[0]) + 1)[::5],
            test_loss[::5],
            "-",
            color=color,
            label=label,
        )

    axes[j].set_ylim(None, 2)
    axes[j].legend(loc="upper right")
    axes[j].set_xlabel("Boosting Iterations")
    axes[j].set_ylabel("Test Set Log Loss")
    axes[j].set_title(model_class.__name__)
