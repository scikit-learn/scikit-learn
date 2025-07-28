"""
===============================================================
Visualizing the probabilistic predictions of a VotingClassifier
===============================================================

.. currentmodule:: sklearn

Plot the predicted class probabilities in a toy dataset predicted by three
different classifiers and averaged by the :class:`~ensemble.VotingClassifier`.

First, three linear classifiers are initialized. Two are spline models with
interaction terms, one using constant extrapolation and the other using periodic
extrapolation. The third classifier is a :class:`~kernel_approximation.Nystroem`
with the default "rbf" kernel.

In the first part of this example, these three classifiers are used to
demonstrate soft-voting using :class:`~ensemble.VotingClassifier` with weighted
average. We set `weights=[2, 1, 3]`, meaning the constant extrapolation spline
model's predictions are weighted twice as much as the periodic spline model's,
and the Nystroem model's predictions are weighted three times as much as the
periodic spline.

The second part demonstrates how soft predictions can be converted into hard
predictions.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# We first generate a noisy XOR dataset, which is a binary classification task.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

n_samples = 500
rng = np.random.default_rng(0)
feature_names = ["Feature #0", "Feature #1"]
common_scatter_plot_params = dict(
    cmap=ListedColormap(["tab:red", "tab:blue"]),
    edgecolor="white",
    linewidth=1,
)

xor = pd.DataFrame(
    np.random.RandomState(0).uniform(low=-1, high=1, size=(n_samples, 2)),
    columns=feature_names,
)
noise = rng.normal(loc=0, scale=0.1, size=(n_samples, 2))
target_xor = np.logical_xor(
    xor["Feature #0"] + noise[:, 0] > 0, xor["Feature #1"] + noise[:, 1] > 0
)

X = xor[feature_names]
y = target_xor.astype(np.int32)

fig, ax = plt.subplots()
ax.scatter(X["Feature #0"], X["Feature #1"], c=y, **common_scatter_plot_params)
ax.set_title("The XOR dataset")
plt.show()

# %%
# Due to the inherent non-linear separability of the XOR dataset, tree-based
# models would often be preferred. However, appropriate feature engineering
# combined with a linear model can yield effective results, with the added
# benefit of producing better-calibrated probabilities for samples located in
# the transition regions affected by noise.
#
# We define and fit the models on the whole dataset.

from sklearn.ensemble import VotingClassifier
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, SplineTransformer, StandardScaler

clf1 = make_pipeline(
    SplineTransformer(degree=2, n_knots=2),
    PolynomialFeatures(interaction_only=True),
    LogisticRegression(C=10),
)
clf2 = make_pipeline(
    SplineTransformer(
        degree=2,
        n_knots=4,
        extrapolation="periodic",
        include_bias=True,
    ),
    PolynomialFeatures(interaction_only=True),
    LogisticRegression(C=10),
)
clf3 = make_pipeline(
    StandardScaler(),
    Nystroem(gamma=2, random_state=0),
    LogisticRegression(C=10),
)
weights = [2, 1, 3]
eclf = VotingClassifier(
    estimators=[
        ("constant splines model", clf1),
        ("periodic splines model", clf2),
        ("nystroem model", clf3),
    ],
    voting="soft",
    weights=weights,
)

clf1.fit(X, y)
clf2.fit(X, y)
clf3.fit(X, y)
eclf.fit(X, y)

# %%
# Finally we use :class:`~inspection.DecisionBoundaryDisplay` to plot the
# predicted probabilities. By using a diverging colormap (such as `"RdBu"`), we
# can ensure that darker colors correspond to `predict_proba` close to either 0
# or 1, and white corresponds to `predict_proba` of 0.5.

from itertools import product

from sklearn.inspection import DecisionBoundaryDisplay

fig, axarr = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(10, 8))
for idx, clf, title in zip(
    product([0, 1], [0, 1]),
    [clf1, clf2, clf3, eclf],
    [
        "Splines with\nconstant extrapolation",
        "Splines with\nperiodic extrapolation",
        "RBF Nystroem",
        "Soft Voting",
    ],
):
    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        response_method="predict_proba",
        plot_method="pcolormesh",
        cmap="RdBu",
        alpha=0.8,
        ax=axarr[idx[0], idx[1]],
    )
    axarr[idx[0], idx[1]].scatter(
        X["Feature #0"],
        X["Feature #1"],
        c=y,
        **common_scatter_plot_params,
    )
    axarr[idx[0], idx[1]].set_title(title)
    fig.colorbar(disp.surface_, ax=axarr[idx[0], idx[1]], label="Probability estimate")

plt.show()

# %%
# As a sanity check, we can verify for a given sample that the probability
# predicted by the :class:`~ensemble.VotingClassifier` is indeed the weighted
# average of the individual classifiers' soft-predictions.
#
# In the case of binary classification such as in the present example, the
# :term:`predict_proba` arrays contain the probability of belonging to class 0
# (here in red) as the first entry, and the probability of belonging to class 1
# (here in blue) as the second entry.

test_sample = pd.DataFrame({"Feature #0": [-0.5], "Feature #1": [1.5]})
predict_probas = [est.predict_proba(test_sample).ravel() for est in eclf.estimators_]
for (est_name, _), est_probas in zip(eclf.estimators, predict_probas):
    print(f"{est_name}'s predicted probabilities: {est_probas}")

# %%
print(
    "Weighted average of soft-predictions: "
    f"{np.dot(weights, predict_probas) / np.sum(weights)}"
)

# %%
# We can see that manual calculation of predicted probabilities above is
# equivalent to that produced by the `VotingClassifier`:

print(
    "Predicted probability of VotingClassifier: "
    f"{eclf.predict_proba(test_sample).ravel()}"
)

# %%
# To convert soft predictions into hard predictions when weights are provided,
# the weighted average predicted probabilities are computed for each class.
# Then, the final class label is then derived from the class label with the
# highest average probability, which corresponds to the default threshold at
# `predict_proba=0.5` in the case of binary classification.

print(
    "Class with the highest weighted average of soft-predictions: "
    f"{np.argmax(np.dot(weights, predict_probas) / np.sum(weights))}"
)

# %%
# This is equivalent to the output of `VotingClassifier`'s `predict` method:

print(f"Predicted class of VotingClassifier: {eclf.predict(test_sample).ravel()}")

# %%
# Soft votes can be thresholded as for any other probabilistic classifier. This
# allows you to set a threshold probability at which the positive class will be
# predicted, instead of simply selecting the class with the highest predicted
# probability.

from sklearn.model_selection import FixedThresholdClassifier

eclf_other_threshold = FixedThresholdClassifier(
    eclf, threshold=0.7, response_method="predict_proba"
).fit(X, y)
print(
    "Predicted class of thresholded VotingClassifier: "
    f"{eclf_other_threshold.predict(test_sample)}"
)
