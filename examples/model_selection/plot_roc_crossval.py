"""
=============================================================
Receiver Operating Characteristic (ROC) with cross validation
=============================================================

This example presents how to estimate and visualize the variance of the Receiver
Operating Characteristic (ROC) metric using cross-validation.

ROC curves typically feature true positive rate (TPR) on the Y axis, and false
positive rate (FPR) on the X axis. This means that the top left corner of the
plot is the "ideal" point - a FPR of zero, and a TPR of one. This is not very
realistic, but it does mean that a larger Area Under the Curve (AUC) is usually
better. The "steepness" of ROC curves is also important, since it is ideal to
maximize the TPR while minimizing the FPR.

This example demonstrates how the classifier's ROC response is influenced by
variations in the training data as obtained through ShuffleSplit cross-validation.
By analyzing all these curves, we can calculate the mean AUC and visualize the
variability of the estimated curves across CV splits via a quantile-based
region.

.. note::

    See :ref:`sphx_glr_auto_examples_model_selection_plot_roc.py` for a
    complement of the present example explaining the averaging strategies to
    generalize the metrics for multiclass classifiers.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Load and prepare data
# =====================
#
# We use :class:`~sklearn.datasets.make_classification` to generate a synthetic
# dataset with 1,000 samples. The generated dataset has two classes by default.
# In this case, we set a class separation factor of 0.5, making the classes
# partially overlapping and not perfectly linearly separable.

from sklearn.datasets import make_classification

X, y = make_classification(
    n_samples=1_000,
    n_features=2,
    n_redundant=0,
    n_informative=2,
    class_sep=0.5,
    random_state=0,
    n_clusters_per_class=1,
)

# %%
# Classification and ROC analysis
# -------------------------------
#
# Here we run :func:`~sklearn.model_selection.cross_validate` on a
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier`, then use the
# computed cross-validation results to plot the ROC curves fold-wise. Notice
# that the baseline to define the chance level (dashed ROC curve) is a
# classifier that would always predict the most frequent class.
#
# In the following plot, quantile coverage is represented in grey, though the
# AUC value is reported in terms of the mean and standard deviation.

import matplotlib.pyplot as plt
import numpy as np

from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import RocCurveDisplay, auc
from sklearn.model_selection import StratifiedShuffleSplit, cross_validate

n_splits = 30
cv = StratifiedShuffleSplit(n_splits=n_splits, random_state=0)
classifier = HistGradientBoostingClassifier(random_state=42)

cv_results = cross_validate(
    classifier, X, y, cv=cv, return_estimator=True, return_indices=True
)

prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]
curve_kwargs_list = [
    dict(
        alpha=0.3,
        lw=1,
        label=None,
        color=colors[fold % len(colors)],
    )
    for fold in range(n_splits)
]
names = [f"ROC fold {idx}" for idx in range(n_splits)]

mean_fpr = np.linspace(0, 1, 100)
interp_tprs = []

_, ax = plt.subplots(figsize=(6, 6))
viz = RocCurveDisplay.from_cv_results(
    cv_results,
    X,
    y,
    ax=ax,
    name=names,
    curve_kwargs=curve_kwargs_list,
    plot_chance_level=True,
)

for idx in range(n_splits):
    interp_tpr = np.interp(mean_fpr, viz.fpr[idx], viz.tpr[idx])
    interp_tpr[0] = 0.0
    interp_tprs.append(interp_tpr)

mean_tpr = np.mean(interp_tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(viz.roc_auc)

ax.plot(
    mean_fpr,
    mean_tpr,
    color="b",
    label=rf"Mean ROC (AUC = {mean_auc:.2f} $\pm$ {std_auc:.2f})",
    lw=2,
    alpha=0.8,
)


upper_quantile = np.quantile(interp_tprs, 0.95, axis=0)
lower_quantile = np.quantile(interp_tprs, 0.05, axis=0)
ax.fill_between(
    mean_fpr,
    lower_quantile,
    upper_quantile,
    color="grey",
    alpha=0.4,
    label="5% to 95% percentile region",
)

ax.set(
    xlabel="False Positive Rate",
    ylabel="True Positive Rate",
    title="Mean ROC curve with variability",
)
ax.legend(loc="lower right")
plt.show()
