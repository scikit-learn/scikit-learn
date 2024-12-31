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
# Here we run a :class:`~sklearn.ensemble.HistGradientBoostingClassifier`
# classifier with cross-validation and plot the ROC curves fold-wise. Notice
# that the baseline to define the chance level (dashed ROC curve) is a
# classifier that would always predict the most frequent class.
#
# In the following plot, quantile coverage is represented in grey, though the
# AUC value is reported in terms of the mean and standar deviation.

import matplotlib.pyplot as plt
import numpy as np

from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import RocCurveDisplay, auc
from sklearn.model_selection import StratifiedShuffleSplit

n_splits = 30
cv = StratifiedShuffleSplit(n_splits=n_splits, random_state=0)
classifier = HistGradientBoostingClassifier(random_state=42)

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

fig, ax = plt.subplots(figsize=(6, 6))
for fold, (train, test) in enumerate(cv.split(X, y)):
    classifier.fit(X[train], y[train])
    viz = RocCurveDisplay.from_estimator(
        classifier,
        X[test],
        y[test],
        label=None,
        alpha=0.3,
        lw=1,
        ax=ax,
        plot_chance_level=(fold == n_splits - 1),
        chance_level_kw={"label": None},
    )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="b",
    label=rf"Mean ROC (AUC = {mean_auc:.2f} $\pm$ {std_auc:.2f})",
    lw=2,
    alpha=0.8,
)


upper_quantile = np.quantile(tprs, 0.95, axis=0)
lower_quantile = np.quantile(tprs, 0.05, axis=0)
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
