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
variations in the training data as obtained through K-fold cross-validation.
By analyzing all these curves, we can calculate the mean AUC and visualize the
variability of the estimated curves across CV folds via a quantile-based
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
# We import the :ref:`iris_dataset` which contains 3 classes, each one
# corresponding to a type of iris plant. One class is linearly separable from
# the other 2; the latter are **not** linearly separable from each other.
#
# In the following we binarize the dataset by dropping the "virginica" class
# (`class_id=2`). This means that the "versicolor" class (`class_id=1`) is
# regarded as the positive class and "setosa" as the negative class
# (`class_id=0`).

import numpy as np

from sklearn.datasets import load_iris

iris = load_iris()
target_names = iris.target_names
X, y = iris.data, iris.target
X, y = X[y != 2], y[y != 2]
n_samples, n_features = X.shape

# %%
# We also add noisy features to make the problem harder.
random_state = np.random.RandomState(0)
X = np.concatenate([X, random_state.randn(n_samples, 200 * n_features)], axis=1)

# %%
# Classification and ROC analysis
# -------------------------------
#
# Here we run a :class:`~sklearn.svm.SVC` classifier with cross-validation and
# plot the ROC curves fold-wise. Notice that the baseline to define the chance
# level (dashed ROC curve) is a classifier that would always predict the most
# frequent class.
#
# In the following plot, quantile coverage is represented in grey, though the
# AUC value is reported in terms of the mean and standar deviation.

import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.metrics import RocCurveDisplay, auc
from sklearn.model_selection import StratifiedKFold

n_splits = 6
cv = StratifiedKFold(n_splits=n_splits)
classifier = svm.SVC(kernel="linear", probability=True, random_state=random_state)

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
        name=f"ROC fold {fold}",
        alpha=0.3,
        lw=1,
        ax=ax,
        plot_chance_level=(fold == n_splits - 1),
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
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
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
