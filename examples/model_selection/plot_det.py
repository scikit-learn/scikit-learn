"""
=================================
Comparing ROC, DET and CAP curves
=================================

In this example, we compare three binary classification multi-threshold metrics:
the Receiver Operating Characteristic (ROC), Detection Error Tradeoff (DET)
and the Cumulative Accuracy Profile (CAP). For such purpose, we evaluate
two different classifiers for the same classification task.

ROC curves feature true positive rate (TPR) on the Y axis, and false positive
rate (FPR) on the X axis. This means that the top left corner of the plot is the
"ideal" point - a FPR of zero, and a TPR of one.

DET curves are a variation of ROC curves where False Negative Rate (FNR) is
plotted on the y-axis instead of the TPR. In this case the origin (bottom left
corner) is the "ideal" point.

CAP curves display the cumulative proportion of true positives on the Y axis
versus the cumulative proportion of the dataset (ranked by predicted probability)
on the X axis. In contrast to ROC and DET curves, the "ideal" classifier does not
correspond to a point, but a curve in the diagram, where each examined case perfectly
leads to a true positive.

.. note::

    - See :func:`sklearn.metrics.roc_curve` for further information about ROC
      curves.

    - See :func:`sklearn.metrics.det_curve` for further information about
      DET curves.

    - This example is loosely based on
      :ref:`sphx_glr_auto_examples_classification_plot_classifier_comparison.py`
      example.

    - See :ref:`sphx_glr_auto_examples_model_selection_plot_roc_crossval.py` for
      an example estimating the variance of the ROC curves and ROC-AUC.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Generate synthetic data
# -----------------------

import numpy as np

from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

X, y = make_classification(
    n_samples=1_000,
    n_features=2,
    n_redundant=0,
    n_informative=2,
    random_state=1,
    n_clusters_per_class=1,
)
# Use string labels for the classes to illustrate the use of the `pos_label` parameter
# to identify which class should be considered as the "positive" class.
y = np.array(["Class A", "Class B"])[y]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)

# %%
# Define the classifiers
# ----------------------
#
# Here we define two different classifiers. The goal is to visually compare their
# statistical performance across thresholds using the ROC, CAP and DET curves. There
# is no particular reason why these classifiers are chosen over other classifiers
# available in scikit-learn.

from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC

classifiers = {
    "Linear SVM": make_pipeline(StandardScaler(), LinearSVC(C=0.025)),
    "Random Forest": RandomForestClassifier(
        max_depth=5, n_estimators=10, max_features=1, random_state=0
    ),
    "Non-informative baseline": DummyClassifier(),
}

# %%
# Compare ROC, DET and CAP curves
# -------------------------------
#
# DET curves are commonly plotted in normal deviate scale. To achieve this the
# DET display transforms the error rates as returned by the
# :func:`~sklearn.metrics.det_curve` and the axis scale using
# `scipy.stats.norm`.

import matplotlib.pyplot as plt

from sklearn.metrics import CAPCurveDisplay, DetCurveDisplay, RocCurveDisplay

fig, [ax_roc, ax_det, ax_cap] = plt.subplots(
    1, 3, figsize=(15, 5), constrained_layout=True
)
pos_label = "Class A"
for clf_idx, (name, clf) in enumerate(classifiers.items()):
    is_last = clf_idx == len(classifiers) - 1
    clf.fit(X_train, y_train)

    RocCurveDisplay.from_estimator(
        clf,
        X_test,
        y_test,
        ax=ax_roc,
        name=name,
        pos_label=pos_label,
        plot_chance_level=is_last,
    )
    DetCurveDisplay.from_estimator(
        clf,
        X_test,
        y_test,
        ax=ax_det,
        name=name,
        pos_label=pos_label,
    )
    CAPCurveDisplay.from_estimator(
        clf,
        X_test,
        y_test,
        ax=ax_cap,
        name=name,
        pos_label=pos_label,
        plot_perfect=is_last,
        plot_chance_level=is_last,
    )

ax_roc.set_title("Receiver Operating Characteristic (ROC) curves")
ax_det.set_title("Detection Error Tradeoff (DET) curves")
ax_cap.set_title("Cumulative Accuracy Profile (CAP) curves")

ax_roc.grid(linestyle="--")
ax_det.grid(linestyle="--")
ax_cap.grid(linestyle="--")

plt.legend()
plt.show()

# %% Analysis
# --------
#
# All curves agree that the Random Forest classifier has more discriminative
# power than the Linear SVM:
#
# - the area under the ROC and CAP curves is larger for the Random Forest
#   classifier,
# - the ROC and CAP curves are uniformly closer to the top left corner for the
#   Random Forest classifier,
# - the DET curve is uniformly closer to the origin for the Random Forest
#   classifier.
#
# Notice that it is easier to visually assess that the Random Forest classifier
# performs better than the Linear SVM classifier using DET curves. As ROC and
# CAP curves are plotted with a linear scale, different classifiers usually appear
# similar for a large part of the plot and differ the most in the top left
# corner of the graph. On the other hand, because DET curves represent straight
# lines in normal deviate scale, they tend to be distinguishable as a whole and
# the area of interest spans a large part of the plot.
#
# DET curves give direct feedback of the detection error tradeoff to aid in
# operating point analysis. The user can then decide the FNR they are willing
# to accept at the expense of the FPR (or vice-versa).
#
# Non-informative classifier baseline
# -----------------------------------
#
# The diagonal black-dotted lines named "chance level" in the plots above
# correspond to a the expected value of a non-informative classifier on an
# infinite evaluation set.
#
# The :class:`~sklearn.dummy.DummyClassifier` model makes constant predictions,
# independently of the input features in `X`, making it a canonical example of
# such a non-informative classifier. We observe that the ROC and DET curves of
# the non-informative classifier exactly match the theoretical chance level
# line. This can be explained as follows. First recall the following
# mathematical definitions:
#
# :math:`\text{FPR} = \frac{\text{FP}}{\text{FP} + \text{TN}}`
#
# :math:`\text{FNR} = \frac{\text{FN}}{\text{TP} + \text{FN}}`
#
# :math:`\text{TPR} = \frac{\text{TP}}{\text{TP} + \text{FN}}`
#
# A classifier that always predict the positive class would have no true
# negatives nor false negatives, giving :math:`\text{FPR} = \text{TPR} = 1` and
# :math:`\text{FNR} = 0`, i.e.:
#
# - a single point in the upper right corner of the ROC plane,
# - a single point in the lower right corner of the DET plane.
#
# Similarly, a classifier that always predict the negative class would have no
# true positives nor false positives, thus :math:`\text{FPR} = \text{TPR} = 0`
# and :math:`\text{FNR} = 1`, i.e.:
#
# - a single point in the lower left corner of the ROC plane,
# - a single point in the upper left corner of the DET plane.
#
# For the non-informative classifier with strategy "prior", the value returned
# by `predict_proba` is the observed frequency of the positive class in the
# training set. If the threshold is above this value, the classifier always
# predicts the negative class, and if the threshold is below this value, the
# classifier always predicts the positive class. Both ROC and DET curves
# linearly interpolate between these two points, hence the diagonal line.
#
# The CAP curve of the non-informative classifier only approximately lies
# on the diagonal. The match would be better in the limit of an infinite
# evaluation set.
