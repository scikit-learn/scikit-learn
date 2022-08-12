"""
=======================================
Receiver Operating Characteristic (ROC)
=======================================

Example of Receiver Operating Characteristic (ROC) metric to evaluate classifier
output quality.

ROC curves typically feature true positive rate on the Y axis, and false
positive rate on the X axis. This means that the top left corner of the plot is
the "ideal" point - a false positive rate of zero, and a true positive rate of
one. This is not very realistic, but it does mean that a larger area under the
curve (AUC) is usually better.

The "steepness" of ROC curves is also important, since it is ideal to maximize
the true positive rate while minimizing the false positive rate.

ROC curves are typically used in binary classification to study the output of a
classifier. In order to extend ROC curve and ROC area to multi-label
classification, it is necessary to binarize the output. One ROC curve can be
drawn per label, but one can also draw a ROC curve by considering each element
of the label indicator matrix as a binary prediction (micro-averaging).

Another evaluation measure for multi-label classification is macro-averaging,
which gives equal weight to the classification of each label.

.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :ref:`sphx_glr_auto_examples_model_selection_plot_roc_crossval.py`

"""

# %%
# Load and prepare data
# =====================
#
# We import the :ref:`iris_dataset` which contains 3 classes, each one
# corresponding to a type of iris plant. One class is linearly separable from
# the other 2; the latter are NOT linearly separable from each other.
#
# Here we binarize the output and add noisy features to make the problem harder.

import numpy as np
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize

(X, y) = load_iris(return_X_y=True)

y = label_binarize(y, classes=[0, 1, 2])
n_classes = y.shape[1]

random_state = np.random.RandomState(0)
n_samples, n_features = X.shape
X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=0)

# %%
# Learn to predict each class against the other
# =============================================

from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc
from sklearn.multiclass import OneVsRestClassifier

classifier = OneVsRestClassifier(
    SVC(kernel="linear", probability=True, random_state=random_state)
)
y_score = classifier.fit(X_train, y_train).decision_function(X_test)


# %%
# Plot ROC curves
# ===============
#
# For a specific class
# --------------------

import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay

class_id = 2

RocCurveDisplay.from_predictions(
    y_test[:, class_id],
    y_score[:, class_id],
    name=f"ROC curve for {load_iris()['target_names'][class_id]}",
    color="darkorange",
)
plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.title(f"Receiver operating characteristic for class label {class_id}")
plt.show()

# %%
# ROC curve showing the "micro-average" of the 3 classes
# ------------------------------------------------------

RocCurveDisplay.from_predictions(
    y_test.ravel(), y_score.ravel(), name="One-vs-Rest SVC", color="darkorange"
)
plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.title("Receiver operating characteristic with micro-average")
plt.show()

# %%
# If we are not interested in the plot but only on the ROC-AUC itself, we can
# reproduce the values of the plot using `roc_auc_score`. By default, the
# computation of this score adds a single point at the maximal false positive
# rate by using linear interpolation and [McClish
# correction](https://pubmed.ncbi.nlm.nih.gov/2668680/).

from sklearn.metrics import roc_auc_score

y_prob = classifier.predict_proba(X_test)

micro_roc_auc_ovr = roc_auc_score(y_test, y_prob, average="micro")
weighted_roc_auc_ovr = roc_auc_score(
    y_test, y_prob, multi_class="ovr", average="weighted"
)

print(
    "One-vs-Rest ROC AUC scores:\n{:.6f} (micro),\n{:.6f} "
    "(weighted by prevalence)".format(micro_roc_auc_ovr, weighted_roc_auc_ovr)
)

# %%
# This is equivalent to computing the ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

# %%
# Plot ROC curves for the multiclass problem
# ..........................................
#
# Compute macro-average ROC curve and ROC area

# First aggregate all false positive rates
all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])

# Finally average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

# %%
# Plot all ROC curves
# -------------------

from itertools import cycle

plt.figure()
plt.plot(
    fpr["micro"],
    tpr["micro"],
    label="micro-average ROC curve (area = {0:0.2f})".format(roc_auc["micro"]),
    color="deeppink",
    linestyle=":",
    linewidth=4,
)

plt.plot(
    fpr["macro"],
    tpr["macro"],
    label="macro-average ROC curve (area = {0:0.2f})".format(roc_auc["macro"]),
    color="navy",
    linestyle=":",
    linewidth=4,
)

colors = cycle(["aqua", "darkorange", "cornflowerblue"])
for i, color in zip(range(n_classes), colors):
    plt.plot(
        fpr[i],
        tpr[i],
        color=color,
        lw=2,
        label="ROC curve of class {0} (area = {1:0.2f})".format(i, roc_auc[i]),
    )

plt.plot([0, 1], [0, 1], "k--", lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Some extension of Receiver operating characteristic to multiclass")
plt.legend(loc="lower right")
plt.show()

# %%
# Area under ROC for the multiclass problem
# .........................................
#
# The :func:`sklearn.metrics.roc_auc_score` function can be used for
# multi-class classification. The multi-class One-vs-One scheme compares every
# unique pairwise combination of classes. In this section, we calculate the AUC
# using the OvR and OvO schemes. We report a macro average, and a
# prevalence-weighted average.

y_prob = classifier.predict_proba(X_test)

macro_roc_auc_ovo = roc_auc_score(y_test, y_prob, multi_class="ovo", average="macro")
weighted_roc_auc_ovo = roc_auc_score(
    y_test, y_prob, multi_class="ovo", average="weighted"
)
macro_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
weighted_roc_auc_ovr = roc_auc_score(
    y_test, y_prob, multi_class="ovr", average="weighted"
)
print(
    "One-vs-One ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
    "(weighted by prevalence)".format(macro_roc_auc_ovo, weighted_roc_auc_ovo)
)
print(
    "One-vs-Rest ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
    "(weighted by prevalence)".format(macro_roc_auc_ovr, weighted_roc_auc_ovr)
)
