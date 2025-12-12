"""
==============================
Probability Calibration curves
==============================

When performing classification, one often wants to predict not only the class
label, but also the associated probability. This probability gives a sense of confidence
on the prediction. This example demonstrates how to visualize how well-calibrated the
predicted probabilities are using calibration curves, also known as reliability
diagrams. We also show how to calibrate an initially uncalibrated classifier.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Dataset
# -------
#
# We will use a synthetic multi-class classification dataset with 3 classes,
# 100,000 samples and 20 features. Of these features, 2 are informative, 10 are
# redundant (random combinations of the informative ones) and the
# remaining 8 are uninformative (random noise). Of the 100,000 samples, 1,000
# are used for model fitting and the rest for testing.
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split

X, y = make_classification(
    n_samples=100_000,
    n_features=20,
    n_informative=2,
    n_redundant=10,
    n_classes=3,
    n_clusters_per_class=1,
    random_state=42,
)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.99, random_state=42
)

# %%
# Reliability Diagrams
# --------------------
#
# Gaussian Naive Bayes
# ^^^^^^^^^^^^^^^^^^^^
#
# First, we will compare:
#
# * :class:`~sklearn.linear_model.LogisticRegression` (used as baseline
#   since very often, properly regularized logistic regression is well
#   calibrated by default thanks to the use of the log-loss)
# * Uncalibrated :class:`~sklearn.naive_bayes.GaussianNB`
# * :class:`~sklearn.naive_bayes.GaussianNB` with isotonic calibration, sigmoid
#   calibration, and temperature scaling (see :ref:`User Guide <calibration>`).
import matplotlib.pyplot as plt

from sklearn.calibration import CalibratedClassifierCV, CalibrationDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB

lr = LogisticRegression(C=1.0)
gnb = GaussianNB()
gnb_isotonic = CalibratedClassifierCV(gnb, cv=2, method="isotonic", ensemble=False)
gnb_sigmoid = CalibratedClassifierCV(gnb, cv=2, method="sigmoid", ensemble=False)
gnb_ts = CalibratedClassifierCV(gnb, cv=2, method="temperature", ensemble=False)

clf_list = [
    (lr, "Logistic"),
    (gnb, "Naive Bayes"),
    (gnb_isotonic, "Naive Bayes + Isotonic"),
    (gnb_sigmoid, "Naive Bayes + Sigmoid"),
    (gnb_ts, "Naive Bayes + Temperature Scaling"),
]

for clf, _ in clf_list:
    clf.fit(X_train, y_train)

classes = clf_list[0][0].classes_
n_classes = len(classes)
n_clf = len(clf_list)

# %%
# Reliability diagrams for each class under all five conditions are plotted below.
# The x-axis shows the average predicted probability for each bin, and the y-axis
# shows the fraction of samples belonging to the class in that bin.
fig, axes = plt.subplots(figsize=(20, 6), ncols=n_classes)

for class_idx, c in enumerate(classes):
    ax = axes[class_idx]

    for clf, name in clf_list:
        CalibrationDisplay.from_predictions(
            y_true=(y_test == c),
            y_prob=clf.predict_proba(X_test)[:, class_idx],
            n_bins=10,
            name=name,
            ax=ax,
        )
    ax.set_title(f"Class {c}")
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.grid(True)
    if c < n_classes - 1:
        ax.get_legend().remove()
    else:
        ax.legend(loc=(1.05, 0.5))

fig.suptitle("Reliability Diagrams per Class (Gaussian Naive Bayes)", fontsize=18)
fig.supxlabel("Mean Predicted Probability")
fig.supylabel("Fraction of Class")
plt.tight_layout()
plt.subplots_adjust(left=0.05)
_ = plt.show()

# %%
# Uncalibrated :class:`~sklearn.naive_bayes.GaussianNB` is poorly-calibrated
# because of the redundant features violates the assumption of feature independence,
# leading to an overly-confident classifier. This is reflected in the typical
# transposed-sigmoid curve.
#
# Calibrating the probabilities of :class:`~sklearn.naive_bayes.GaussianNB` with
# :ref:`isotonic regression <isotonic_regressor>` can fix this issue, as seen
# in the nearly diagonal calibration curve.
#
# :ref:`Sigmoid regression <sigmoid_regressor>` also improves calibration,
# but not as strongly as the non-parametric
# :ref:`isotonic regression <isotonic_regressor>`.
# This can be attributed to the fact that we have plenty of calibration data,
# allowing the more flexibile non-parametric model to perform better.
#
# :ref:`Temperature scaling <temperature_scaling>` offers a lightweight, parametric
# alternative that preserves probability rankings via monotone transformation, while
# pushing them away from 0.0 and 1.0. On this dataset, it has the best calibration
# improvement among the methods.
#
# We also plot the distribution of predicted probabilities below as histograms to
# illustrate how the calibration methods differ in practice.
fig, axes = plt.subplots(
    figsize=(20, 10),
    nrows=n_clf,
    ncols=n_classes,
)

colors = plt.get_cmap("Dark2")
lines = []
labels = []

for col in range(n_classes):
    for row, (clf, name) in enumerate(clf_list):
        ax = axes[row, col]

        probs = clf.predict_proba(X_test)[:, col]

        ax.hist(
            probs,
            range=(0, 1),
            bins=10,
            label=name,
            color=colors(row),
        )

        ax.set_xlabel(None)
        ax.set_ylabel(None)

        if col == 0:
            Line, Label = ax.get_legend_handles_labels()
            lines.extend(Line)
            labels.extend(Label)


for col, c in enumerate(classes):
    axes[0, col].set_title(f"Class {c}")

fig.supxlabel("Mean Predicted Probability")
fig.supylabel("Count")
fig.suptitle("Predicted Probability Distributions (Gaussian Naive Bayes)", fontsize=18)
fig.legend(lines, labels, loc="center right")
plt.subplots_adjust(right=0.82)
_ = plt.show()

# %%
# :class:`~sklearn.naive_bayes.GaussianNB` predictions are heavily concentrated near
# 0.0 and 1.0; whereas the calibrated variants, particularly
# :ref:`isotonic regression <isotonic_regressor>`, redistribute probability mass,
# and produce distributions that more closely match empirical frequencies.
#
# Below we analyse several classification metrics: :ref:`brier_score_loss`,
# :ref:`log_loss`,
# :ref:`precision, recall, F1 score <precision_recall_f_measure_metrics>`, and
# :ref:`ROC AUC <roc_metrics>`.
from collections import defaultdict
from functools import partial

import pandas as pd

from sklearn.metrics import (
    brier_score_loss,
    f1_score,
    log_loss,
    precision_score,
    recall_score,
    roc_auc_score,
)

scores = defaultdict(list)
for i, (clf, name) in enumerate(clf_list):
    clf.fit(X_train, y_train)
    y_prob = clf.predict_proba(X_test)
    y_pred = clf.predict(X_test)
    scores["Classifier"].append(name)

    for metric, score_name in [
        (brier_score_loss, "Brier Score"),
        (log_loss, "Log Loss"),
        (partial(roc_auc_score, multi_class="ovr"), "ROC AUC"),
    ]:
        scores[score_name].append(metric(y_test, y_prob))

    for metric, score_name in [
        (partial(precision_score, average="weighted"), "Precision"),
        (partial(recall_score, average="weighted"), "Recall"),
        (partial(f1_score, average="weighted"), "F1"),
    ]:
        scores[score_name].append(metric(y_test, y_pred))

score_df = pd.DataFrame(scores).set_index("Classifier")
score_df = score_df.round(decimals=3)

score_df

# %%
# Note that although calibration improves the :ref:`brier_score_loss` (which includes
# both a calibration and refinement term) and :ref:`log_loss`,
# it does not significantly change the prediction accuracy metrics
# (:ref:`precision, recall, F1 score <precision_recall_f_measure_metrics>`).
# This is because calibration should not significantly change prediction
# probabilities at the location of the decision threshold (at x = 0.5).
# Instead, calibration aims to make the predicted probabilities more
# accurate and thus, more useful for allocation decisions under
# uncertainty.
# Furthermore, :ref:`ROC AUC <roc_metrics>` should remain unchanged, because calibration
# applies a monotone transformation, and rank-based metrics are unaffected by such
# transformations.
#
# Linear Support Vector Classifier
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Next, we will compare:
#
# * :class:`~sklearn.linear_model.LogisticRegression` (baseline)
# * Uncalibrated :class:`~sklearn.svm.LinearSVC`. Since SVC does not output
#   probabilities by default, we naively scale the output of the
#   :term:`decision_function` into [0, 1] by applying min-max scaling.
# * :class:`~sklearn.svm.LinearSVC` with isotonic calibration, sigmoid
#   calibration, and temperature scaling (see :ref:`User Guide <calibration>`).
import numpy as np

from sklearn.svm import LinearSVC


class NaivelyCalibratedLinearSVC(LinearSVC):
    """LinearSVC with `predict_proba` method that naively scales
    `decision_function` output for binary classification."""

    def fit(self, X, y):
        super().fit(X, y)
        df = self.decision_function(X)
        self.df_min_ = df.min()
        self.df_max_ = df.max()

    def predict_proba(self, X):
        df = self.decision_function(X)
        df = np.atleast_2d(df)

        # Min-max scale per class
        scaled = (df - self.df_min_) / (self.df_max_ - self.df_min_)
        scaled = np.clip(scaled, 0, 1)

        # Normalise across classes so rows sum to 1
        row_sums = scaled.sum(axis=1, keepdims=True)
        # Avoid divide-by-zero (degenerate cases)
        row_sums[row_sums == 0] = 1
        proba = scaled / row_sums

        return proba


# %%

lr = LogisticRegression(C=1.0)
svc = NaivelyCalibratedLinearSVC(max_iter=10_000)
svc_isotonic = CalibratedClassifierCV(svc, cv=2, method="isotonic")
svc_sigmoid = CalibratedClassifierCV(svc, cv=2, method="sigmoid")
svc_ts = CalibratedClassifierCV(svc, cv=2, method="temperature")

clf_list = [
    (lr, "Logistic"),
    (svc, "SVC"),
    (svc_isotonic, "SVC + Isotonic"),
    (svc_sigmoid, "SVC + Sigmoid"),
    (svc_ts, "SVC + Temperature Scaling"),
]

for clf, _ in clf_list:
    clf.fit(X_train, y_train)

classes = clf_list[0][0].classes_
n_classes = len(classes)
n_clf = len(clf_list)

# %%
# The reliability diagram for each class is given below.
fig, axes = plt.subplots(figsize=(20, 6), ncols=n_classes)

for class_idx, c in enumerate(classes):
    ax = axes[class_idx]

    for clf, name in clf_list:
        CalibrationDisplay.from_predictions(
            y_true=(y_test == c),
            y_prob=clf.predict_proba(X_test)[:, class_idx],
            n_bins=10,
            name=name,
            ax=ax,
        )
    ax.set_title(f"Class {c}")
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.grid(True)
    if c < n_classes - 1:
        ax.get_legend().remove()
    else:
        ax.legend(loc=(1.05, 0.5))

fig.suptitle("Reliability Diagrams per Class (SVC)", fontsize=18)
fig.supxlabel("Mean Predicted Probability")
fig.supylabel("Fraction of Class")
plt.tight_layout()
plt.subplots_adjust(left=0.05)
_ = plt.show()

# %%
# We also plot the distribution of the predicted probabilities as histograms below,
# to demonstrate how these calibration methods behave differently.
fig, axes = plt.subplots(
    figsize=(20, 10),
    nrows=n_clf,
    ncols=n_classes,
)

colors = plt.get_cmap("Dark2")
lines = []
labels = []

for col in range(n_classes):
    for row, (clf, name) in enumerate(clf_list):
        ax = axes[row, col]

        probs = clf.predict_proba(X_test)[:, col]

        ax.hist(
            probs,
            range=(0, 1),
            bins=10,
            label=name,
            color=colors(row),
        )

        ax.set_xlabel(None)
        ax.set_ylabel(None)

        if col == 0:
            Line, Label = ax.get_legend_handles_labels()
            lines.extend(Line)
            labels.extend(Label)


for col, c in enumerate(classes):
    axes[0, col].set_title(f"Class {c}")

fig.supxlabel("Mean Predicted Probability")
fig.supylabel("Count")
fig.suptitle("Predicted Probability Distributions (SVC)", fontsize=18)
fig.legend(lines, labels, loc="center right")
plt.subplots_adjust(right=0.82)
_ = plt.show()

# %%
# :class:`~sklearn.svm.LinearSVC` shows the opposite
# behavior to :class:`~sklearn.naive_bayes.GaussianNB`. Its calibration curves have a
# sigmoid shape, which is typical for an under-confident classifier. This is caused by
# the margin property of the :ref:`hinge loss <_hinge_loss>`, which focuses on samples
# that are close to the decision boundary (support vectors). Samples that are far away
# from the decision boundary do not influence the :ref:`hinge loss <_hinge_loss>`.
# Therefore, the model does not attempt to separate samples in the high-confidence
# regions. This produces flatter calibration curves near 0 and 1, as shown empirically
# across many datasets in Niculescu-Mizil & Caruana [1]_.
#
# All three calibration methods can fix this and yield similar results.
# :ref:`Temperature scaling <temperature_scaling>` again produces the most diagonal
# calibration curves across all classes.
#
# As before, we show the :ref:`brier_score_loss`, :ref:`log_loss`,
# :ref:`precision, recall, F1 score <precision_recall_f_measure_metrics>`, and
# :ref:`ROC AUC <roc_metrics>`.

scores = defaultdict(list)
for i, (clf, name) in enumerate(clf_list):
    clf.fit(X_train, y_train)
    y_prob = clf.predict_proba(X_test)
    y_pred = clf.predict(X_test)
    scores["Classifier"].append(name)

    for metric, score_name in [
        (brier_score_loss, "Brier Score"),
        (log_loss, "Log Loss"),
        (partial(roc_auc_score, multi_class="ovr"), "ROC AUC"),
    ]:
        scores[score_name].append(metric(y_test, y_prob))

    for metric, score_name in [
        (partial(precision_score, average="weighted"), "Precision"),
        (partial(recall_score, average="weighted"), "Recall"),
        (partial(f1_score, average="weighted"), "F1"),
    ]:
        scores[score_name].append(metric(y_test, y_pred))

score_df = pd.DataFrame(scores).set_index("Classifier")
score_df = score_df.round(decimals=3)

score_df

# %%
# As with :class:`~sklearn.naive_bayes.GaussianNB`, calibration improves
# both :ref:`brier_score_loss` and :ref:`log_loss` for :class:`~sklearn.svm.LinearSVC`,
# but does not alter the prediction accuracy measures
# (:ref:`precision, recall, F1 score <precision_recall_f_measure_metrics>`) much.
#
# Summary
# -------
#
# Parametric :ref:`sigmoid calibration <sigmoid_regressor>` works well when the base
# classifier's calibration curve is sigmoid (e.g., for :class:`~sklearn.svm.LinearSVC`),
# but not where it is transposed-sigmoid
# (e.g., :class:`~sklearn.naive_bayes.GaussianNB`).
# Non-parametric :ref:`isotonic calibration <isotonic_regressor>` can handle both shapes
# but may require more data to perform well.
# :ref:`Temperature scaling <temperature_scaling>` handles both situations even with
# limited data, and often produce the best-calibrated probabilities.
#
# References
# ----------
#
# .. [1] `Predicting Good Probabilities with Supervised Learning
#        <https://dl.acm.org/doi/pdf/10.1145/1102351.1102430>`_,
#        A. Niculescu-Mizil & R. Caruana, ICML 2005
