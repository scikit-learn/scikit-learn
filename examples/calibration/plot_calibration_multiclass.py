"""
==================================================
Probability Calibration for 3-class Classification
==================================================

This example illustrates how :ref:`sigmoid calibration <sigmoid_regressor>` and
:ref:`temperature scaling <temperature_scaling>` change
predicted probabilities for a 3-class classification problem. Illustrated is
the standard 2-simplex, where the three corners correspond to the three
classes. Arrows point from the probability vectors predicted by an uncalibrated
classifier to the probability vectors predicted by the same classifier after
sigmoid calibration on a hold-out validation set. Colors indicate the true
class of an instance (red: class 1, green: class 2, blue: class 3).

"""

# %%
# Generating the Dataset and Training the Classifier and Calibrators
# -------------------------------------------------------------------
# Below, we generate a classification dataset with 2000 samples, 2 features
# and 3 target classes. We then split the data as follows:
#
# * train: 600 samples (for training the classifier)
# * valid: 400 samples (for calibrating predicted probabilities)
# * test: 1000 samples
#
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.datasets import make_blobs

np.random.seed(0)

X, y = make_blobs(
    n_samples=2000, n_features=2, centers=3, random_state=42, cluster_std=5.0
)
X_train, y_train = X[:600], y[:600]
X_valid, y_valid = X[600:1000], y[600:1000]
X_test, y_test = X[1000:], y[1000:]

# %%
# To train the calibrated classifier, we start with the same
# :class:`~sklearn.ensemble.RandomForestClassifier` but train it using only
# the train data subset (600 samples):

from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.frozen import FrozenEstimator

clf = RandomForestClassifier(n_estimators=25)
clf.fit(X_train, y_train)

# %%
# Then calibrate, with sigmoid calibration
# `method="sigmoid"` and temperature scaling `method="temperature"` separately, using
# the valid data subset (400 samples) in a 2-stage process.
#
# Sigmoid calibration:

cal_clf_sig = CalibratedClassifierCV(FrozenEstimator(clf), method="sigmoid")
cal_clf_sig.fit(X_valid, y_valid)

# %%
# Temperature scaling:

cal_clf_ts = CalibratedClassifierCV(FrozenEstimator(clf), method="temperature")
cal_clf_ts.fit(X_valid, y_valid)

# %%
# Comparing Probabilities Before and After Calibration
# ----------------------------------------------------
# To compare the probabilities before and after calibration on the testing set,
# we plot a probability vector field on a 2-simplex using the following helper
# functions:


def _plot_simplex_boundary(ax):
    """Plot the boundary of the standard 2-simplex"""
    ax.plot([0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], "k", label="Simplex")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-0.25, 1.1)
    ax.set_ylim(-0.1, 1.1)

    ax.set_xlabel("Probability of Class 1", fontsize=14)
    ax.set_ylabel("Probability of Class 2", fontsize=14)


def _plot_perfect_predictions(ax):
    """Plot perfect predictions for each class"""
    ax.plot([1.0], [0.0], "ro", ms=10, label="Class 1")
    ax.plot([0.0], [1.0], "go", ms=10, label="Class 2")
    ax.plot([0.0], [0.0], "bo", ms=10, label="Class 3")

    labels = [
        ((-0.08, -0.08), r"$\left( 0, 0, 1 \right)$"),
        ((0.9, -0.08), r"$\left( 1, 0, 0 \right)$"),
        ((-0.08, 1.03), r"$\left( 0, 1, 0 \right)$"),
    ]
    for (lx, ly), label in labels:
        ax.text(lx, ly, label, fontsize=10)


def _annotate_midpoints(ax):
    """Plot the edge mid-points and the mid-point of the simplex"""

    # Edge mid-points
    xs = [1.0 / 3, 0.5, 0.0, 0.5]
    ys = [1.0 / 3, 0.0, 0.5, 0.5]
    ax.scatter(xs, ys, s=100, c="black")

    # Edge mid-point labels
    labels = [
        ((0.25, 0.17), r"$\left( \frac{1}{3}, \frac{1}{3}, \frac{1}{3} \right)$"),
        ((0.4, -0.08), r"$\left( \frac{1}{2}, 0, \frac{1}{2} \right)$"),
        ((-0.23, 0.5), r"$\left( 0, \frac{1}{2}, \frac{1}{2} \right)$"),
        ((0.55, 0.5), r"$\left( \frac{1}{2}, \frac{1}{2}, 0 \right)$"),
    ]
    for (lx, ly), label in labels:
        ax.text(lx, ly, label, fontsize=10)

    # Plot the mid-point
    ax.plot([1.0 / 3], [1.0 / 3], "ko", ms=5)


def _add_grid(ax):
    """Add grid to the simplex"""
    ax.grid(False)

    for x in np.linspace(0.0, 1.0, num=10):
        ax.plot([0, x], [x, 0], "k", alpha=0.2)
        ax.plot([0, 0 + (1 - x) / 2], [x, x + (1 - x) / 2], "k", alpha=0.2)
        ax.plot([x, x + (1 - x) / 2], [0, 0 + (1 - x) / 2], "k", alpha=0.2)


def plot_simplex(ax, annotate=False):
    """Plot the standard 2-simplex"""

    _plot_simplex_boundary(ax)

    if annotate:
        _plot_perfect_predictions(ax)
        _annotate_midpoints(ax)

    _add_grid(ax)

    ax.legend()


def plot_calibration_shifts(ax, clf_probs, cal_clf_probs, y_test):
    """Plot predicted probability shifts induced by calibration."""
    colors = ["r", "g", "b"]

    for i in range(clf_probs.shape[0]):
        ax.arrow(
            clf_probs[i, 0],
            clf_probs[i, 1],
            cal_clf_probs[i, 0] - clf_probs[i, 0],
            cal_clf_probs[i, 1] - clf_probs[i, 1],
            color=colors[y_test[i]],
            head_width=1e-2,
        )


def plot_calibrated_predictor(ax, cal_clf):
    """Plot the calibrated predictor map on the simplex."""
    p1d = np.linspace(0, 1, 20)
    p0, p1 = np.meshgrid(p1d, p1d)
    p2 = 1 - p0 - p1
    p = np.c_[p0.ravel(), p1.ravel(), p2.ravel()]
    p = p[p[:, 2] >= 0]

    # Compute calibrated probabilities
    calibrated_classifier = cal_clf.calibrated_classifiers_[0]
    if cal_clf.method == "temperature":
        prediction = calibrated_classifier.calibrators[0].predict(p)
    else:
        # Use the three class-wise calibrators to compute calibrated probabilities
        prediction = np.vstack(
            [
                calibrator.predict(this_p)
                for calibrator, this_p in zip(calibrated_classifier.calibrators, p.T)
            ]
        ).T

    # Re-normalize the calibrated predictions to make sure they stay inside the
    # simplex. This same renormalization step is performed internally by the
    # predict method of CalibratedClassifierCV on multiclass problems.
    prediction /= prediction.sum(axis=1)[:, None]

    # Plot changes in predicted probabilities induced by the calibrators
    colors = ["r", "g", "b"]

    for i in range(prediction.shape[0]):
        ax.arrow(
            p[i, 0],
            p[i, 1],
            prediction[i, 0] - p[i, 0],
            prediction[i, 1] - p[i, 1],
            head_width=1e-2,
            color=colors[np.argmax(p[i])],
        )


# %%
# Probability Flow From Predicted Probabilities to Calibrated Probabilities
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# We now plot the probability vector field on the 2-simplex to show how different
# calibration method moves the predicted probabilities among classes.

import matplotlib.pyplot as plt

clf_probs = clf.predict_proba(X_test)
cal_clf_probs_sig = cal_clf_sig.predict_proba(X_test)
cal_clf_probs_ts = cal_clf_ts.predict_proba(X_test)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

plot_simplex(axes[0], annotate=True)
plot_calibration_shifts(axes[0], clf_probs, cal_clf_probs_sig, y_test)

plot_simplex(axes[1], annotate=True)
plot_calibration_shifts(axes[1], clf_probs, cal_clf_probs_ts, y_test)

axes[0].legend().remove()

axes[0].set_title("Sigmoid Calibration", fontsize=15)
axes[1].set_title("Temperature Scaling", fontsize=15)
fig.suptitle(
    "Calibration of Random Forest Predicted Probabilities on the Simplex", fontsize=18
)

_ = plt.show()

# %%
# In the figure above:
#
# * Each vertex of the simplex represents a perfectly predicted class (e.g., 1, 0, 0).
# * The mid point inside the simplex represents predicting the three classes with equal
#   probability (i.e., 1/3, 1/3, 1/3).
# * Each arrow starts at the predicted probabilities and end with the arrow head at the
#   calibrated probability.
# * The color of the arrow represents the true class of that test sample.
#
# The uncalibrated classifier is overly confident in its predictions and
# incurs a large :ref:`log loss <log_loss>`, which we will see below.
#
# For sigmoid calibration, the calibrated classifier incurs a lower
# :ref:`log loss <log_loss>` due to two factors. First, notice in the
# sigmoid calibration figure that the arrows generally point away from the edges of the
# simplex, where the probability of one class is 0. Second, a large proportion
# of the arrows point towards the true class, e.g., green arrows (samples where
# the true class is 'green') generally point towards the green vertex. This
# results in fewer over-confident, 0 predicted probabilities and at the same
# time an increase in the predicted probabilities of the correct class.
# Thus, the calibrated classifier produces more accurate predicted probabilities
# that incur a lower :ref:`log loss <log_loss>`
#
#
# Log Loss and Brier Score
# ^^^^^^^^^^^^^^^^^^^^^^^^
# We can show this objectively by comparing the :ref:`log loss <log_loss>` of
# the uncalibrated and calibrated classifiers on the predictions of the 1000
# test samples. Note that an alternative would have been to increase the number
# of base estimators (trees) of the
# :class:`~sklearn.ensemble.RandomForestClassifier` which would have resulted
# in a similar decrease in :ref:`log loss <log_loss>`.

from sklearn.metrics import log_loss

loss = log_loss(y_test, clf_probs)
cal_loss_sig = log_loss(y_test, cal_clf_probs_sig)
cal_loss_ts = log_loss(y_test, cal_clf_probs_ts)

print("Log-loss of:")
print(f" - uncalibrated classifier: {loss:.3f}")
print(f" - sigmoid calibrated classifier: {cal_loss_sig:.3f}")
print(f" - calibrated classifier by temperature scaling: {cal_loss_ts:.3f}")

# %%
# We can also assess calibration with the Brier score for probabilistics predictions
# (lower is better, possible range is [0, 2]):

from sklearn.metrics import brier_score_loss

loss = brier_score_loss(y_test, clf_probs)
cal_loss_sig = brier_score_loss(y_test, cal_clf_probs_sig)
cal_loss_ts = brier_score_loss(y_test, cal_clf_probs_ts)

print("Brier score of")
print(f" - uncalibrated classifier: {loss:.3f}")
print(f" - sigmoid calibrated classifier: {cal_loss_sig:.3f}")
print(f" - calibrated classifier by temperature scaling: {cal_loss_ts:.3f}")

# %%
# According to the Brier score, the calibrated classifier is not better than
# the original model.
#
# Probability Flow From Arbitrary Probabilities to Calibrated Probabilities
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Finally we generate a grid of possible uncalibrated probabilities over
# the 2-simplex, compute the corresponding calibrated probabilities and
# plot arrows for each. The arrows are colored according the highest
# uncalibrated probability. This illustrates the learned calibration map:

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

plot_simplex(axes[0])
plot_calibrated_predictor(axes[0], cal_clf_sig)

plot_simplex(axes[1])
plot_calibrated_predictor(axes[1], cal_clf_ts)

axes[0].legend().remove()

axes[0].set_title("Sigmoid Calibration", fontsize=15)
axes[1].set_title("Temperature Scaling", fontsize=15)
fig.suptitle("Calibration of Arbitrary Probabilities on the Simplex", fontsize=18)

_ = plt.show()

# %%
# One can observe that, on average, the calibrator is pushing highly confident
# predictions away from the boundaries of the simplex while simultaneously
# moving uncertain predictions towards one of three modes, one for each class.
# We can also observe that the mapping is not symmetric. Furthermore some
# arrows seem to cross class assignment boundaries which is not necessarily
# what one would expect from a calibration map as it means that some predicted
# classes will change after calibration.
#
# All in all, the One-vs-Rest multiclass-calibration strategy implemented in
# `CalibratedClassifierCV` should not be trusted blindly.
