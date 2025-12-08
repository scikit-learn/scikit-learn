"""
================================================
Probability Calibration with Temperature Scaling
================================================

:class:`~sklearn.calibration.CalibratedClassifierCV` uses a one-vs-rest scheme to
calibrate multi-class probabilities. When using
:ref:`sigmoid regressor <sigmoid_regressor>` or
:ref:`isotonic regressor <isotonic_regressor>`, this means fitting
`n_classes` separate binary calibrators. This increases both the time and memory
required for calibration. In contrast,
:ref:`temperature scaling <temperature_scaling>` fits a single scalar parameter
and computes multi-class probabilities using the softmax function.

In this example, we illustrate how this simple one-parameter model can produce
better-calibrated probabilities in a multi-class setting.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Generate a Synthetic Dataset
# ----------------------------
# We begin by generating a synthetic dataset for multi-class classification.

from sklearn.datasets import make_classification

X, y = make_classification(n_classes=3, n_informative=8, random_state=42)

# %%
# Gaussian Naive Bayes and Calibration
# ------------------------------------
# Next, we apply :ref:`temperature scaling <temperature_scaling>` to calibrate
# the probabilities predicted by a :class:`~sklearn.naive_bayes.GaussianNB`, and
# compare the result with probabilities calibrated using the
# :ref:`sigmoid regressor <sigmoid_regressor>`.

from sklearn.calibration import CalibratedClassifierCV
from sklearn.naive_bayes import GaussianNB

gnb = GaussianNB().fit(X, y)
sig = CalibratedClassifierCV(gnb, method="sigmoid", ensemble=False).fit(X, y)
ts = CalibratedClassifierCV(gnb, method="temperature", ensemble=False).fit(X, y)

# %%
# As noted above, :ref:`sigmoid regressor <sigmoid_regressor>` fits one
# calibrator per class (three for this dataset), whereas
# :ref:`temperature scaling <temperature_scaling>` fits a single shared
# temperature parameter.

sig_calibrators_num = len(sig.calibrated_classifiers_[0].calibrators)
ts_calibrators_num = len(ts.calibrated_classifiers_[0].calibrators)

print(f"Number of calibrators in Sigmoid Calibration: {sig_calibrators_num}")
print(f"Number of calibrators in Temperature Scaling: {ts_calibrators_num}")

# %%
# The optimal temperature learned by :ref:`temperature scaling <temperature_scaling>`
# can be retrieved as:

print(
    f"Optimal temperature: {1 / ts.calibrated_classifiers_[0].calibrators[0].beta_:.3f}"
)

# %%
# Evaluation
# ----------
#
# Reliability Diagrams
# ^^^^^^^^^^^^^^^^^^^^
# To assess calibration, we plot a reliability diagram for each class.

import matplotlib.pyplot as plt

from sklearn.calibration import CalibrationDisplay

fig, axes = plt.subplots(figsize=(16, 9), ncols=3, sharey=True)
for i, c in enumerate(ts.classes_):
    CalibrationDisplay.from_predictions(
        y == c,
        gnb.predict_proba(X)[:, i],
        name="Gaussian Naive Bayes",
        ax=axes[i],
        linestyle=":",
        marker=".",
    )

    CalibrationDisplay.from_predictions(
        y == c,
        ts.predict_proba(X)[:, i],
        name="Temperature scaling",
        ax=axes[i],
        n_bins=10,
        marker="o",
    )

    CalibrationDisplay.from_predictions(
        y == c,
        sig.predict_proba(X)[:, i],
        name="Sigmoid",
        ax=axes[i],
        n_bins=10,
        linestyle=":",
        marker=".",
    )

    axes[i].set_title(f"Class {c}")
    axes[i].set_xlabel(None)
    axes[i].set_ylabel(None)
    if i == 2:
        axes[i].legend(loc="lower right")
    else:
        axes[i].get_legend().remove()

fig.suptitle("Reliability Diagrams per Class")
fig.supxlabel("Mean Predicted Probability")
fig.supylabel("Fraction of Class")
_ = fig.show()

# %%
# Interpretation
# ^^^^^^^^^^^^^^
#
# Across all three classes, we observe the following patterns:
#
# 1. **Class 0**:
#
#    - :class:`~sklearn.naive_bayes.GaussianNB` is under-confident in the
#      mid-probability region and becomes over-confident near 1.0.
#
#    - :ref:`Sigmoid calibration <sigmoid_regressor>` can overcorrect and
#      shows noticeable fluctuations.
#
#    - :ref:`Temperature scaling <temperature_scaling>` stays closer to the
#      diagonal throughout, indicating better calibration.
#
# 2. **Class 1**:
#
#    - :class:`~sklearn.naive_bayes.GaussianNB` again shows systematic
#      under-confidence in mid-ranges and over-confidence at high predicted
#      probabilities.
#
#    - :ref:`Sigmoid calibration <sigmoid_regressor>` oscillates between
#      over- and under-correction, resulting in a non-monotonic curve.
#
#    - :ref:`Temperature scaling <temperature_scaling>` produces a smoother,
#      more monotonic reliability curve that lies closest to the diagonal.
#
# 3. **Class 2**:
#
#    - :class:`~sklearn.naive_bayes.GaussianNB` underestimates the true class
#      frequency for mid-range probabilities.
#
#    - :ref:`Sigmoid calibration <sigmoid_regressor>` continues to fluctuate,
#      deviating from the diagonal in both directions.
#
#    - :ref:`Temperature scaling <temperature_scaling>` remains the most
#      consistent, providing the curve closest to on-diagonal behaviour.
#
# Overall, :ref:`temperature scaling <temperature_scaling>` yields the most
# reliable probability estimates among the three in this example.
