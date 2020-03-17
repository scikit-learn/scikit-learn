"""
==============================
Probability Calibration curves
==============================

When performing classification one often wants to predict not only the class
label, but also the associated probability. This probability gives some
kind of confidence on the prediction. This example demonstrates how to display
how well calibrated the predicted probabilities are and how to calibrate an
uncalibrated classifier.

The experiment is performed on an artificial dataset for binary classification
with 100,000 samples (1,000 of them are used for model fitting) with 20
features. Of the 20 features, only 2 are informative and 10 are redundant. The
first figure shows the estimated probabilities obtained with logistic
regression, Gaussian naive Bayes, and Gaussian naive Bayes with both isotonic
calibration and sigmoid calibration. The calibration performance is evaluated
with Brier score, reported in the legend (the smaller the better). One can
observe here that logistic regression is well calibrated while raw Gaussian
naive Bayes performs very badly. This is because of the redundant features
which violate the assumption of feature-independence and result in an overly
confident classifier, which is indicated by the typical transposed-sigmoid
curve.

Calibration of the probabilities of Gaussian naive Bayes with isotonic
regression can fix this issue as can be seen from the nearly diagonal
calibration curve. Sigmoid calibration also improves the brier score slightly,
albeit not as strongly as the non-parametric isotonic regression. This can be
attributed to the fact that we have plenty of calibration data such that the
greater flexibility of the non-parametric model can be exploited.

The second figure shows the calibration curve of a linear support-vector
classifier (LinearSVC). LinearSVC shows the opposite behavior as Gaussian
naive Bayes: the calibration curve has a sigmoid curve, which is typical for
an under-confident classifier. In the case of LinearSVC, this is caused by the
margin property of the hinge loss, which lets the model focus on hard samples
that are close to the decision boundary (the support vectors).

Both kinds of calibration can fix this issue and yield nearly identical
results. This shows that sigmoid calibration can deal with situations where
the calibration curve of the base classifier is sigmoid (e.g., for LinearSVC)
but not where it is transposed-sigmoid (e.g., Gaussian naive Bayes).
"""
print(__doc__)

# Author: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD Style.

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, f1_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC

matplotlib.style.use("classic")

# Create dataset of classification task with many redundant and few
# informative features
X, y = datasets.make_classification(
    n_samples=100000, n_features=20, n_informative=2, n_redundant=10, random_state=42
)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.99, random_state=42
)


def plot_calibration_curve(estimator, estimator_name, fig_index):
    """Plot calibration curve for est w/o and with calibration. """
    plt.figure(fig_index, figsize=(8, 11))
    # ax_cali is for plotting the calibration plot
    ax_cali = plt.subplot2grid((4, 2), (0, 0), rowspan=2, colspan=2)
    axes = []
    # axes are for plotting distributions of predicted probabilities
    for i in range(2):
        for j in range(2):
            if i == 0 and j == 0:
                _ax = plt.subplot2grid((4, 2), (i + 2, j))
                ax_ref = _ax  # reference ax for sharing X and Y axes
            else:
                _ax = plt.subplot2grid((4, 2), (i + 2, j), sharex=ax_ref, sharey=ax_ref)
            axes.append(_ax)

    classifiers = [
        # Logistic regression with no calibration as baseline,
        LogisticRegression(C=1.0),
        # the raw estimator without calibraition
        estimator,
        # Calibrated with isotonic calibration
        CalibratedClassifierCV(estimator, cv=2, method="isotonic"),
        # Calibrated with sigmoid calibration
        CalibratedClassifierCV(estimator, cv=2, method="sigmoid"),
    ]
    labels = [
        "Logistic",
        estimator_name,
        estimator_name + " + Isotonic",
        estimator_name + " + Sigmoid",
    ]
    markers = ["o", "^", "s", "d"]
    colors = ["blue", "red", "orange", "magenta"]

    ax_cali.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
    for k, (clf, label, marker, color, ax) in enumerate(
        zip(classifiers, labels, markers, colors, axes)
    ):
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        if hasattr(clf, "predict_proba"):
            prob_pos = clf.predict_proba(X_test)[:, 1]
        else:  # use decision function
            prob_pos = clf.decision_function(X_test)
            _min = prob_pos.min()
            _max = prob_pos.max()
            prob_pos = (prob_pos - _min) / (_max - _min)

        brier_score = brier_score_loss(y_test, prob_pos, pos_label=y.max())
        print(label)
        print(f"\tBrier: {brier_score:.4f}")
        print(f"\tPrecision: {precision_score(y_test, y_pred):.4f}")
        print(f"\tRecall: {recall_score(y_test, y_pred):.4f}")
        print(f"\tF1: {f1_score(y_test, y_pred):.4f}\n")

        prob_true, prob_pred = calibration_curve(y_test, prob_pos, n_bins=20)

        ax_cali.plot(
            prob_pred,
            prob_true,
            marker=marker,
            color=color,
            markeredgecolor="none",
            label=f"{label}, Brier: {brier_score:.4f}",
        )

        ax.hist(
            prob_pos,
            bins=np.arange(0, 1.01, 0.04),
            density=False,
            color=color,
            edgecolor="none",
        )

        ax.set_xlabel("Predicted P(Y=1)")
        ax.set_ylabel("Count")
        ax.set_xlim(-0.02, 1.02)
        ax.set_title(label)
        ax.grid()

    ax_cali.set_xlabel("Mean predicted value per bin")
    ax_cali.set_ylabel("Fraction of positives per bin")
    ax_cali.set_xlim(-0.02, 1.02)
    ax_cali.set_ylim([-0.05, 1.05])
    ax_cali.legend(loc="lower right", fontsize=12)
    ax_cali.set_title("Calibration plots  (reliability curve)")
    ax_cali.grid()

    plt.tight_layout()


# Plot calibration curve for Gaussian Naive Bayes
plot_calibration_curve(GaussianNB(), "Naive Bayes", 1)

# Plot calibration curve for Linear SVC
plot_calibration_curve(LinearSVC(max_iter=10000), "SVC", 2)

plt.show()
