"""
========================================
Comparison of Calibration of Classifiers
========================================

Well calibrated classifiers are probabilistic classifiers for which the output
of the predict_proba method can be directly interpreted as a confidence level.
For instance a well calibrated (binary) classifier should classify the samples
such that among the samples to which it gave a predict_proba value close to
0.8, approx. 80% actually belong to the positive class.

LogisticRegression returns well calibrated predictions as it directly
optimizes log-loss. In contrast, the other methods return biased probabilities,
with different biases per method:

* GaussianNaiveBayes tends to push probabilities to 0 or 1 (note the counts in
  the histograms). This is mainly because it makes the assumption that features
  are conditionally independent given the class, which is not the case in this
  dataset which contains 2 redundant features.

* RandomForestClassifier shows the opposite behavior: the histograms show peaks
  at approx. 0.2 and 0.9 probability, while probabilities close to 0 or 1 are
  very rare. An explanation for this is given by Niculescu-Mizil and Caruana
  [1]_: "Methods such as bagging and random forests that average predictions
  from a base set of models can have difficulty making predictions near 0 and 1
  because variance in the underlying base models will bias predictions that
  should be near zero or one away from these values. Because predictions are
  restricted to the interval [0,1], errors caused by variance tend to be one-
  sided near zero and one. For example, if a model should predict p = 0 for a
  case, the only way bagging can achieve this is if all bagged trees predict
  zero. If we add noise to the trees that bagging is averaging over, this noise
  will cause some trees to predict values larger than 0 for this case, thus
  moving the average prediction of the bagged ensemble away from 0. We observe
  this effect most strongly with random forests because the base-level trees
  trained with random forests have relatively high variance due to feature
  subsetting." As a result, the calibration curve shows a characteristic
  sigmoid shape, indicating that the classifier could trust its "intuition"
  more and return probabilities closer to 0 or 1 typically.

* Support Vector Classification (SVC) shows an even more sigmoid curve as
  the  RandomForestClassifier, which is typical for maximum-margin methods
  (compare Niculescu-Mizil and Caruana [1]_), which focus on hard samples
  that are close to the decision boundary (the support vectors).

.. topic:: References:

    .. [1] Predicting Good Probabilities with Supervised Learning,
          A. Niculescu-Mizil & R. Caruana, ICML 2005

"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD Style.

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets
from sklearn.calibration import calibration_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC

np.random.seed(0)

matplotlib.style.use("classic")

X, y = datasets.make_classification(
    n_samples=100000, n_features=20, n_informative=2, n_redundant=2
)

train_samples = 100  # Samples used for training the models

X_train = X[:train_samples]
X_test = X[train_samples:]
y_train = y[:train_samples]
y_test = y[train_samples:]

# Create classifiers

classifiers = [
    LogisticRegression(),
    GaussianNB(),
    LinearSVC(C=1.0),
    RandomForestClassifier(),
]
labels = [
    "Logistic",
    "Naive Bayes",
    "Support Vector Classification",
    "Random Forest",
]
markers = ["o", "^", "s", "d"]
colors = ["blue", "red", "orange", "magenta"]


# #############################################################################
# Plot calibration plots

plt.figure(figsize=(8, 11))
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
            _ax = plt.subplot2grid(
                (4, 2), (i + 2, j), sharex=ax_ref, sharey=ax_ref
            )
        axes.append(_ax)


ax_cali.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
for k, (clf, label, marker, color, ax) in enumerate(
    zip(classifiers, labels, markers, colors, axes)
):
    clf.fit(X_train, y_train)
    if hasattr(clf, "predict_proba"):
        prob_pos = clf.predict_proba(X_test)[:, 1]
    else:  # use decision function
        prob_pos = clf.decision_function(X_test)
        _min = prob_pos.min()
        _max = prob_pos.max()
        prob_pos = (prob_pos - _min) / (_max - _min)
    prob_true, prob_pred = calibration_curve(y_test, prob_pos, n_bins=20)

    ax_cali.plot(
        prob_pred,
        prob_true,
        marker=marker,
        color=color,
        markeredgecolor="none",
        label=label,
    )

    ax.hist(
        prob_pos,
        bins=np.arange(0, 1.01, 0.04),
        density=False,
        color=color,
        edgecolor="none",
        label=label,
    )

    ax.set(
        title=label,
        xlabel="Predicted P(Y=1)",
        ylabel="Count",
        xlim=(-0.02, 1.02),
    )
    ax.grid()

ax_cali.set(
    title="Calibration plots  (reliability curve)",
    xlabel="Mean predicted value per bin",
    ylabel="Fraction of positives per bin",
    xlim=(-0.02, 1.02),
    ylim=(-0.05, 1.05),
)
ax_cali.legend(loc="lower right", fontsize=12)
ax_cali.grid()

plt.tight_layout()

plt.show()
