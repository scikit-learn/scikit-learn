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

import numpy as np

np.random.seed(0)

import matplotlib
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.calibration import calibration_curve

matplotlib.style.use("classic")


def init_data():
    X, y = datasets.make_classification(
        n_samples=100000, n_features=20, n_informative=2, n_redundant=2
    )

    train_samples = 100  # Samples used for training the models

    X_train = X[:train_samples]
    X_test = X[train_samples:]
    y_train = y[:train_samples]
    y_test = y[train_samples:]
    return X_train, X_test, y_train, y_test


def init_classifiers():
    """Create classifiers"""
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
    return list(zip(labels, classifiers))


def predict_prob(clf, X_test):
    if hasattr(clf, "predict_proba"):
        prob_pos = clf.predict_proba(X_test)[:, 1]
    else:  # use decision function
        prob_pos = clf.decision_function(X_test)
        _min = prob_pos.min()
        _max = prob_pos.max()
        prob_pos = (prob_pos - _min) / (_max - _min)
    return prob_pos


def init_axes():
    """prepare axes for plotting"""
    plt.figure(figsize=(8, 11))
    ax_calibration = plt.subplot2grid((4, 2), (0, 0), rowspan=2, colspan=2)
    axes_hist = []  # for distribution plots
    for i in range(2):
        for j in range(2):
            if i == 0 and j == 0:
                _ax = plt.subplot2grid((4, 2), (i + 2, j))
                ax_ref = _ax
            else:
                _ax = plt.subplot2grid((4, 2), (i + 2, j), sharex=ax_ref, sharey=ax_ref)
            axes_hist.append(_ax)
    return ax_calibration, axes_hist


def main():
    ax_cali, axes_hist = init_axes()
    ax_cali.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")

    markers = ["o", "^", "s", "d"]
    colors = ["blue", "red", "orange", "magenta"]

    X_train, X_test, y_train, y_test = init_data()
    for k, (label, clf) in enumerate(init_classifiers()):
        clf.fit(X_train, y_train)
        prob_pos = predict_prob(clf, X_test)
        prob_true, prob_pred = calibration_curve(y_test, prob_pos, n_bins=20)

        ax, marker, color = axes_hist[k], markers[k], colors[k]

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
            color=color,
            edgecolor="none",
            label=label,
        )

        ax.set_xlabel("Predicted P(Y=1)")
        ax.set_ylabel("Density")
        ax.set_xlim(-0.02, 1.02)
        ax.set_title(label)
        ax.grid()

        if k in [1, 3]:
            ax.set_ylabel("")
        if k in [0, 1]:
            ax.set_xlabel("")

    ax_cali.set_xlabel("Mean predicted value per bin")
    ax_cali.set_ylabel("Fraction of positives per bin")
    ax_cali.set_ylim([-0.05, 1.05])
    ax_cali.legend(loc="lower right", fontsize=12)
    ax_cali.set_title("Cali plots  (reliability curve)")
    ax_cali.grid()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
