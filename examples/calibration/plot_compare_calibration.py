"""
========================================
Comparison of Calibration of Classifiers
========================================

Well calibrated classifiers are probabilistic classifiers for which the output
of :term:`predict_proba` can be directly interpreted as a confidence level.
For instance, a well calibrated (binary) classifier should classify the samples
such that for the samples to which it gave a `predict_proba` value close to
0.8, approximately 80% actually belong to the positive class.
"""

# %%

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD Style.

import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.calibration import plot_calibration_curve

# %%
# Dataset
# -------
#
# We will use a synthetic binary classification dataset with 100,000 samples
# and 20 features. Of the 20 features, only 2 are informative, 2 are
# redundant (random combinations of the informative features) and the
# remaining 16 are 'useless' (random numbers). Of the 100,000 samples, only
# 100 will be used for model fitting.

X, y = make_classification(
  n_samples=100000, n_features=20, n_informative=2, n_redundant=2,
  random_state=42
)

train_samples = 100  # Samples used for training the models
X_train = X[:train_samples]
X_test = X[train_samples:]
y_train = y[:train_samples]
y_test = y[train_samples:]

# %%
# Calibration curves
# ------------------
#
# :class:`~sklearn.linear_model.LogisticRegression` returns well calibrated
# predictions as it directly optimizes log-loss. In contrast, the other methods
# return biased probabilities, with different biases per method:
#
# * :class:`~sklearn.naive_bayes.GaussianNaiveBayes` tends to push
#   probabilities to 0 or 1 (see histogram). This is mainly
#   because it makes the assumption that features are conditionally independent
#   given the class, which is not the case in this dataset which contains 2
#   redundant features.
#
# * :class:`~sklearn.ensemble.RandomForestClassifier` shows the opposite
#   behavior: the histograms show peaks at approx. 0.2 and 0.9 probability,
#   while probabilities close to 0 or 1 are very rare. An explanation for this
#   is given by Niculescu-Mizil and Caruana [1]_: "Methods such as bagging and
#   random forests that average predictions from a base set of models can have
#   difficulty making predictions near 0 and 1 because variance in the
#   underlying base models will bias predictions that should be near zero or
#   one away from these values. Because predictions are restricted to the
#   interval [0,1], errors caused by variance tend to be one- sided near zero
#   and one. For example, if a model should predict p = 0 for a case, the only
#   way bagging can achieve this is if all bagged trees predict zero. If we add
#   noise to the trees that bagging is averaging over, this noise will cause
#   some trees to predict values larger than 0 for this case, thus moving the
#   average prediction of the bagged ensemble away from 0. We observe this
#   effect most strongly with random forests because the base-level trees
#   trained with random forests have relatively high variance due to feature
#   subsetting." As a result, the calibration curve shows a characteristic
#   sigmoid shape, indicating that the classifier could trust its "intuition"
#   more and return probabilities closer to 0 or 1 typically.
#
# * :class:`~sklearn.svm.LinearSVC` shows an even more sigmoid curve as the
#   :class:`~sklearn.ensemble.RandomForestClassifier`, which is typical for
#   maximum-margin methods [1]_, which focus on difficult to classify samples
#   that are close to the decision boundary (the support vectors).

# Create classifiers
lr = LogisticRegression()
gnb = GaussianNB()
svc = LinearSVC(C=1.0)
rfc = RandomForestClassifier()

clf_list = [(lr, 'Logistic'),
            (gnb, 'Naive Bayes'),
            (svc, 'SVC'),
            (rfc, 'Random forest')]

fig = plt.figure(figsize=(10, 10))
ax1 = plt.subplot2grid((4, 2), (0, 0), rowspan=2, colspan=2)

viz_objects = {}
for clf, name in clf_list:
    clf.fit(X_train, y_train)
    viz = plot_calibration_curve(
        clf, X_test, y_test, n_bins=10, name=name, ax=ax1
    )
    viz_objects[name] = viz

ax1.grid()
ax1.set_title('Calibration plots')
ax1.set(xlabel="")

colors = plt.cm.get_cmap('tab10')
# Add histogram
for i, (_, name) in enumerate(clf_list):
    if i <= 1:
        ax = plt.subplot2grid((4, 2), (2, i))
    else:
        ax = plt.subplot2grid((4, 2), (3, i - 2))

    ax.hist(
        viz_objects[name].y_prob, range=(0, 1), bins=10, label=name,
        color=colors(i)
    )
    ax.set(title=name, xlabel="Mean predicted probability", ylabel="Count")

plt.tight_layout()
plt.show()

# %%
# References
# ----------
#
# .. [1] `Predicting Good Probabilities with Supervised Learning
#        <https://dl.acm.org/doi/pdf/10.1145/1102351.1102430>`_,
#        A. Niculescu-Mizil & R. Caruana, ICML 2005
