"""
==============================
Probability Calibration curves
==============================

When performing classification one often wants to predict not only the class
label, but also the associated probability. This probability gives some
kind of confidence on the prediction. This example demonstrates how to
visualize how well calibrated the predicted probabilities are and how to
calibrate an uncalibrated classifier.
"""
print(__doc__)

# Author: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD Style.

import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (brier_score_loss, plot_calibration_curve,
                             CalibrationDisplay)
# from sklearn.metrics import (brier_score_loss, precision_score, recall_score,
#                              f1_score, plot_calibration_curve,
#                              CalibrationDisplay)
from sklearn.calibration import CalibratedClassifierCV

# %%
# Dataset
# -------
#
# We will use a synthetic binary classification dataset with 100,000 samples
# and 20 features. Of the 20 features, only 2 are informative, 10 are
# redundant (random combinations of the informative features) and the
# remaining 8 are 'useless' (random numbers). Of the 100,000 samples, 1,000
# will be used for model fitting.

X, y = make_classification(n_samples=100000, n_features=20, n_informative=2,
                           n_redundant=10, random_state=42)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.99,
                                                    random_state=42)

# %%
# Calibration curves
# ------------------
#
# First, we will compare:
#
# * :class:`~sklearn.linear_model.LogisticRegression` (used as baseline as it
#   returns well calibrated predictions as it directly optimizes log-loss)
# * Raw :class:`~sklearn.naive_bayes.GaussianNB`
# * :class:`~sklearn.naive_bayes.GaussianNB` calibrated with isotonic and
#   sigmoid calibration (see :ref:`User Guide <calibration>`)
#
# Calibration curves for all 4 conditions are plotted below, with the average
# predicted probability for each bin on the x-axis and the fraction of positive
# classes in each bin on the y-axis.
#
# Raw :class:`~sklearn.naive_bayes.GaussianNB` performs very badly because of
# the redundant features which violate the assumption of feature-independence
# and result in an overly confident classifier, which is indicated by the
# typical transposed-sigmoid curve. Calibration of the probabilities of
# :class:`~sklearn.naive_bayes.GaussianNB` with isotonic regression can fix
# this issue as can be seen from the nearly diagonal calibration curve. Sigmoid
# calibration also improves the Brier score slightly, albeit not as strongly as
# the non-parametric isotonic regression. This can be attributed to the fact
# that we have plenty of calibration data such that the greater flexibility of
# the non-parametric model can be exploited.

lr = LogisticRegression(C=1.)
gnb = GaussianNB()
gnb_isotonic = CalibratedClassifierCV(gnb, cv=2, method='isotonic')
gnb_sigmoid = CalibratedClassifierCV(gnb, cv=2, method='sigmoid')

clf_list = [(lr, 'Logistic'),
            (gnb, 'Naive Bayes'),
            (gnb_isotonic, 'Naive Bayes + Isotonic'),
            (gnb_sigmoid, 'Naive Bayes + Sigmoid')]
viz_objects = {}
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
for clf, name in clf_list:
    clf.fit(X_train, y_train)
    viz = plot_calibration_curve(
        clf, X_test, y_test, n_bins=10, name=name, ax=ax1
    )
    viz_objects[name] = viz

# Add histogram
for _, name in clf_list:
    ax2.hist(viz_objects[name].prob_pred, range=(0, 1), bins=10, label=name,
             histtype="step", lw=2)

ax2.legend(loc="upper center", ncol=2)
ax2.set(xlabel="Mean predicted probability",
        ylabel="Count")

# %%
# Next, we will compare:
#
# * :class:`~sklearn.linear_model.LogisticRegression` (used as baseline as it
#   returns well calibrated predictions as it directly optimizes log-loss)
# * :class:`~sklearn.svm.LinearSVC` with linear, isotonic and sigmoid
#   calibration (note as `LinearSVC` only has a :term:`decision_function`
#   method and thus always needs calibration to probability in [0,1])
#
# Linearly calibrated :class:`~sklearn.svm.LinearSVC` shows the opposite
# behavior as :class:`~sklearn.naive_bayes.GaussianNB`; the calibration
# curve has a sigmoid curve, which is typical for an under-confident
# classifier. In the case of :class:`~sklearn.svm.LinearSVC`, this is caused
# by the margin property of the hinge loss, which lets the model focus on
# difficult samples that are close to the decision boundary (the support
# vectors).
#
# Both kinds of calibration (sigmoid and isotonic) can fix this issue and
# yield nearly identical results.

lr = LogisticRegression(C=1.)
svc = LinearSVC(max_iter=10000)
svc_isotonic = CalibratedClassifierCV(svc, cv=2, method='isotonic')
svc_sigmoid = CalibratedClassifierCV(svc, cv=2, method='sigmoid')

clf_list = [(lr, 'Logistic'),
            (svc, 'SVC'),
            (svc_isotonic, 'SVC + Isotonic'),
            (svc_sigmoid, 'SVC + Sigmoid')]
viz_objects = {}
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
for clf, name in clf_list:
    clf.fit(X_train, y_train)
    if name != 'SVC':
        viz = plot_calibration_curve(
            clf, X_test, y_test, n_bins=10, name=name, ax=ax1
        )
        viz_objects[name] = viz
    else:
        # As LinearSVC has no `predict_proba` method
        y_prob = clf.decision_function(X_test)
        y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min())
        brier = brier_score_loss(y_test, y_prob)
        viz = CalibrationDisplay(
            y_test, y_prob, brier_score_value=brier, estimator_name=name
        )
        viz.plot(ax=ax1)
        viz_objects[name] = viz

# Add histogram
for _, name in clf_list:
    ax2.hist(viz_objects[name].prob_pred, range=(0, 1), bins=10, label=name,
             histtype="step", lw=2)

ax2.legend(loc="upper center", ncol=2)
ax2.set(xlabel="Mean predicted probability",
        ylabel="Count")

# %%
# Summary
# -------
#
# Parametric sigmoid calibration can deal with situations where the calibration
# curve of the base classifier is sigmoid (e.g., for
# :class:`~sklearn.svm.LinearSVC`) but not where it is transposed-sigmoid
# (e.g., :class:`~sklearn.naive_bayes.GaussianNB`). Non-parametric
# isotonic calibration can deal with both situations but may require more
# data.
