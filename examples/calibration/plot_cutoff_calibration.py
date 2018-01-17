"""
=======================================================================
Decision threshold (cutoff point) calibration for binary classification
=======================================================================

Most machine learning classifiers that offer probability estimation do so by
optimizing for accuracy (minimizing the classification error). The class with
the highest probability, which can also be interpreted as confidence or score,
is the predicted one. For a binary classification task that sets the decision
threshold arbitrarily to 0.5.

Depending on the classification task and the cost of error per class using an
arbitrary decision threshold of 0.5 can be elusive. Calibrating the decision
threshold to achieve better true positive rate or better true negative rate
or both can be a valid way for increasing the classifiers trustworthiness.

This example illustrates how the decision threshold calibration can be used on
a binary classification task with imbalanced classes for finding a decision
threshold for a logistic regression and an AdaBoost classifier with the goal
to improve the sum each their respective true positive and true negative rates.
"""

# Author: Prokopios Gryllos <prokopis.gryllos@sentiance.com>
#
# License: BSD 3 clause

from __future__ import division

import numpy as np

from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import confusion_matrix
from sklearn.calibration import CutoffClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import make_classification
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split


print(__doc__)


n_samples = 20000
calibration_samples = 2000

X, y = make_classification(n_samples=n_samples, n_features=30, random_state=42,
                           n_classes=2, shuffle=True, flip_y=0.17,
                           n_informative=6)

# unbalance dataset by removing 50% of the samples that belong to class 0
indexes_to_delete = np.random.choice(
    np.where(y == 0)[0], size=int((n_samples / 2) * 0.5)
)

X = np.delete(X, indexes_to_delete, axis=0)
y = np.delete(y, indexes_to_delete, axis=0)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=0.4, random_state=42
)

# we hold out a part of the training dataset to use for calibration
clf_lr = LogisticRegression().fit(
    X_train[:-calibration_samples], y_train[:-calibration_samples]
)

clf_ada = AdaBoostClassifier().fit(
    X_train[:-calibration_samples], y_train[:-calibration_samples]
)

clf_lr_roc = CutoffClassifier(clf_lr, method='roc', cv='prefit').fit(
    X_train[calibration_samples:], y_train[calibration_samples:]
)

clf_ada_roc = CutoffClassifier(clf_ada, method='roc', cv='prefit').fit(
    X_train[calibration_samples:], y_train[calibration_samples:]
)

y_pred_lr = clf_lr.predict(X_test)
y_pred_ada = clf_ada.predict(X_test)
y_pred_lr_roc = clf_lr_roc.predict(X_test)
y_pred_ada_roc = clf_ada_roc.predict(X_test)

tn_lr, fp_lr, fn_lr, tp_lr = confusion_matrix(y_test, y_pred_lr).ravel()
tn_ada, fp_ada, fn_ada, tp_ada = confusion_matrix(y_test, y_pred_ada).ravel()
tn_lr_roc, fp_lr_roc, fn_lr_roc, tp_lr_roc = \
    confusion_matrix(y_test, y_pred_lr_roc).ravel()
tn_ada_roc, fp_ada_roc, fn_ada_roc, tp_ada_roc = \
    confusion_matrix(y_test, y_pred_ada_roc).ravel()

print('\n')
print('Calibrated threshold')
print('Logistic Regression classifier: {}'.format(clf_lr_roc.threshold_))
print('AdaBoost classifier: {}'.format(clf_ada_roc.threshold_))

print('\n')
print('Sum of true positive and true negative rate before calibration')

tpr_lr = tp_lr / (tp_lr + fn_lr)
tnr_lr = tn_lr / (tn_lr + fp_lr)

print('Logistic Regression classifier: tpr + tnr = {} + {} = {}'.format(
    tpr_lr, tnr_lr, tpr_lr + tnr_lr
))

tpr_ada = tp_ada / (tp_ada + fn_ada)
tnr_ada = tn_ada / (tn_ada + fp_ada)

print('AdaBoost classifier: tpr + tnr = {} + {} = {}'.format(
    tpr_ada, tnr_ada, tpr_ada + tnr_ada
))

print('\n')
print('Sum of true positive and true negative rate after calibration')

tpr_lr_roc = tp_lr_roc / (tp_lr_roc + fn_lr_roc)
tnr_lr_roc = tn_lr_roc / (tn_lr_roc + fp_lr_roc)

print('Logistic Regression classifier: tpr + tnr = {} + {} = {}'.format(
    tpr_lr_roc, tnr_lr_roc, tpr_lr_roc + tnr_lr_roc
))

tpr_ada_roc = tp_ada_roc / (tp_ada_roc + fn_ada_roc)
tnr_ada_roc = tn_ada_roc / (tn_ada_roc + fp_ada_roc)

print('AdaBoost classifier: tpr + tnr = {} + {} = {}'.format(
    tpr_ada_roc, tnr_ada_roc, tpr_ada_roc + tnr_ada_roc
))

#######
# plots
#######
bar_width = 0.2
opacity = 0.35

index = np.asarray([1, 2, 3, 4])
plt.bar(index, [tpr_lr, tnr_lr, tpr_ada, tnr_ada],
        bar_width, alpha=opacity, color='b', label='Before')

plt.bar(index + bar_width, [tpr_lr_roc, tnr_lr_roc, tpr_ada_roc, tnr_ada_roc],
        bar_width, alpha=opacity, color='g', label='After')

plt.xticks(
    index + bar_width / 2,
    ('true positive rate logistic regression',
     'true negative rate logistic regression',
     'true positive rate adaboost',
     'true negative rate adaboost')
)
plt.ylabel('scores')
plt.title('Classifiers tpr and tnr before and after calibration')
plt.legend()
plt.show()

