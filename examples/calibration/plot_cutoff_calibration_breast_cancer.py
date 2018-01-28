"""
======================================================================
Decision threshold (cutoff point) calibration on breast cancer dataset
======================================================================

Most machine learning classifiers that offer probability estimation do so by
optimizing for accuracy (minimizing the classification error). The class with
the highest probability, which can also be interpreted as confidence or score,
is the predicted one. For a binary classification task that sets the decision
threshold arbitrarily to 0.5.

Depending on the classification task and the cost of error per class using an
arbitrary decision threshold of 0.5 can be elusive. Calibrating the decision
threshold to achieve better true positive rate or better true negative rate
or both can be a valid way for increasing the classifiers trustworthiness.

In this example the decision threshold calibration is applied on the breast
cancer dataset to maximize the true positive and true negative rate
respectively
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
from sklearn.datasets import load_breast_cancer
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split


print(__doc__)

# percentage of the training set that will be used for calibration
calibration_samples_percentage = 0.3

X, y = load_breast_cancer(return_X_y=True)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=0.6, random_state=42
)

calibration_samples = int(len(X_train) * calibration_samples_percentage)

clf_lr = LogisticRegression().fit(
    X_train[:-calibration_samples], y_train[:-calibration_samples]
)

clf_ada = AdaBoostClassifier().fit(
    X_train[:-calibration_samples], y_train[:-calibration_samples]
)

# we want to maximize the true positive rate while the true negative rate is at
# least 0.5
clf_lr_max_tpr = CutoffClassifier(
    clf_lr, method='max_tpr', cv=3, min_val_tnr=0.5
).fit(X_train[calibration_samples:], y_train[calibration_samples:])

clf_ada_max_tpr = CutoffClassifier(
    clf_ada, method='max_tpr', cv=3, min_val_tnr=0.5
).fit(X_train[calibration_samples:], y_train[calibration_samples:])

y_pred_lr = clf_lr.predict(X_test)
y_pred_ada = clf_ada.predict(X_test)
y_pred_lr_max_tpr = clf_lr_max_tpr.predict(X_test)
y_pred_ada_max_tpr = clf_ada_max_tpr.predict(X_test)

tn_lr, fp_lr, fn_lr, tp_lr = confusion_matrix(y_test, y_pred_lr).ravel()
tn_ada, fp_ada, fn_ada, tp_ada = confusion_matrix(y_test, y_pred_ada).ravel()

tn_lr_max_tpr, fp_lr_max_tpr, fn_lr_max_tpr, tp_lr_max_tpr = \
    confusion_matrix(y_test, y_pred_lr_max_tpr).ravel()
tn_ada_max_tpr, fp_ada_max_tpr, fn_ada_max_tpr, tp_ada_max_tpr = \
    confusion_matrix(y_test, y_pred_ada_max_tpr).ravel()

print('\n')
print('Calibrated threshold')
print('Logistic Regression classifier: {}'.format(clf_lr_max_tpr.threshold_))
print('AdaBoost classifier: {}'.format(clf_ada_max_tpr.threshold_))

print('\n')
print('true positive and true negative rates before calibration')

tpr_lr = tp_lr / (tp_lr + fn_lr)
tnr_lr = tn_lr / (tn_lr + fp_lr)

print('Logistic Regression classifier: tpr = {}, tnr = {}'.format(
    tpr_lr, tnr_lr
))

tpr_ada = tp_ada / (tp_ada + fn_ada)
tnr_ada = tn_ada / (tn_ada + fp_ada)

print('AdaBoost classifier: tpr = {}, tpn = {}'.format(tpr_ada, tnr_ada))

print('\n')
print('true positive and true negative rates after calibration')

tpr_lr_max_tpr = tp_lr_max_tpr / (tp_lr_max_tpr + fn_lr_max_tpr)
tnr_lr_max_tpr = tn_lr_max_tpr / (tn_lr_max_tpr + fp_lr_max_tpr)

print('Logistic Regression classifier: tpr = {}, tnr = {}'.format(
    tpr_lr_max_tpr, tnr_lr_max_tpr
))

tpr_ada_max_tpr = tp_ada_max_tpr / (tp_ada_max_tpr + fn_ada_max_tpr)
tnr_ada_max_tpr = tn_ada_max_tpr / (tn_ada_max_tpr + fp_ada_max_tpr)

print('AdaBoost classifier: tpr = {}, tnr = {}'.format(
    tpr_ada_max_tpr, tnr_ada_max_tpr
))

#######
# plots
#######
bar_width = 0.2
opacity = 0.35

index = np.asarray([1, 2, 3, 4])
plt.bar(index, [tpr_lr, tnr_lr, tpr_ada, tnr_ada],
        bar_width, alpha=opacity, color='b', label='Before')

plt.bar(index + bar_width,
        [tpr_lr_max_tpr, tnr_lr_max_tpr, tpr_ada_max_tpr, tnr_ada_max_tpr],
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
