"""
======================================================================
Decision threshold (cutoff point) calibration on breast cancer dataset
======================================================================

Machine learning classifiers often base their predictions on real-valued
decision functions that don't always have accuracy as their objective. Moreover
the learning objective of a model can differ from the user's needs hence using
an arbitrary decision threshold as defined by the model can be not ideal.

The CutoffClassifier can be used to calibrate the decision threshold of a model
in order to increase the classifier's trustworthiness. Optimization objectives
during the decision threshold calibration can be the true positive and / or
the true negative rate as well as the f beta score.

In this example the decision threshold calibration is applied on two
classifiers trained on the breast cancer dataset. The goal in the first case is
to maximize the f1 score of the classifiers whereas in the second the goal is
to maximize the true positive rate while maintaining a minimum true negative
rate.

As you can see after calibration the f1 score of the LogisticRegression
classifiers has increased slightly whereas the accuracy of the
AdaBoostClassifier classifier has stayed the same.

For the second goal as seen after calibration both classifiers achieve better
true positive rate while their respective true negative rates have decreased
slightly or remained stable.
"""

# Author: Prokopios Gryllos <prokopis.gryllos@sentiance.com>
#
# License: BSD 3 clause

from __future__ import division

import numpy as np

from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import confusion_matrix, f1_score
from sklearn.calibration import CutoffClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import load_breast_cancer
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split


print(__doc__)

# percentage of the training set that will be used for calibration
calibration_samples_percentage = 0.2

X, y = load_breast_cancer(return_X_y=True)

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.6,
                                                    random_state=42)

calibration_samples = int(len(X_train) * calibration_samples_percentage)

lr = LogisticRegression().fit(
    X_train[:-calibration_samples], y_train[:-calibration_samples])

y_pred_lr = lr.predict(X_test)
tn_lr, fp_lr, fn_lr, tp_lr = confusion_matrix(y_test, y_pred_lr).ravel()
tpr_lr = tp_lr / (tp_lr + fn_lr)
tnr_lr = tn_lr / (tn_lr + fp_lr)
f_one_lr = f1_score(y_test, y_pred_lr)

ada = AdaBoostClassifier().fit(
    X_train[:-calibration_samples], y_train[:-calibration_samples])

y_pred_ada = ada.predict(X_test)
tn_ada, fp_ada, fn_ada, tp_ada = confusion_matrix(y_test, y_pred_ada).ravel()
tpr_ada = tp_ada / (tp_ada + fn_ada)
tnr_ada = tn_ada / (tn_ada + fp_ada)
f_one_ada = f1_score(y_test, y_pred_ada)

# objective 1: we want to calibrate the decision threshold in order to achieve
# better f1 score
lr_f_beta = CutoffClassifier(
    lr, strategy='f_beta', method='predict_proba', beta=1, cv='prefit').fit(
    X_train[calibration_samples:], y_train[calibration_samples:])

y_pred_lr_f_beta = lr_f_beta.predict(X_test)
f_one_lr_f_beta = f1_score(y_test, y_pred_lr_f_beta)

ada_f_beta = CutoffClassifier(
    ada, strategy='f_beta', method='predict_proba', beta=1, cv='prefit'
).fit(X_train[calibration_samples:], y_train[calibration_samples:])

y_pred_ada_f_beta = ada_f_beta.predict(X_test)
f_one_ada_f_beta = f1_score(y_test, y_pred_ada_f_beta)

# objective 2: we want to maximize the true positive rate while the true
# negative rate is at least 0.7
lr_max_tpr = CutoffClassifier(
    lr, strategy='max_tpr', method='predict_proba', threshold=0.7, cv='prefit'
).fit(X_train[calibration_samples:], y_train[calibration_samples:])

y_pred_lr_max_tpr = lr_max_tpr.predict(X_test)
tn_lr_max_tpr, fp_lr_max_tpr, fn_lr_max_tpr, tp_lr_max_tpr = \
    confusion_matrix(y_test, y_pred_lr_max_tpr).ravel()
tpr_lr_max_tpr = tp_lr_max_tpr / (tp_lr_max_tpr + fn_lr_max_tpr)
tnr_lr_max_tpr = tn_lr_max_tpr / (tn_lr_max_tpr + fp_lr_max_tpr)

ada_max_tpr = CutoffClassifier(
    ada, strategy='max_tpr', method='predict_proba', threshold=0.7, cv='prefit'
).fit(X_train[calibration_samples:], y_train[calibration_samples:])

y_pred_ada_max_tpr = ada_max_tpr.predict(X_test)
tn_ada_max_tpr, fp_ada_max_tpr, fn_ada_max_tpr, tp_ada_max_tpr = \
    confusion_matrix(y_test, y_pred_ada_max_tpr).ravel()
tpr_ada_max_tpr = tp_ada_max_tpr / (tp_ada_max_tpr + fn_ada_max_tpr)
tnr_ada_max_tpr = tn_ada_max_tpr / (tn_ada_max_tpr + fp_ada_max_tpr)

print('Calibrated threshold')
print('Logistic Regression classifier: {}'.format(
    lr_max_tpr.decision_threshold_))
print('AdaBoost classifier: {}'.format(ada_max_tpr.decision_threshold_))
print('before calibration')
print('Logistic Regression classifier: tpr = {}, tnr = {}, f1 = {}'.format(
    tpr_lr, tnr_lr, f_one_lr))
print('AdaBoost classifier: tpr = {}, tpn = {}, f1 = {}'.format(
    tpr_ada, tnr_ada, f_one_ada))

print('true positive and true negative rates after calibration')
print('Logistic Regression classifier: tpr = {}, tnr = {}, f1 = {}'.format(
    tpr_lr_max_tpr, tnr_lr_max_tpr, f_one_lr_f_beta))
print('AdaBoost classifier: tpr = {}, tnr = {}, f1 = {}'.format(
    tpr_ada_max_tpr, tnr_ada_max_tpr, f_one_ada_f_beta))

#########
# plots #
#########
bar_width = 0.2

plt.subplot(2, 1, 1)
index = np.asarray([1, 2])
plt.bar(index, [f_one_lr, f_one_ada], bar_width, color='r',
        label='Before calibration')

plt.bar(index + bar_width, [f_one_lr_f_beta, f_one_ada_f_beta], bar_width,
        color='b', label='After calibration')

plt.xticks(index + bar_width / 2, ('f1 logistic', 'f1 adaboost'))

plt.ylabel('scores')
plt.title('f1 score')
plt.legend(bbox_to_anchor=(.5, -.2), loc='center', borderaxespad=0.)

plt.subplot(2, 1, 2)
index = np.asarray([1, 2, 3, 4])
plt.bar(index, [tpr_lr, tnr_lr, tpr_ada, tnr_ada],
        bar_width, color='r', label='Before calibration')

plt.bar(index + bar_width,
        [tpr_lr_max_tpr, tnr_lr_max_tpr, tpr_ada_max_tpr, tnr_ada_max_tpr],
        bar_width, color='b', label='After calibration')

plt.xticks(
    index + bar_width / 2,
    ('tpr logistic', 'tnr logistic', 'tpr adaboost', 'tnr adaboost'))
plt.ylabel('scores')
plt.title('true positive & true negative rate')

plt.subplots_adjust(hspace=0.6)
plt.show()
