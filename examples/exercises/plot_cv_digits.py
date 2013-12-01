"""
=============================================
Cross-validation on Digits Dataset Exercise
=============================================

A tutorial excercise using Cross-validation with an SVM on the Digits dataset.

This exercise is used in the :ref:`cv_generators_tut` part of the
:ref:`model_selection_tut` section of the :ref:`stat_learn_tut_index`.
"""
print(__doc__)


import numpy as np
from sklearn import cross_validation, datasets, svm

digits = datasets.load_digits()
X = digits.data
y = digits.target

svc = svm.SVC(kernel='linear')
C_s = np.logspace(-10, 0, 10)

scores = list()
scores_std = list()
for C in C_s:
    svc.C = C
    this_scores = cross_validation.cross_val_score(svc, X, y, n_jobs=1)
    scores.append(np.mean(this_scores))
    scores_std.append(np.std(this_scores))

# Do the plotting
import pylab as pl
pl.figure(1, figsize=(4, 3))
pl.clf()
pl.semilogx(C_s, scores)
pl.semilogx(C_s, np.array(scores) + np.array(scores_std), 'b--')
pl.semilogx(C_s, np.array(scores) - np.array(scores_std), 'b--')
locs, labels = pl.yticks()
pl.yticks(locs, map(lambda x: "%g" % x, locs))
pl.ylabel('CV score')
pl.xlabel('Parameter C')
pl.ylim(0, 1.1)
pl.show()
