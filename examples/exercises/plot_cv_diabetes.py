"""
===============================================
Cross-validation on diabetes Dataset Exercise
===============================================

This exercise is used in the :ref:`cv_estimators_tut` part of the
:ref:`model_selection_tut` section of the :ref:`stat_learn_tut_index`.
"""
print __doc__

import numpy as np
import pylab as pl

from sklearn import cross_validation, datasets, linear_model

diabetes = datasets.load_diabetes()
X = diabetes.data[:150]
y = diabetes.target[:150]

lasso = linear_model.Lasso()

alphas = np.logspace(-4, -.5, 30)

scores = list()
scores_std = list()

for alpha in alphas:
    lasso.alpha = alpha
    this_scores = cross_validation.cross_val_score(lasso, X, y, n_jobs=1)
    scores.append(np.mean(this_scores))
    scores_std.append(np.std(this_scores))

pl.figure(figsize=(4, 3))
pl.semilogx(alphas, scores)
# plot error lines showing +/- std. errors of the scores
pl.semilogx(alphas, np.array(scores) + np.array(scores_std) / np.sqrt(len(X)),
    'b--')
pl.semilogx(alphas, np.array(scores) - np.array(scores_std) / np.sqrt(len(X)),
    'b--')
pl.ylabel('CV score')
pl.xlabel('alpha')
pl.axhline(np.max(scores), linestyle='--', color='.5')

##############################################################################
# Bonus: how much can you trust the selection of alpha?

# To answer this question we use the LassoCV object that sets its alpha
# parameter automatically from the data by internal cross-validation (i.e. it
# performs cross-validation on the training data it receives).
# We use external cross-validation to see how much the automatically obtained
# alphas differ across different cross-validation folds.
lasso_cv = linear_model.LassoCV(alphas=alphas)
k_fold = cross_validation.KFold(len(X), 3)

print "Answer to the bonus question: how much can you trust"
print "the selection of alpha?"
print
print "alpha parameter maximising the generalization score on different"
print "subsets of the data:"
print [lasso_cv.fit(X[train], y[train]).alpha_ for train, _ in k_fold]
print
print "Answer: Not very much since alphas for different subsets of data differ"
print "quite a lot."

pl.show()
