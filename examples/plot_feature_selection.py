"""
===============================
Univariate Feature Selection
===============================

An example showing univariate feature selection.

Noisy (non informative) features are added to the iris data and
univariate feature selection is applied. For each feature, we plot the
p-values for the univariate feature selection and the corresponding
weights of an SVM. We can see that univariate feature selection
selects the informative features and that these have larger SVM weights.

In the total set of features, only the 4 first ones are significant. We
can see that they have the highest score with univariate feature
selection. The SVM attributes small weights to these features, but these
weight are non zero. Applying univariate feature selection before the SVM
increases the SVM weight attributed to the significant features, and will
thus improve classification.
"""
print __doc__

import numpy as np
import pylab as pl

################################################################################
# import some data to play with

# The IRIS dataset
from scikits.learn import datasets, svm
iris = datasets.load_iris()

# Some noisy data not correlated
E = np.random.normal(size=(len(iris.data), 35))

# Add the noisy data to the informative features
x = np.hstack((iris.data, E))
y = iris.target

################################################################################
pl.figure(1)
pl.clf()

x_indices = np.arange(x.shape[-1])

################################################################################
# Univariate feature selection
from scikits.learn.feature_selection import SelectFpr, f_classif
# As a scoring function, we use a F test for classification
# We use the default selection function: the 10% most significant
# features

selector = SelectFpr(f_classif, alpha=0.1)
selector.fit(x, y)
scores = -np.log10(selector._pvalues)
scores /= scores.max()
pl.bar(x_indices-.45, scores, width=.3,
        label=r'Univariate score ($-Log(p_{value})$)',
        color='g')

################################################################################
# Compare to the weights of an SVM
clf = svm.SVC(kernel='linear')
clf.fit(x, y)

svm_weights = (clf.coef_**2).sum(axis=0)
svm_weights /= svm_weights.max()
pl.bar(x_indices-.15, svm_weights, width=.3, label='SVM weight',
        color='r')


# ################################################################################
# # Now fit an SVM with added feature selection
# selector = univ_selection.Univ(
#                 score_func=univ_selection.f_classif)

# selector.fit(x, clf.predict(x))
# svm_weights = (clf.support_**2).sum(axis=0)
# svm_weights /= svm_weights.max()
# full_svm_weights = np.zeros(selector.support_.shape)
# full_svm_weights[selector.support_] = svm_weights
# pl.bar(x_indices+.15, full_svm_weights, width=.3,
#         label='SVM weight after univariate selection',
#         color='b')

pl.title("Comparing feature selection")
pl.xlabel('Feature number')
pl.yticks(())
pl.axis('tight')
pl.legend(loc='upper right')
pl.show()

