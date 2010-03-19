"""
An example showing feature selection.
"""

import numpy as np
import pylab as pl


################################################################################
# import some data to play with

# The IRIS dataset
from scikits.learn import datasets, svm
iris = datasets.load('iris')

# Some noisy data not correlated
E1, E2 = np.random.normal(size=(2, len(iris.data)))

x = iris.data
y = iris.target

################################################################################
pl.figure(1)
pl.clf()

################################################################################
# Univariate feature selection
from scikits.learn.feature_selection import univ_selection 
selector = univ_selection.UnivSelection(
                score_func=univ_selection.f_classif)

selector.fit(x, y)
scores = -np.log(selector.p_values_)
pl.plot(scores/scores.max(), label='Univariate score (p values)')

################################################################################
# Compare to the weights of an SVM
clf = svm.SVC(kernel='linear')
clf.fit(x, y)

svm_weights = (clf.support_**2).sum(axis=0)
pl.plot(svm_weights/svm_weights.max(), label='SVM weight')

pl.title("Comparing feature selection")
pl.xlabel('Feature number')
pl.yticks(())
pl.legend()
pl.show()
 
