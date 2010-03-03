"""
An example showing feature selection.
"""

import numpy as np
import pylab as pl


################################################################################
# import some data to play with

# The IRIS dataset
from scikits.learn.datasets.iris import load
SP, SW, PL, PW, LABELS = load()

# Some noisy data not correlated
E1, E2 = np.random.normal(size=(2, len(SP)))

x = np.c_[SP, SW, PL, PW, E1, E2]
y = LABELS

################################################################################
pl.figure(1)
pl.clf()

################################################################################
# Univariate feature selection
from scikits.learn.feature_select import univ_selection 
selector = univ_selection.UnivSelection(
                score_func=univ_selection.f_classif)

selector.fit(x, y)
scores = -np.log(selector.p_values_)
pl.plot(scores/scores.max(), label='Univariate score (p values)')

################################################################################
# Compare to the weights of an SVM
from scikits.learn.svm import SVM
svm = SVM(kernel_type='linear')
svm.fit(x, y)

svm_weights = (svm.support_**2).sum(axis=0)
pl.plot(svm_weights/svm_weights.max(), label='SVM weight')

pl.title("Comparing feature selection")
pl.xlabel('Feature number')
pl.yticks(())
pl.legend()
pl.show()
 
