"""
==================================================
Plot Label Propagation versus SVM with iris dataset
==================================================

Performance comparison between Label Propagation in the semi-supervised setting
to SVM in the supervised setting in the iris dataset. 
"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import svm, datasets

X = iris.data
Y = iris.target

lp = label_propagation.LabelPropagation()
svc = svm.SVC(kernel='linear').fit(X, Y)

# title for the plots
titles = ['Label Propagation 20% labeled data',
          'Support Vector Machine 80% labeled data']
# THE REST IS TODO !!!
