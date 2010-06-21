"""
Simple Gaussian Mixture model plotting example

TODO: use the faithful dataset
"""

import numpy as np
from scikits.learn import gmm, shortcuts
from datetime import datetime

n, m = 442, 2

np.random.seed(0)

C = np.array([[0., -0.7], [3.5, .7]]) # rotation and scale matrix

X = np.r_[np.dot(np.random.randn(n, 2), C),
          np.random.randn(n, 2) + np.array([10,10])]

clf = gmm.GMM(2, cvtype='full')
clf.fit(X)

shortcuts.plot_gmixture (X, clf.means, clf.covars, clf=clf)
