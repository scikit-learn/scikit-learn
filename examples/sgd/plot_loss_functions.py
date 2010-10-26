"""
==========================
SGD: Convex Loss Functions
==========================

Plot the convex loss functions supported by `scikits.learn.sgd`. 
"""

import numpy as np
import pylab as pl
from scikits.learn import svm
from scikits.learn.sgd.sparse.sgd_fast_sparse import Hinge, \
     Log, ModifiedHuber

xmin, xmax = -3, 3
hinge = Hinge()
log_loss = lambda z, p: np.log2(1.0 + np.exp(-z))
modified_huber = ModifiedHuber()
xx = np.linspace(xmin, xmax, 100)
pl.plot([xmin, 0, 0, xmax], [1, 1, 0, 0], 'k-', label="Zero-one loss")
pl.plot(xx, [hinge.loss(x,1) for x in xx], 'g-', label="Hinge loss")
pl.plot(xx, [log_loss(x,1) for x in xx], 'r-', label="Log loss")
pl.plot(xx, [modified_huber.loss(x,1) for x in xx], 'y-', label="Modified huber loss")
pl.ylim((0, 5))
pl.legend(loc="upper right")
pl.xlabel(r"$y \cdot f(x)$")
pl.ylabel("$L(y, f(x))$")
#pl.axis('tight')
pl.show()

