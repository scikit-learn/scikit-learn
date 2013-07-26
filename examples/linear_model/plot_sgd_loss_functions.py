"""
==========================
SGD: Convex Loss Functions
==========================

An example that compares various convex loss functions.


All of the above loss functions are supported by
:class:`sklearn.linear_model.stochastic_gradient` .
"""
print(__doc__)

import numpy as np
import pylab as pl
from sklearn.linear_model.sgd_fast import SquaredHinge
from sklearn.linear_model.sgd_fast import Hinge
from sklearn.linear_model.sgd_fast import ModifiedHuber
from sklearn.linear_model.sgd_fast import SquaredLoss

###############################################################################
# Define loss functions
xmin, xmax = -4, 4
hinge = Hinge(1)
squared_hinge = SquaredHinge()
perceptron = Hinge(0)
log_loss = lambda z, p: np.log2(1.0 + np.exp(-z))
modified_huber = ModifiedHuber()
squared_loss = SquaredLoss()


###############################################################################
# Plot loss funcitons
xx = np.linspace(xmin, xmax, 100)
pl.plot([xmin, 0, 0, xmax], [1, 1, 0, 0], 'k-',
        label="Zero-one loss")
pl.plot(xx, [hinge.loss(x, 1) for x in xx], 'g-',
        label="Hinge loss")
pl.plot(xx, [perceptron.loss(x, 1) for x in xx], 'm-',
        label="Perceptron loss")
pl.plot(xx, [log_loss(x, 1) for x in xx], 'r-',
        label="Log loss")
#pl.plot(xx, [2 * squared_loss.loss(x, 1) for x in xx], 'c-',
#        label="Squared loss")
pl.plot(xx, [squared_hinge.loss(x, 1) for x in xx], 'b-',
        label="Squared hinge loss")
pl.plot(xx, [modified_huber.loss(x, 1) for x in xx], 'y--',
        label="Modified huber loss")
pl.ylim((0, 8))
pl.legend(loc="upper right")
pl.xlabel(r"$y \cdot f(x)$")
pl.ylabel("$L(y, f(x))$")
pl.show()
