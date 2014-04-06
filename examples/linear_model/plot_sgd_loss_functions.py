"""
==========================
SGD: convex loss functions
==========================

A plot that compares the various convex loss functions supported by
:class:`sklearn.linear_model.SGDClassifier` .
"""
print(__doc__)

import numpy as np
import pylab as pl


def modified_huber_loss(y_true, y_pred):
    z = y_pred * y_true
    loss = -4 * z
    loss[z >= -1] = (1 - z[z >= -1]) ** 2
    loss[z >= 1.] = 0
    return loss


xmin, xmax = -4, 4
xx = np.linspace(xmin, xmax, 100)
pl.plot([xmin, 0, 0, xmax], [1, 1, 0, 0], 'k-',
        label="Zero-one loss")
pl.plot(xx, np.where(xx < 1, 1 - xx, 0), 'g-',
        label="Hinge loss")
pl.plot(xx, -np.minimum(xx, 0), 'm-',
        label="Perceptron loss")
pl.plot(xx, np.log2(1 + np.exp(-xx)), 'r-',
        label="Log loss")
pl.plot(xx, np.where(xx < 1, 1 - xx, 0) ** 2, 'b-',
        label="Squared hinge loss")
pl.plot(xx, modified_huber_loss(xx, 1), 'y--',
        label="Modified Huber loss")
pl.ylim((0, 8))
pl.legend(loc="upper right")
pl.xlabel(r"Decision function $f(x)$")
pl.ylabel("$L(y, f(x))$")
pl.show()
