"""
===================
Isotonic Regression
===================

An illustration of the isotonic regression on generated data.
"""

# Author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
# Licence: BSD

import numpy as np
import pylab as pl
from matplotlib.collections import LineCollection

from sklearn.linear_model import IsotonicRegression

n = 100
x = np.arange(n)
y = np.random.randint(-50, 50, size=(n,)) + 50. * np.log(1 + np.arange(n))

###############################################################################
# Fit IsotonicRegression model

ir = IsotonicRegression()
y_ = ir.fit_transform(x, y)

###############################################################################
# plot result

segments = [[[i, y[i]], [i, y_[i]]] for i in range(n)]
lc = LineCollection(segments, zorder=0)
lc.set_array(np.ones(len(y)))
lc.set_linewidths(0.5 * np.ones(n))

fig = pl.figure()
pl.plot(x, y, 'r.', markersize=12)
pl.plot(x, y_, 'g.-', markersize=12)
pl.gca().add_collection(lc)
pl.legend(('Data', 'Fit'), loc='lower right')
pl.title('Isotonic regression')
pl.show()
