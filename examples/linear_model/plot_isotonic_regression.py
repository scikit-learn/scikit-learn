"""
===================
Isotonic Regression
===================

An illustration of the isotonic regression on generated data.
"""

# Author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# Licence: BSD

import numpy as np
import pylab as pl
from matplotlib.collections import LineCollection

from sklearn.linear_model import IsotonicRegression

a = 50.
n = 100
x = np.arange(n)
y = np.random.randint(-50, 50, size=(n,)) + a * np.log(1 + np.arange(n))

###############################################################################
# Fit IsotonicRegression model

ir = IsotonicRegression()
y_ = ir.fit_transform(x, y)

# Show case extrapolation outside of fit interval
x_extrapolated = np.arange(-5, n + 5)
y_extrapolated_ = ir.transform(x_extrapolated)

###############################################################################
# plot result

segments = [[[i, y[i]], [i, y_[i]]] for i in range(n)]
lc = LineCollection(segments, zorder=0)
lc.set_array(np.ones(len(y)))
lc.set_linewidths(0.5 * np.ones(n))

fig = pl.figure()
pl.plot(x, y, 'r.', markersize=12)
pl.plot(x, y_, 'g.-', markersize=12)
pl.plot(x_extrapolated, y_extrapolated_, 'k-')
pl.gca().add_collection(lc)
pl.legend(('Data', 'Fit', 'Fit with extrapolation'), loc='lower right')
pl.title('Isotonic regression')
pl.show()
