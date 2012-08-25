"""
===================
Isotonic Regression
===================

An illustration of the isotonic regression on generated data.
"""

# Author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# Licence: BSD

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn.linear_model import IsotonicRegression

a = 50.
n = 100
X = np.random.randint(-50, 50, size=(n,)) + a * np.log(np.arange(n))

Y = IsotonicRegression().fit(X).X_
segments = [[[i, X[i]], [i, Y[i]]] for i in range(n)]
lc = LineCollection(
        segments,
        zorder=0)
lc.set_array(np.ones(len(X)))
lc.set_linewidths(0.5 * np.ones(n))

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(X, 'r.')
ax.plot(Y, 'g.-')
ax.add_collection(lc)
ax.legend(('Data', 'Fit'))
ax.set_title('Isotonic regression')
