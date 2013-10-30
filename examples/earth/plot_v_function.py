'''
======================================================
Fitting an EarthRegressor model to a v-shaped function
======================================================


In this example, a simple piecewise linear model is used to generate an artificial data set.  An :class:`Earth` model
is then fitted to that data set and the resulting predictions are plotted against the original data.

'''
print(__doc__)

import numpy as np
from sklearn.earth import EarthRegressor
import pylab as pl

# Create some fake data
np.random.seed(2)
m = 1000
n = 10
X = 80 * np.random.uniform(size=(m, n)) - 40
y = np.abs(X[:, 6] - 4.0) + 5 * np.random.normal(size=m)

# Fit an EarthRegressor model
model = EarthRegressor(max_degree=1)
model.fit(X, y)

# Print the model
print(model.trace())
print(model.summary())

# Plot the model
y_hat = model.predict(X)
pl.figure()
pl.plot(X[:, 6], y, 'r.')
pl.plot(X[:, 6], y_hat, 'b.')
pl.show()
