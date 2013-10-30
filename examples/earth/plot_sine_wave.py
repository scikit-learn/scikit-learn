'''
=====================================
Fitting an EarthRegressor model to a sine wave
=====================================


In this example, a simple sine model is used to generate an artificial data set.  An :class:`EarthRegressor` model
is then fitted to that data set and the resulting predictions are plotted against the original data.

'''
print(__doc__)

import numpy as np
import pylab as pl
from sklearn.earth import EarthRegressor

# Create some fake data
np.random.seed(2)
m = 10000
n = 10
X = 80 * np.random.uniform(size=(m, n)) - 40
y = 100 * \
    np.abs(np.sin((X[:, 6]) / 10) - 4.0) + \
    20 * np.random.normal(size=m)

# Fit an EarthRegressor model
model = EarthRegressor(max_degree=3, minspan_alpha=.5)
model.fit(X, y)

# Print the model
print(model.trace())
print(model.summary())

# Plot the model
pl.figure()
y_hat = model.predict(X)
pl.plot(X[:, 6], y, 'r.')
pl.plot(X[:, 6], y_hat, 'b.')
pl.show()
