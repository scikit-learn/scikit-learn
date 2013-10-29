'''
=====================================
Fitting an EarthRegressor model to a sine wave
=====================================


In this example, a simple sine model is used to generate an artificial data set.  An :class:`EarthRegressor` model
is then fitted to that data set and the resulting predictions are plotted against the original data.

'''
from __future__ import print_function
import numpy
from sklearn.earth import EarthRegressor
from matplotlib import pyplot

print(__doc__)

# Create some fake data
numpy.random.seed(2)
m = 10000
n = 10
X = 80 * numpy.random.uniform(size=(m, n)) - 40
y = 100 * \
    numpy.abs(numpy.sin((X[:, 6]) / 10) - 4.0) + \
    20 * numpy.random.normal(size=m)

# Fit an EarthRegressor model
model = EarthRegressor(max_degree=3, minspan_alpha=.5)
model.fit(X, y)

# Print the model
print(model.trace())
print(model.summary())

# Plot the model
y_hat = model.predict(X)
pyplot.figure()
pyplot.plot(X[:, 6], y, 'r.')
pyplot.plot(X[:, 6], y_hat, 'b.')
pyplot.show()
