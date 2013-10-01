'''
=============================================
Fitting an Earth model to a v-shaped function
=============================================


In this example, a simple piecewise linear model is used to generate an artificial data set.  An :class:`Earth` model
is then fitted to that data set and the resulting predictions are plotted against the original data.

'''
from __future__ import print_function
import numpy
from sklearn.earth import Earth
from matplotlib import pyplot

print(__doc__)

# Create some fake data
numpy.random.seed(2)
m = 1000
n = 10
X = 80 * numpy.random.uniform(size=(m, n)) - 40
y = numpy.abs(X[:, 6] - 4.0) + 5 * numpy.random.normal(size=m)

# Fit an Earth model
model = Earth(max_degree=1)
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
