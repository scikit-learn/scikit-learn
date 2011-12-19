"""
=========================================
Feature importances with forests of trees
=========================================

This example shows the use of forests of trees to evaluate the importance
of the pixels in an image classification task. The larger, the more important
the pixel.
"""
print __doc__

from sklearn.datasets import load_digits
from sklearn.ensemble import RandomForestClassifier

# Loading the digits dataset
digits = load_digits()
X = digits.images.reshape((len(digits.images), -1))
y = digits.target

# Build a forest and compute the pixel importances
forest = RandomForestClassifier(n_estimators=50)
forest.fit(X, y)
importances = forest.feature_importances()
importances = importances.reshape(digits.images[0].shape)

# Plot pixel importances
import pylab as pl
pl.matshow(importances)
pl.colorbar()
pl.title("Pixel importances with forests of trees")
pl.show()
