"""
Sample usage of Support Vector Machines to classify a sample.
It will plot the decision surface and the support vectors.
"""
import numpy as np
import pylab as pl
from scikits.learn.svm import SVM

# import some data to play with
from scikits.learn.datasets.iris import load
SP, SW, PL, PW, LABELS = load()
X = np.c_[SP, SW]
Y = LABELS

h=.05
kernel_type='linear'

# we create an instance of SVM and fit out data
clf = SVM(kernel_type='linear')
clf.fit(X, Y)

# Plot the decision boundary. For that, we will asign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = SP.min()-1, SP.max()+1
y_min, y_max = SW.min()-1, SW.max()+1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.pcolormesh(xx, yy, Z)

# Plot also the training points
pl.scatter(SP, SW, c=Y)
# and the support vectors
pl.scatter(clf.support_[:,0], clf.support_[:, 1], marker='+')
pl.axis('tight')
pl.show()
