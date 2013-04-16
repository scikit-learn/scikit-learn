"""
==========================
Naive Bayes Classification
==========================

Sample usage of Naive Bayes classification.
It will plot the decision boundaries for each class.
"""
print __doc__

import numpy as np
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import naive_bayes, datasets

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2]  # keep only two features
y = iris.target

h = .02  # step size in the mesh

# c1reate color maps
cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
cmap_bold = ListedColormap(['#FF0000', '#00FF00', '#0000FF'])

for i in range(3):
    sample_weight = None
    class_prior_override = None
    if i == 0:
        subtitle = "(default settings)"
    elif i == 1:
        subtitle = "(sample weights tripling the importance of red points)"
        sample_weight = [3.] * 50 + [1.] * 100
    elif i == 2:
        class_prior_override = [0.15, 0.15, 0.7]
        subtitle = "(class prior override = " + \
            str(class_prior_override) + ")"
    # create an instance of GaussianNB and fit the data
    clf = naive_bayes.GaussianNB(class_prior_override=class_prior_override)
    clf.fit(X, y, sample_weight=sample_weight)
    # plot the decision boundary
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

    # put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.figure()
    pl.pcolormesh(xx, yy, Z, cmap=cmap_light)

    # plot also the training points
    pl.scatter(X[:, 0], X[:, 1], c=y, cmap=cmap_bold)
    pl.title("3-Class Gaussian Naive Bayes classification\n" + subtitle)
    pl.axis('tight')

pl.show()
