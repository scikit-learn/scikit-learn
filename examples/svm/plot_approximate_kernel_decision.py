"""
========================================================
Plot decision surface of exact vs approximate kernel SVM
========================================================

Comparison of rbf kernel SVM decision surface against
linear SVM with approximate rbf kernel feature map.
We train on a projection of the iris dataset.

"""
print __doc__

import numpy as np
import pylab as pl
from sklearn import svm, datasets
from sklearn.kernel_approximation import RBFSampler
from sklearn.pipeline import Pipeline

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2]
Y = iris.target

h=.02 # step size in the mesh

# we create an instance of SVM and fit out data.
rbf_svc = svm.SVC(kernel='rbf', gamma=2, C=10).fit(X, Y)

feature_map = RBFSampler(gamma=2, n_components=100)
approx_kernel_svm = Pipeline([("feature_map", feature_map),
    ("svm", svm.SVC(kernel='linear', C=10))])
approx_kernel_svm.fit(X, Y)

# create a mesh to plot in
x_min, x_max = X[:,0].min()-1, X[:,0].max()+1
y_min, y_max = X[:,1].min()-1, X[:,1].max()+1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                     np.arange(y_min, y_max, h))

# title for the plots
titles = ['SVC with rbf kernel',
          'SVC (linear kernel) with rbf feature map\n n_components=100']

pl.figure(figsize=(12,5))
pl.set_cmap(pl.cm.Paired)

for i, clf in enumerate((rbf_svc, approx_kernel_svm)):
    # Plot the decision boundary. For that, we will asign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    pl.subplot(1, 2, i+1)
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.set_cmap(pl.cm.Paired)
    pl.contourf(xx, yy, Z)
    pl.axis('off')

    # Plot also the training points
    pl.scatter(X[:,0], X[:,1], c=Y)

    pl.title(titles[i])

pl.show()
