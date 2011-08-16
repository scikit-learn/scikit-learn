"""
===================================================
Plot Label Propagation versus SVM with iris dataset
===================================================

Comparison for decision boundary generated on iris dataset
between Label Propagation and SVM.

This demonstrates Label Propagation learning a good boundary
even with a small amount of labeled data.

"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import datasets
from scikits.learn import svm
from scikits.learn import label_propagation

iris = datasets.load_iris()

X = iris.data[:, :2]
Y = iris.target

# step size in the mesh
h = .02

Y_30 = np.copy(Y)
Y_30[np.random.rand(len(Y)) < 0.3] = 3
Y_50 = np.copy(Y)
Y_50[np.random.rand(len(Y)) < 0.5] = 3
# we create an instance of SVM and fit out data. We do not scale our
# data since we want to plot the support vectors
ls30 = (label_propagation.LabelSpreading().fit(X, Y_30,\
        unlabeled_identifier=3), Y_30)
ls50 = (label_propagation.LabelSpreading().fit(X, Y_50,\
        unlabeled_identifier=3), Y_50)
ls100 = (label_propagation.LabelSpreading().fit(X, Y), Y)
rbf_svc = (svm.SVC(kernel='rbf').fit(X, Y), Y)

# create a mesh to plot in
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                     np.arange(y_min, y_max, h))

# title for the plots
titles = ['Label Spreading 30% data',
          'Label Spreading 50% data',
          'Label Spreading 100% data',
          'SVC with rbf kernel']


pl.set_cmap(pl.cm.Paired)

for i, (clf, Y_train) in enumerate((ls30, ls50, ls100, rbf_svc)):
    # Plot the decision boundary. For that, we will asign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    pl.subplot(2, 2, i + 1)
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.set_cmap(pl.cm.Paired)
    pl.contourf(xx, yy, Z)
    pl.axis('off')

    # Plot also the training points
    pl.scatter(X[:, 0], X[:, 1], c=(Y_train))

    pl.title(titles[i])

pl.show()
