"""
========================================
Plot multi-class SGD on the iris dataset
========================================

Plot decision surface of SGD one-against-all (OVA) classifier.
The hyperplanes corresponding to the three OVA classifiers
are represented by the dashed lines. 

"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import sgd, datasets

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2] # we only take the first two features. We could
                     # avoid this ugly slicing by using a two-dim dataset
y = iris.target

# shuffle
idx = np.arange(X.shape[0])
np.random.seed(13)
np.random.shuffle(idx)
X = X[idx]
y = y[idx]

# standardize
mean = X.mean(axis=0)
std = X.std(axis=0)
X = (X-mean) / std

h=.02 # step size in the mesh

clf = sgd.SGD(alpha=0.01, n_iter=100).fit(X, y)

# create a mesh to plot in
x_min, x_max = X[:,0].min()-1, X[:,0].max()+1
y_min, y_max = X[:,1].min()-1, X[:,1].max()+1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                     np.arange(y_min, y_max, h))

pl.set_cmap(pl.cm.Paired)

# Plot the decision boundary. For that, we will asign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.set_cmap(pl.cm.Paired)
cs = pl.contourf(xx, yy, Z)
pl.axis('tight')

# Plot also the training points
pl.scatter(X[:,0], X[:,1], c=y)
pl.title("Decision surface of multi-class SGD")
pl.axis('tight')

# Plot the three one-against-all classifiers
xmin, xmax = pl.xlim()
ymin, ymax = pl.ylim()
coef = clf.coef_
intercept = clf.intercept_

def plot_hyperplane(i):
    def line(x1):
        return (-(x1 * coef[i, 1]) - intercept[i]) / coef[i, 0]
    pl.plot([xmin, xmax], [line(xmin), line(xmax)],
            ls="--", color=pl.cm.Paired(clf.classes[i]))

for i in range(3):
    plot_hyperplane(i)

pl.show()
