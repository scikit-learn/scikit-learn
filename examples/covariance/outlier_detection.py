"""
==========================================
Outlier detection with several methods.
==========================================

This example illustrates two ways of performing outliers detection when
the amount of contamination is known:
 - using the One-Class SVM and its ability to capture the shape of the
   data set, hence performing better when two clusters are identifiables;
 - based on a robust estimator of covariance, which is assuming that the
   data are Gaussian distributed and performs better than the One-Class SVM
   in that case.

The ground truth about inliers and outliers is given by the points colors
while the orange-filled area indicates which points are reported as outliers
by each method.

"""
print __doc__

import numpy as np
import pylab as pl
from sklearn import svm
from sklearn.covariance import EllipticData

n_samples = 200

xx, yy = np.meshgrid(np.linspace(-7, 7, 500), np.linspace(-7, 7, 500))
clusters_distances = [2, 1, 0]
for i, offset in enumerate(clusters_distances):
    # Data generation
    X1 = 0.3 * np.random.randn(int(0.45 * n_samples), 2) - offset
    X2 = 0.3 * np.random.randn(int(0.45 * n_samples), 2) + offset
    X = np.r_[X1, X2]
    # Add 10 % of outliers (leads to nu=0.1)
    n_outliers = int(0.1 * n_samples)
    X = np.r_[
        X, np.random.uniform(low=-6, high=6, size=(n_outliers, 2))]

    # Fit the model with the One-Class SVM
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
    clf.fit(X)
    # plot the levels lines and the points
    Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    y_pred = clf.predict(X)
    pl.figure(figsize=(12, 5))
    subplot = pl.subplot(1, 2, 1)
    subplot.set_title("Outlier detection")
    subplot.contourf(
        xx, yy, Z, levels=np.linspace(Z.min(), 0, 7), cmap=pl.cm.Blues_r)
    a = subplot.contour(xx, yy, Z, levels=[0], linewidths=2, colors='red')
    subplot.contourf(xx, yy, Z, levels=[0, Z.max()], colors='orange')
    b = subplot.scatter(X[:-n_outliers, 0], X[:-n_outliers, 1], c='white')
    c = subplot.scatter(X[-n_outliers:, 0], X[-n_outliers:, 1], c='black')
    subplot.axis('tight')
    subplot.legend(
        [a.collections[0], b, c],
        ['learnt decision function', 'true inliers', 'true outliers'])
    subplot.set_xlabel("1. with the One-Class SVM")
    subplot.set_xlim((-7, 7))
    subplot.set_ylim((-7, 7))

    # Fit the model with an estimator of covariance
    clf = EllipticData(contamination=.1)
    clf.fit(X)
    outliers, threshold = clf.predict(X)
    # plot the levels lines and the points
    Z = (clf.decision_function(np.c_[xx.ravel(), yy.ravel()])) ** 0.33
    Z = Z.reshape(xx.shape)
    subplot = pl.subplot(1, 2, 2)
    subplot.set_title("Outlier detection")
    subplot.contourf(xx, yy, Z,
                     levels=np.linspace(threshold ** 0.33, Z.max(), 7),
                     cmap=pl.cm.Blues)
    a = subplot.contour(
        xx, yy, Z, levels=[threshold ** 0.33], linewidths=2, colors='red')
    subplot.contourf(
        xx, yy, Z, levels=[0., threshold ** 0.33], colors='orange')
    b = subplot.scatter(X[:-n_outliers, 0], X[:-n_outliers, 1], c='white')
    c = subplot.scatter(X[-n_outliers:, 0], X[-n_outliers:, 1], c='black')
    subplot.axis('tight')
    subplot.legend(
        [a.collections[0], b, c],
        ['learnt decision function', 'true inliers', 'true outliers'])
    subplot.set_xlabel("2. from a robust covariance estimator")
    subplot.set_xlim((-7, 7))
    subplot.set_ylim((-7, 7))

pl.show()
