"""
============================================================
Plot decision surfaces of stacked models on the iris dataset
============================================================

This plot compares the decision surfaces for different individual classifiers
and then combined using a Linear SVM. To make it easier to visualize, only
sepal width and sepal length features were used.

The cross_val_score for each classifier was, in decreasing order::

    stack_estimators(base_models, LinearSVC()) # 0.820
    RandomForestClassifier(random_state=1) # 0.738
    BaggingClassifier(KNeighborsClassifier(), random_state=7, # 0.740
                      max_samples=0.5, max_features=0.5)
    SVC(kernel='rbf') # 0.799

"""

print(__doc__)

import matplotlib.pyplot as plt
import numpy as np

from sklearn.ensemble import stack_estimators
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn.ensemble import BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier

from sklearn import datasets

iris = datasets.load_iris()
X = iris.data[:, :2]
y = iris.target

base_models = [SVC(kernel='rbf'),
               RandomForestClassifier(random_state=1),
               BaggingClassifier(KNeighborsClassifier(), random_state=7,
                                 max_samples=0.5, max_features=0.5)]

eclf = stack_estimators(base_models, LinearSVC())

eclf.fit(X, y)

num_plots = len(base_models) + 1

titles = ["SVC (RBF kernel)", "Random Forest", "Bag of KNNs"]

xmin, xmax = X[:, 0].min() - .5, X[:, 0].max() + .5
ymin, ymax = X[:, 1].min() - .5, X[:, 1].max() + .5

xx, yy = np.meshgrid(np.arange(xmin, xmax, .02),
                     np.arange(ymin, ymax, .02))

plt.figure(1, figsize=(9, 4))

for i, m in enumerate(base_models):
    plt.subplot(1, num_plots, i + 1)
    plt.title(titles[i])
    z = m.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
    plt.contourf(xx, yy, z, cmap=plt.cm.Paired, alpha=.5)

plt.subplot(1, num_plots, num_plots)
plt.title("Stack")
z = eclf.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
plt.contourf(xx, yy, z, cmap=plt.cm.Paired, alpha=.5)

for i in range(num_plots):
    plt.subplot(1, num_plots, i + 1)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel('Sepal length')
    plt.ylabel('Sepal width')
    plt.xticks(())
    plt.yticks(())

    plt.scatter(X[:, 0], X[:, 1], c=y, cmap=plt.cm.Paired)

plt.suptitle("Classifiers on feature subsets of the Iris dataset")
plt.show()
