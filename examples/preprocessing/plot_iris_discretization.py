"""

================================================================
Using KBinsDiscretizer to discretize continuous features
================================================================

The example utilizes two features from the iris dataset (i.e., petal length
and petal width). First, data before and after discretization is plotted to
visualize the effect of discretization. Next, several classifiers are trained
using the original data and the transformed data. The results through cross
validation shows that discretization leads to some increase in accuracy.

"""

# Author: Hanmin Qin <qinhanmin2005@sina.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_iris
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier

print(__doc__)

# load the data and take two features
iris = load_iris()
X = iris.data
y = iris.target
X = X[:, [2, 3]]

# discretize the traning data
Xt = KBinsDiscretizer(n_bins=10, encode='ordinal').fit_transform(X)

# plot the data before and after discretization
plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.scatter(X[:, 0], X[:, 1], c=y, alpha=.3, cmap=plt.cm.RdYlBu,
            edgecolor='black')
plt.title('before discretization')
plt.xlabel('petal length (cm)')
plt.ylabel('petal width (cm)')
plt.subplot(122)
plt.scatter(Xt[:, 0], Xt[:, 1], c=y, cmap=plt.cm.RdYlBu, edgecolor='black')
plt.title('after discretization')
plt.xlabel('petal length (group)')
plt.ylabel('petal width (group)')
plt.show()

# train several classifiers and compare the results through cross validation
clf1 = DecisionTreeClassifier(random_state=0)
cv_result_1 = cross_val_score(clf1, X, y, cv=5)
print("DecisionTreeClassifier average score before discretization : {}"
      .format(np.mean(cv_result_1)))
print("DecisionTreeClassifier score std before discretization : {}"
      .format(np.std(cv_result_1)))
cv_result_2 = cross_val_score(clf1, Xt, y, cv=5)
print("DecisionTreeClassifier average score after discretization : {}"
      .format(np.mean(cv_result_2)))
print("DecisionTreeClassifier score std after discretization : {}"
      .format(np.std(cv_result_2)))

clf2 = SVC(random_state=0)
cv_result_1 = cross_val_score(clf2, X, y, cv=5)
print("SVC average score before discretization : {}"
      .format(np.mean(cv_result_1)))
print("SVC score std before discretization : {}"
      .format(np.std(cv_result_1)))
cv_result_2 = cross_val_score(clf2, Xt, y, cv=5)
print("SVC average score after discretization : {}"
      .format(np.mean(cv_result_2)))
print("SVC score std after discretization : {}"
      .format(np.std(cv_result_2)))
