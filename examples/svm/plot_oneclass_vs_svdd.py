"""
=====================================
Comparison of One-class SVM with SVDD
=====================================

.. currentmodule:: sklearn.svm

This example illustrates the difference between :class:`OneClassSVM` and
:class:`SVDD`. For easier interpretation both algorithms use linear kernel,
i.e. build decision boundaries in the original space.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm


X = np.array([
    [1, 3],
    [1, 1],
    [2, -1],
    [2, 1],
    [3, -2],
    [3, 4],
    [3, 5],
    [4, 3],
    [4, 7],
    [5, 5],
    [5, 9],
    [6, 2],
    [7, 4],
    [7, 6],
    [4, 8],
    [9, 4],
    [9, 6],
    [12, 10],
    [14, 2],
    [11, 2]
])


oc_svm = svm.OneClassSVM(kernel='linear', nu=0.3)
oc_svm.fit(X)

w = oc_svm.coef_[0]
a = -w[0] / w[1]
b = -oc_svm.intercept_ / w[1]
x = np.linspace(-1, 15)
y = a * x + b

svdd = svm.SVDD(kernel='linear', C=0.1)
svdd.fit(X)

xx, yy = np.meshgrid(np.linspace(-1, 15), np.linspace(-3, 11))
X_test = np.c_[xx.ravel(), yy.ravel()]
Z = svdd.decision_function(X_test).reshape(xx.shape)

plt_1 = plt.scatter(X[:, 0], X[:, 1], c='black')
plt_2 = plt.plot(x, y, linewidth=2)
plt_3 = plt.contour(xx, yy, Z, levels=[0], linewidths=2, colors='red')

plt.axhline(0, color='black')
plt.axvline(0, color='black')

plt.ylim(np.min(X[:, 1]) - 1, np.max(X[:, 1]) + 1)
plt.title("Comparison of One-class SVM and SVDD")
plt.legend([plt_1, plt_2[0], plt_3.collections[0]],
           ['Training data', 'OCSVM boundary', 'SVDD boundary'])
plt.show()
