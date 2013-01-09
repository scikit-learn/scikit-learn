"""Figure 10.2 from Elements of Statistical Learning, Ed. 2.

.. author:: Peter Prettenhofer <peter.prettenhofer@gmail.com>
"""
import numpy as np
from sklearn import datasets
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier


X, y = datasets.make_hastie_10_2(n_samples=12000, random_state=1)

X_test, y_test = X[2000:], y[2000:]
X_train, y_train = X[:2000], y[:2000]

dt_stump = DecisionTreeClassifier(max_depth=1, min_samples_leaf=1)
dt_stump.fit(X_train, y_train)
dt_stump_err = 1.0 - dt_stump.score(X_test, y_test)

dt = DecisionTreeClassifier(max_depth=9, min_samples_leaf=1)
dt.fit(X_train, y_train)
dt_err = 1.0 - dt.score(X_test, y_test)

ada = AdaBoostClassifier(base_estimator=dt_stump,
                         learning_rate=1,
                         n_estimators=400)
ada.fit(X_train, y_train)

# plot test error
import pylab as plt

fig = plt.figure(figsize=(12, 7), facecolor='w')
ax = fig.add_subplot(111)

ax.plot([1, 400], [dt_stump_err] * 2, 'k-', label='Decision Stump Error')
ax.plot([1, 400], [dt_err] * 2, 'k--', label='Decision Tree Error')

ada_err = np.zeros((ada.n_estimators,))
for i, y_pred in enumerate(ada.staged_predict(X_test)):
    ada_err[i] = (y_pred != y_test).mean()

ada_err_t = np.zeros((ada.n_estimators,))
for i, y_pred in enumerate(ada.staged_predict(X_train)):
    ada_err_t[i] = (y_pred != y_train).mean()

ax.plot(np.arange(ada.n_estimators) + 1, ada_err,
        label='AdaBoost Test Error',
        color='orange')
ax.plot(np.arange(ada.n_estimators) + 1, ada_err_t,
        label='AdaBoost Train Error',
        color='blue')
ax.plot(np.arange(ada.n_estimators) + 1, ada.weights_,
        label='AdaBoost Weights',
        color='red')
ax.plot(np.arange(ada.n_estimators) + 1, ada.errors_,
        label='AdaBoost Errors',
        color='green')

ax.set_ylim((0.0, 0.5))
ax.set_xlabel('n_estimators')

# shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0 * 0.5, box.y0, box.width * 0.85, box.height])

# Put a legend to the right of the current axis
l = ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
l.draw_frame(False)

plt.show()
