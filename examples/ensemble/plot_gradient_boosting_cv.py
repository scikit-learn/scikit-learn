"""
===================================
Early stopping of Gradient Boosting
===================================

Gradient boosting works by fitting a regression tree at each iteration to
improve its prediction over the training data. The improvement in prediction
succesively reduces as more regression trees are added. The
:class:`sklearn.ensemble.GradientBoostingClassifierCV` class halts the
training when there is no significant improvement in accuracy on the validation
data during the last few iterations. This example illustrates how the
:class:`sklearn.ensemble.GradientBoostingClassifierCV` class can acheive almost
the same accuracy as compared to the
:class:`sklearn.ensemble.GradientBoostingClassifier` class with significantly
lesser number of estimators. This can save memory and prediction time.
Moreover, the class can be used as a drop-in replacement for the
:class:`sklearn.ensemble.GradientBoostingClassifier` class.

"""

print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn import ensemble
from sklearn import datasets
from sklearn.model_selection import train_test_split

data_list = [datasets.load_iris(), datasets.load_digits()]
data_list = [(d.data, d.target) for d in data_list]
data_list += [datasets.make_hastie_10_2()]
names = ['Iris Data', 'Digits Data', 'Hastie Data']

n_gb = []
score_gb = []
n_gbcv = []
score_gbcv = []

for X, y in data_list:

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33,
                                                        random_state=13)

    gbcv = ensemble.GradientBoostingClassifierCV()
    gb = ensemble.GradientBoostingClassifier()

    gb.fit(X_train, y_train)
    gbcv.fit(X_train, y_train)

    score_gb.append(gb.score(X_test, y_test))
    score_gbcv.append(gbcv.best_estimator_.score(X_test, y_test))

    n_gb.append(gb.n_estimators)
    n_gbcv.append(gbcv.best_estimator_.n_estimators)

plt.figure(figsize=(9, 5))
bar_width = 0.2
n = len(data_list)
index = np.arange(0, n*bar_width, bar_width)*2.5
index = index[0:n]

plt.subplot(121)

plt.bar(index, n_gb, bar_width, label='Gradient Boosting', color='teal')
plt.bar(index + bar_width, n_gbcv, bar_width,
        label='Gradient Boosting CV', color='cyan')

max_y = np.amax(np.maximum(n_gb, n_gbcv))

plt.xticks(index + bar_width, names)
plt.yticks(np.arange(0, 130, 10))

plt.ylim([0, 130])
plt.legend(loc='best')
plt.grid(True)

plt.xlabel('Datasets')
plt.ylabel('Number of Estimators')

plt.subplot(122)

plt.bar(index, score_gb, bar_width, label='Gradient Boosting', color='crimson')
plt.bar(index + bar_width, score_gbcv, bar_width,
        label='Gradient Boosting CV', color='coral')
max_y = np.amax(np.maximum(score_gb, score_gbcv))

plt.xticks(index + bar_width, names)
plt.yticks(np.arange(0, 1.3, 0.1))

plt.ylim([0, 1.3])
plt.legend(loc='best')
plt.grid(True)

plt.xlabel('Datasets')
plt.ylabel('Score')

plt.show()
