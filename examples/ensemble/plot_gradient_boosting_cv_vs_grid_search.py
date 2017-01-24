"""
===============================
Gradient Boosting Classifier CV
===============================

Gradient boosting is an ensembling technique where several weak learners
(regression trees) are combined to yield a powerful single model, in an
iterative fashion.

:class:`sklearn.ensemble.GradientBoostingClassifierCV` enables us to
efficiently search for the best number of boosting stages. This example is a
comparison between the grid search and
:class:`sklearn.ensemble.GradientBoostingClassifierCV`.
"""

# Authors: Raghav RV <rvraghav93@gmail.com>
#          Vighnesh Birodkar <vighneshbirodkar@nyu.edu>
# License: BSD 3 clause

import time

import numpy as np
import matplotlib.pyplot as plt

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingClassifierCV
from sklearn import datasets
from sklearn.model_selection import GridSearchCV

print(__doc__)

data_list = [datasets.load_iris(), datasets.load_digits()]
data_list = [(d.data, d.target) for d in data_list]
data_list += [datasets.make_hastie_10_2()]
names = ['Iris Data', 'Digits Data', 'Hastie Data']

search_n_estimators = range(1, 20)

gbcv_times = []
gb_gs_times = []

for X, y in data_list:
    start = time.time()
    for _ in range(3):
        gb_gs = GridSearchCV(
            GradientBoostingClassifier(random_state=42),
            param_grid={'n_estimators': search_n_estimators},
            scoring='f1_micro', cv=3, refit=True).fit(X, y)
    gb_gs_times.append((time.time() - start) / 3.)

    start = time.time()
    for _ in range(3):
        gbcv = GradientBoostingClassifierCV(
            cv_n_estimators=search_n_estimators, scoring='f1_micro', cv=3,
            random_state=42).fit(X, y)
    gbcv_times.append((time.time() - start) / 3.)

plt.figure(figsize=(9, 5))
bar_width = 0.2
n = len(data_list)
index = np.arange(0, n*bar_width, bar_width)*2.5
index = index[0:n]


plt.bar(index, gbcv_times, bar_width, label='GradientBoostingClassifierCV',
        color='green')
plt.bar(index + bar_width, gb_gs_times, bar_width,
        label='GridSearchCV', color='red')

plt.xticks(index + bar_width, names)

plt.legend(loc='best')
plt.grid(True)

plt.xlabel('Datasets')
plt.ylabel('Mean fit time')
plt.show()
