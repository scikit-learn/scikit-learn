"""
===========================================
Efficienct GridSearchCV with use_warm_start
===========================================

A number of estimators are able to reuse a previously fit model as certain
parameters change.  This is facilitated by a ``warm_start`` parameter.  For
:class:`ensemble.GradientBoostingClassifier`, for instance, with
``warm_start=True``, fit can be called repeatedly with the same data while
increasing its ``n_estimators`` parameter.

:class:`model_selection.GridSearchCV` can efficiently search over such
warm-startable parameters through its ``use_warm_start`` parameter.  This
example compares ``GridSearchCV`` performance for searching over
``n_estimators`` in :class:`ensemble.GradientBoostingClassifier` with
and without ``use_warm_start='n_estimators'``.  """

# Authors: Vighnesh Birodkar <vighneshbirodkar@nyu.edu>
#          Raghav RV <rvraghav93@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn.ensemble import GradientBoostingClassifier
from sklearn import datasets
from sklearn.model_selection import GridSearchCV

print(__doc__)

data_list = [datasets.load_iris(return_X_y=True), datasets.make_hastie_10_2()]
names = ["Iris Data", "Hastie Data"]

search_n_estimators = range(1, 20)

times = []

for use_warm_start in [None, "n_estimators"]:
    for X, y in data_list:
        gb_gs = GridSearchCV(
            GradientBoostingClassifier(random_state=42, warm_start=True),
            param_grid={
                "n_estimators": search_n_estimators,
                "min_samples_leaf": [1, 5],
            },
            scoring="f1_micro",
            cv=3,
            refit=True,
            verbose=True,
            use_warm_start=use_warm_start,
        ).fit(X, y)
        times.append(gb_gs.cv_results_["mean_fit_time"].sum())


plt.figure(figsize=(9, 5))
bar_width = 0.2
n = len(data_list)
index = np.arange(0, n * bar_width, bar_width) * 2.5
index = index[0:n]

true_times = times[len(times) // 2 :]
false_times = times[: len(times) // 2]


plt.bar(
    index, true_times, bar_width, label='use_warm_start="n_estimators"', color="green"
)
plt.bar(
    index + bar_width, false_times, bar_width, label="use_warm_start=None", color="red"
)

plt.xticks(index + bar_width, names)

plt.legend(loc="best")
plt.grid(True)

plt.xlabel("Datasets")
plt.ylabel("Mean fit time")
plt.show()
