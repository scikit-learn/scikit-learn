from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import time
from sklearn.preprocessing.data import CountFeaturizer
from collections import OrderedDict
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier

# with reference to auto_examples/ensemble/plot_ensemble_oob


RANDOM_STATE = 123

n_datapoints = 500 # 500
n_informative = 15 # 15
n_features = 25 # 25

# Generate a binary classification dataset.
X, y = make_classification(n_samples=n_datapoints, n_features=n_features,
                           n_clusters_per_class=1, n_informative=n_informative,
                           random_state=RANDOM_STATE)

# X_count is X with an additional feature column 'count'
cf = CountFeaturizer()
X_count = cf.fit_transform(X)

# NOTE: Setting the `warm_start` construction parameter to `True` disables
# support for parallelized ensembles but is necessary for tracking the OOB
# error trajectory during training.
ensemble_clfs = None


def init_clfs():
    global ensemble_clfs
    ensemble_clfs = [
        ("RandomForestClassifier, max_features='sqrt'",
            RandomForestClassifier(warm_start=True, oob_score=True,
                                   max_features="sqrt",
                                   random_state=RANDOM_STATE)),
        ("RandomForestClassifier, max_features='log2'",
            RandomForestClassifier(warm_start=True, max_features='log2',
                                   oob_score=True,
                                   random_state=RANDOM_STATE)),
        ("RandomForestClassifier, max_features=None",
            RandomForestClassifier(warm_start=True, max_features=None,
                                   oob_score=True,
                                   random_state=RANDOM_STATE))
    ]

classifiers_used = [1]
init_clfs()

label_list = ["RandomForestClassifier, max_features='sqrt', with_count",
              "RandomForestClassifier, max_features='log2', with_count",
              "RandomForestClassifier, max_features=None, with_count",
              "RandomForestClassifier, max_features='sqrt', no_count",
              "RandomForestClassifier, max_features='log2', no_count",
              "RandomForestClassifier, max_features=None, no_count"]

# Map a classifier name to a list of (<n_estimators>, <error rate>) pairs.
error_rate = OrderedDict((label, []) for index, label in
                         enumerate(label_list)
                         if index in classifiers_used
                         or (index - 3) in classifiers_used)

# Range of `n_estimators` values to explore.
min_estimators = 15
max_estimators = 175 # 175

index = 0
time_start = time.time()

for label, clf in ensemble_clfs:
    if index in classifiers_used:
        for i in range(min_estimators, max_estimators + 1):
            clf.set_params(n_estimators=i)
            clf.fit(X_count, y)
            # Record the OOB error for each `n_estimators=i` setting.
            oob_error = 1 - clf.oob_score_
            error_rate[label + ", with_count"].append((i, oob_error))
    index = index + 1

index = 0
print("Time taken using CountFeaturizer: ", (time.time() - time_start))
time_start = time.time()

# reinitialize the clfs for another test run with no counts
init_clfs()

for label, clf in ensemble_clfs:
    if index in classifiers_used:
        for i in range(min_estimators, max_estimators + 1):
            clf.set_params(n_estimators=i)
            clf.fit(X, y)
            # Record the OOB error for each `n_estimators=i` setting.
            oob_error = 1 - clf.oob_score_
            error_rate[label + ", no_count"].append((i, oob_error))
    index = index + 1

print("Time taken without CountFeaturizer: ", (time.time() - time_start))

# Generate the "OOB error rate" vs. "n_estimators" plot.
for label, clf_err in error_rate.items():
    xs, ys = zip(*clf_err)
    plt.plot(xs, ys, label=label)

plt.xlim(min_estimators, max_estimators)
plt.xlabel("n_estimators")
plt.ylabel("OOB error rate")
plt.legend(loc="upper right")
plt.show()
