"""
=========================================================
Using CountFeaturizer to featurize frequencies
=========================================================

Shows how to use CountFeaturizer to transform some categorical variables
into a frequency feature. CountFeaturizer can often be used to reduce
training time, classification time, and classification error.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import time
from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing.data import CountFeaturizer
from sklearn.preprocessing.data import OneHotEncoder
from collections import OrderedDict
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.pipeline import FeatureUnion

RANDOM_STATE = 30

n_datapoints = 1000
n_informative = 30
n_features = 30
n_redundant = 0

# Generate a binary classification dataset.
X, y = make_classification(n_samples=n_datapoints, n_features=n_features,
                           n_clusters_per_class=1, n_informative=n_informative,
                           n_redundant=n_redundant, random_state=RANDOM_STATE)

# only make these selected features "categorical"
discretized_features = [0, 1]
non_discretized_features = \
    list(set(range(n_features)) - set(discretized_features))
non_discretized_features_count = \
    list(set(range(n_features + 1)) - set(discretized_features))


def select_non_discrete(X, count=False):
    """Selects the non-discrete features."""
    if count:
        return X[:, non_discretized_features_count]
    return X[:, non_discretized_features]


def select_discrete(X):
    """Selects the discrete features."""
    return X[:, discretized_features]


def process_discrete(X):
    """Processes discrete features to make them categorical."""
    for feature in discretized_features:
        X_transform_col = (X[:, feature]).astype(int)
        col_min = min(np.amin(X_transform_col), 0)
        X[:, feature] = X_transform_col - col_min
    return X


time_start = time.time()
pipeline_cf = make_pipeline(
    FunctionTransformer(func=process_discrete),
    CountFeaturizer(inclusion=discretized_features),
    FunctionTransformer(func=lambda X: select_non_discrete(X, count=True)))
X_count = pipeline_cf.fit_transform(X, y=y)
cf_time_preprocessing = time.time() - time_start

time_start = time.time()
pipeline_ohe_nd = make_pipeline(FunctionTransformer(func=select_non_discrete))
pipeline_ohe_d = make_pipeline(
    FunctionTransformer(func=select_discrete),
    FunctionTransformer(func=process_discrete),
    OneHotEncoder())
pipeline_ohe = FeatureUnion(
    [("discrete", pipeline_ohe_d), ("nondiscrete", pipeline_ohe_nd)])
X_one_hot = pipeline_ohe.fit_transform(X, y=y).todense()
ohe_time_preprocessing = time.time() - time_start


def get_classifier():
    return RandomForestClassifier(warm_start=True, max_features="log2",
                                  oob_score=True, random_state=RANDOM_STATE)


clf = get_classifier()
labels = ["CountFeaturizer + RandomForestClassifier",
          "OneHotEncoder + RandomForestClassifier",
          "Only RandomForestClassifier"]
error_rate = OrderedDict((label, []) for label in labels)

min_estimators = (15 * n_datapoints // 500)
max_estimators = (175 * n_datapoints // 500)
time_start = time.time()

for i in range(min_estimators, max_estimators + 1):
    clf.set_params(n_estimators=i)
    clf.fit(X_count, y)
    oob_error = 1 - clf.oob_score_
    error_rate[labels[0]].append((i, oob_error))

print("Time taken on CountFeaturizer: ",
      (time.time() - time_start + cf_time_preprocessing))
clf = get_classifier()
time_start = time.time()

for i in range(min_estimators, max_estimators + 1):
    clf.set_params(n_estimators=i)
    clf.fit(X_one_hot, y)
    oob_error = 1 - clf.oob_score_
    error_rate[labels[1]].append((i, oob_error))

print("Time taken on OneHotEncoder: ",
      (time.time() - time_start + ohe_time_preprocessing))
clf = get_classifier()
time_start = time.time()

for i in range(min_estimators, max_estimators + 1):
    clf.set_params(n_estimators=i)
    clf.fit(X, y)
    oob_error = 1 - clf.oob_score_
    error_rate[labels[2]].append((i, oob_error))

print("Time taken on No Encoding: ", (time.time() - time_start))

# Generate the "OOB error rate" vs. "n_estimators" plot.
for label, clf_err in error_rate.items():
    xs, ys = zip(*clf_err)
    plt.plot(xs, ys, label=label)

plt.xlim(min_estimators, max_estimators)
plt.xlabel("n_estimators")
plt.ylabel("OOB error rate")
plt.legend(loc="upper right")
plt.show()
