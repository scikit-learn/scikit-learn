"""
==============================================================
IsolationForest features selection with feature weight example
==============================================================

An example of features selection using feature weighting functionality of
:class:`sklearn.ensemble.IsolationForest` for anomaly detection.

This example follows closely with that in the original isolation forest paper
[1]_.

.. [1] Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. "Isolation-based
    anomaly detection." ACM Transactions on Knowledge Discovery from Data
    (TKDD) 6.1 (2012): 3.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kurtosis
from sklearn.datasets import fetch_kddcup99
from sklearn.ensemble import IsolationForest
from sklearn.metrics import roc_auc_score

rng = np.random.RandomState(2)

data = fetch_kddcup99("SF")

# Selecting numerical features
X = data.data[:, [0, 2, 3]]
# Take classes to be normal and others.
y = data.target == b'normal.'

# fit the model
clf_wo_noise = IsolationForest(random_state=rng)
clf_wo_noise.fit(X)
roc_auc_wo_noise = roc_auc_score(y, clf_wo_noise.score_samples(X))

# Add noise to data
# Get range of X
min_val = float("inf")
max_val = float("-inf")
for i in range(3):
    min_val = min(min_val, X[:, i].min())
    max_val = max(max_val, X[:, i].max())
# Adding 10 dimensions of uniform distributed noise
# We only need to cover the range of X because isolated forest
# doesn't do any handling beyond the min/max values
n_noise = 10
X_w_noise = np.hstack(
    [X, rng.uniform(min_val, max_val, size=(X.shape[0], n_noise))]
)
# Calculate feature weight using kurtosis
# Unset fisher so that normal distribution has a value of 3
# and uniform distribution has a value of 1.2, a smaller value
feature_weight = kurtosis(X_w_noise, fisher=False)
feature_weight /= feature_weight.sum()

roc_vals = []
max_features_range = range(1, 14)
for max_features in max_features_range:
    # fit the model
    clf = IsolationForest(
        random_state=rng,
        feature_weight=feature_weight,
        max_features=max_features,
    )
    clf.fit(X_w_noise)
    roc_vals.append(roc_auc_score(y, clf.score_samples(X_w_noise)))

# Plot
plt.axhline(
    roc_auc_wo_noise,
    ls='dotted',
    label="X data without noise (only has 3 features)",
)
plt.plot(
    max_features_range, roc_vals, label="X data with 10 noise features"
)
plt.axvline(
    max_features_range[np.array(roc_vals).argmax()],
    ls='dashed',
    label="Best number of features selected with feature weight",
)
plt.xlabel("Number of features selected")
plt.ylabel("AUROC")
plt.title("Isolation forest feature selection with feature weight")
plt.xlim(1, 13)
plt.ylim(0.5, 1)
plt.legend()
plt.show()
