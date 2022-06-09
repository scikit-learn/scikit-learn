"""
=============================================
Cheap Cross Validation using Out of Bag Error
=============================================

Out-of-bag error is an estimation of the generalisation error of an estimator, much
like the cross validation error. However, it is much cheaper to calculate the out-of-bag
error, since we can do so with a single training pass of the dataset, as opposed to
:math:`k` passes in :math:`k`-fold cross validation [1]_.

In this example we perform greedy sequential feature selection with a small number of
features, and compare the runtime for out-of-bag error and cross validation error.

.. [1] D.H. Wolpert and W.G. Macready, "An Efficient Method To Estimate Bagging's
        Generalization Error.", 1999
"""
# Author: Michael Milton <michael.r.milton@gmail.com>
#
# License: BSD 3 Clause

from timeit import default_timer

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SequentialFeatureSelector

# Generate a binary classification dataset.
X, y = make_classification(
    n_samples=500,
    n_features=5,
    n_clusters_per_class=1,
    n_informative=4,
    n_redundant=1,
    random_state=0,
)


# Create functions to allow us to use the oob score for cross validation
def oob_split(X: np.ndarray):
    """
    "Split" the data into a training set with all the data, and an empty validation set
    """
    yield np.arange(X.shape[0]), np.empty(0, dtype=int)


def oob_scorer(estimator, X, y):
    """
    Score using the out-of-bag score
    """
    return estimator.oob_score_


# We use the same estimator and feature selection params for each selector
kwargs = {
    "estimator": RandomForestClassifier(oob_score=True, n_estimators=300),
    "tol": 1,
    "n_features_to_select": "auto",
}

# Define 3 cross validators to compare
selectors = {
    "OOB": SequentialFeatureSelector(
        cv=list(oob_split(X)), scoring=oob_scorer, **kwargs
    ),
    "3-fold": SequentialFeatureSelector(cv=3, **kwargs),
    "5-fold": SequentialFeatureSelector(**kwargs, cv=5),
}

# Perform feature selection, and plot the results
times = []
for name, cv in selectors.items():
    start = default_timer()
    cv.fit(X, y)
    end = default_timer()
    times.append(end - start)
plt.bar(x=list(selectors.keys()), height=times)
plt.xlabel("Cross Validator")
plt.ylabel("Grid Search Time (s)")
plt.tight_layout()
plt.show()

# Check that we achieved the same results with all cross validators
assert np.array_equal(selectors["OOB"].support_, selectors["3-fold"].support_)
assert np.array_equal(selectors["3-fold"].support_, selectors["5-fold"].support_)
