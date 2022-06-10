"""
=============================================
Cheap Cross Validation using Out of Bag Error
=============================================

Out-of-bag error is an estimation of the generalisation error of an estimator, much
like the cross validation error. However, it is much cheaper to calculate the out-of-bag
error, since we can do so with a single training pass of the dataset, as opposed to
:math:`k` passes in :math:`k`-fold cross validation [1]_.

In this example we perform grid search to select the optimal `max_depth` for the trees
in the forest. In doing so, we compare the runtime and the predicted error for the
out-of-bag error and 3-fold cross validation error.

Note however that while this approach is fairly sensible for tuning parameters such as
the tree depth, or parameter selection, it does not make sense to tune ``n_estimators``,
the number of trees in the forest in this way. OOB error also tends to overestimate
the true error, particularly when there are few observations or a large number of
predictors [2]_.

.. [1] D.H. Wolpert and W.G. Macready, "An Efficient Method To Estimate Bagging's
        Generalization Error.", 1999
.. [2] S. Janitza and R. Hornung, “On the overestimation of random forest’s out-of-bag
        error”, 2018
"""
# Author: Michael Milton <michael.r.milton@gmail.com>
#
# License: BSD 3 Clause

from timeit import default_timer

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV

# Generate a binary classification dataset.
X, y = make_classification(
    n_samples=200,
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


# We use the same estimator and parameter grid for each cross validator
kwargs = {
    "estimator": RandomForestClassifier(oob_score=True, n_estimators=100),
    "param_grid": {"max_depth": [1, 2, 3, 4]},
}

# Define the cross validators to compare
cvs = {
    "OOB": GridSearchCV(cv=list(oob_split(X)), scoring=oob_scorer, **kwargs),
    "3-Fold": GridSearchCV(cv=3, **kwargs),
}

# Perform feature selection, and plot the results
fig, (ax1, ax2) = plt.subplots(2)
times = []
colors = []
for name, cv in cvs.items():
    start = default_timer()
    cv.fit(X, y)
    end = default_timer()
    line = ax1.plot(cv.param_grid["max_depth"], cv.cv_results_["mean_test_score"])
    times.append(end - start)
    colors.append(line[0].get_color())
ax1.set_xlabel("Maximum Tree Depth")
ax1.set_ylabel("Cross Validation Accuracy")
ax2.bar(x=list(cvs.keys()), height=times, color=colors)
ax2.set_xlabel("Cross Validator")
ax2.set_ylabel("Grid Search Time (s)")
plt.tight_layout()
plt.show()
