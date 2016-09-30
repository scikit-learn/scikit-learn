"""Demonstration of GridSearchCV with multiple metrics.

Multiple metric grid search can be done by setting the ``scoring`` parameter
to a list of metric strings or a dict mapping the scorer names to the scorer
callables.

The scores of all the scorers are available in the ``cv_results_`` dict at keys
ending in ``'_<scorer_name/metric_name>'`` (``'mean_test_precision'``,
``'rank_test_precision'``, etc...)

The ``best_estimator_`` and ``best_index_`` attributes are now
a dict mapping the scorer name to the best estimator and the best candidate
index corresponding to the best score for that metric which is stored in
``best_score_[<scorer_name/metric_name>]``.
"""

import numpy as np

from sklearn.model_selection import GridSearchCV
from sklearn.metrics import precision_score, recall_score, accuracy_score
from sklearn.metrics.scorer import make_scorer
from sklearn.datasets import make_classification
from sklearn.svm import LinearSVC

svc = LinearSVC(random_state=0)
X, y = make_classification()

scoring = {'p': make_scorer(precision_score),
           'r': make_scorer(recall_score),
           'a': make_scorer(accuracy_score)}

gs = GridSearchCV(svc,
                  param_grid={'C': np.linspace(0.1, 1)},
                  scoring=scoring, cv=5)
print("GridSearch instance:\n", gs.fit(X, y), "\n")
print("Best Estimators for all the scorers: \n", gs.best_estimator_, "\n")
print("Best scores for all the scorers: \n", gs.best_score_, "\n")
print("Best candidate index for all the scorers: \n", gs.best_index_, "\n")
print("The cv_results_ dict keys: \n", list(gs.cv_results_.keys()), "\n")
