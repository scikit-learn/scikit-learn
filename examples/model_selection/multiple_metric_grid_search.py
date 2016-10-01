"""Demonstration of multimetric evaluation on cross_val_score and GridSearchCV

Multiple metric grid search (or random search) can be done by setting the
``scoring`` parameter to a list of metric strings or a dict mapping the scorer
names to the scorer callables.

The scores of all the scorers are available in the ``cv_results_`` dict at keys
ending in ``'_<scorer_name/metric_name>'`` (``'mean_test_precision'``,
``'rank_test_precision'``, etc...)

The ``best_estimator_`` and ``best_index_`` attributes are now
a dict mapping the scorer name to the best estimator and the best candidate
index corresponding to the best score for that metric which is stored in
``best_score_[<scorer_name/metric_name>]``.

Similarly the :func:`sklearn.model_selection.cross_val_score`,
:func:`sklearn.model_selection.learning_curve` and
:func:`sklearn.model_selection.validation_curve` all support multiple metric
evaluation when the scoring parameter is a list/tuple of strings denoting
predefined metrics or a dict mapping the scorer names to the scorer callables.
"""
import numpy as np

from sklearn.datasets import make_classification
from sklearn.datasets import make_regression

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score

from sklearn.metrics import mean_absolute_error
from sklearn.metrics import median_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics.scorer import make_scorer

from sklearn.svm import LinearSVC
from sklearn.tree import DecisionTreeRegressor

svc = LinearSVC(random_state=0)
X, y = make_classification()

scoring = {'p': make_scorer(precision_score),
           'r': make_scorer(recall_score),
           'a': make_scorer(accuracy_score)}

# Multiple metric GridSearchCV
gs = GridSearchCV(svc,
                  param_grid={'C': np.linspace(0.1, 1)},
                  scoring=scoring, cv=5)
print("GridSearch instance:\n", gs.fit(X, y), "\n")
print("Best Estimators for all the scorers: \n", gs.best_estimator_, "\n")
print("Best scores for all the scorers: \n", gs.best_score_, "\n")
print("Best candidate index for all the scorers: \n", gs.best_index_, "\n")
print("The cv_results_ dict keys: \n", list(gs.cv_results_.keys()), "\n")


# Multiple metric cross_val_score
dtc = DecisionTreeRegressor()
X, y = make_regression(n_samples=100, random_state=42)
# For multiple metric - as list of metrics
print("Result of cross_val_score when multiple scorers are given as a list:")
print(cross_val_score(dtc, X, y, cv=2, scoring=['neg_mean_absolute_error',
                                                'neg_mean_squared_error',
                                                'neg_median_absolute_error']))
print("")

# For multiple metric - as dict of callables
neg_uae_scorer = make_scorer(mean_absolute_error)
neg_mse_scorer = make_scorer(mean_squared_error)
neg_mae_scorer = make_scorer(median_absolute_error)
print("Multiple scorer callables via dict:")
print(cross_val_score(dtc, X, y, cv=2,
                      scoring={'neg_mean_absolute_error': neg_uae_scorer,
                               'neg_mean_squared_error': neg_mse_scorer,
                               'neg_median_absolute_error': neg_mae_scorer}))
print("")

# For single metric
print("Single metric evaluation:")
print(cross_val_score(dtc, X, y, cv=2, scoring='neg_mean_absolute_error'))
