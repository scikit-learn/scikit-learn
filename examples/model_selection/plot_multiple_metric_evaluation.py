"""Demonstration of multimetric evaluation on cross_val_score and GridSearchCV

Multiple metric grid search (or random search) can be done by setting the
``scoring`` parameter to a list of metric scorer names or a dict mapping the
scorer names to the scorer callables.

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
predefined scorers or a dict mapping the scorer names to the scorer callables.
"""

# Author: Raghav RV <rvraghav93@gmail.com>
# License: BSD

import numpy as np

from sklearn.datasets import make_classification
from sklearn.model_selection import GridSearchCV

from sklearn.metrics import precision_score
from sklearn.metrics.scorer import make_scorer

from sklearn.ensemble import RandomForestClassifier

from matplotlib import pyplot as plt

print(__doc__)


X, y = make_classification(n_samples=1000, random_state=42)

# The scorers can either be a std scorer referenced by its name or one wrapped
# by sklearn.metrics.scorer.make_scorer
scoring = {'Precision': make_scorer(precision_score), 'AUC Score': 'roc_auc'}

# Multiple metric GridSearchCV, best_* attributes are exposed for the scorer
# with key 'auc' ('roc_auc')
gs = GridSearchCV(RandomForestClassifier(random_state=42),
                  param_grid={'n_estimators': np.arange(1, 11)},
                  scoring=scoring, cv=3, refit='AUC Score')
gs.fit(X, y)

results = gs.cv_results_

plt.figure()
plt.title("GridSearchCV evaluating using multiple scorers simultaneously")

plt.xlabel("Values for number of estimators")
plt.ylabel("Score")
plt.grid()

ax = plt.axes()

# Get the regular numpy array from the MaskedArray
X_axis = np.array(results['param_n_estimators'].data, dtype=float)

for scorer, color in (('Precision', 'g'), ('AUC Score', 'r')):
    for sample, style in (('train', '--'), ('test', '-')):
        sample_score_mean = results['mean_%s_%s' % (sample, scorer)]
        sample_score_std = results['std_%s_%s' % (sample, scorer)]
        ax.fill_between(X_axis, sample_score_mean - sample_score_std,
                        sample_score_mean + sample_score_std,
                        alpha=0.1 if sample == 'test' else 0.05, color=color)
        ax.plot(X_axis, sample_score_mean, style, color=color,
                alpha=1 if sample == 'test' else 0.7,
                label="Mean %s on the %sing set" % (scorer, sample))

plt.legend(loc="lower right")
plt.show()
