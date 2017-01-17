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
from matplotlib import pyplot as plt

from sklearn.datasets import make_hastie_10_2
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import precision_score
from sklearn.metrics.scorer import make_scorer
from sklearn.tree import DecisionTreeClassifier

print(__doc__)


X, y = make_hastie_10_2()

# The scorers can be a std scorer referenced by its name or one wrapped
# by sklearn.metrics.scorer.make_scorer
scoring = {'AUC Score': 'roc_auc', 'Precision': make_scorer(precision_score),
           'recall': 'recall', 'F1 Score': 'f1'}

# Multiple metric GridSearchCV, best_* attributes are exposed for the scorer
# with key 'AUC Score' ('roc_auc')
gs = GridSearchCV(DecisionTreeClassifier(random_state=42),
                  param_grid={'min_samples_split': range(2, 403, 10)},
                  scoring=scoring, cv=5, refit='AUC Score')
gs.fit(X, y)

results = gs.cv_results_

plt.figure().set_size_inches(15, 15)
plt.title("GridSearchCV evaluating using multiple scorers simultaneously",
          fontsize=16)

plt.xlabel("min_samples_split", fontsize=14)
plt.ylabel("Score", fontsize=14)
plt.grid()

ax = plt.axes()
ax.set_xlim(0, 402)

# Get the regular numpy array from the MaskedArray
X_axis = np.array(results['param_min_samples_split'].data, dtype=float)

for scorer, color in (('AUC Score', 'g'), ('F1 Score', 'grey'),
                      ('Precision', 'r'), ('recall', 'blue')):
    for sample, style in (('train', '--'), ('test', '-')):
        sample_score_mean = results['mean_%s_%s' % (sample, scorer)]
        sample_score_std = results['std_%s_%s' % (sample, scorer)]
        ax.fill_between(X_axis, sample_score_mean - sample_score_std,
                        sample_score_mean + sample_score_std,
                        alpha=0.1 if sample == 'test' else 0, color=color)
        ax.plot(X_axis, sample_score_mean, style, color=color,
                alpha=1 if sample == 'test' else 0.7,
                label="Mean %s on the %sing set" % (scorer, sample))

plt.legend(loc="best", fontsize=12)
plt.show()
