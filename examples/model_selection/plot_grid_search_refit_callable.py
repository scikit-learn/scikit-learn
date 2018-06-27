"""
=======================================================================
Balance model complexity and cross-validated score using refit=callable
=======================================================================

A simple example demonstrates the usage of `refit=callable` interface in
`sklearn.model_selection.GridSearchCV`. It shows this interface adds certain
amount of flexibility in identifying the "best" estimator. The function
passed to paramter `refit` incorporate of which metric(s) to optimise. This
interface can also be used in multiple metrics evaluation.

This example balances model complexity and cross-validated score by
finding a decent accuracy within 1 standard deviation of the best accuracy
score while minimising the number of PCA components.

The figure shows the trade of between cross-validated score and number of
PCA components. The balanced case is when n_components=6 and accuracy=0.77,
which falls into the range within 1 standard deviation of the best accuracy
score.
"""
# Author: Wenhao Zhang <wenhaoz@ucla.edu>

print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC


def score_bounds(scores):
    """
    Calculate the upper/lower bounds within 1 standard deviation
    of the best `mean_test_scores`.

    Args:
        scores: 'mean_test_score' here

    Return:
        upper/lower bounds within 1 standard deviation of the
        best `mean_test_scores`.
    """
    std_test_score = np.std(scores)
    best_test_score = max(scores)
    score_upper = np.minimum(best_test_score + std_test_score, 1)
    score_lower = np.maximum(best_test_score - std_test_score, 0)
    return score_upper, score_lower


def refit_callable(cv_results):
    """
    Balance model complexity with cross-validated score.

    Args:
        cv_results: see attribute cv_results_ of `GridSearchCV`.

    Return:
        Index of a model that has the fewest PCA components
        while has test score within 1 standard deviation of the best
        `mean_test_score`.
    """
    test_score_upper, test_score_lower = score_bounds(cv_results[
                                                      'mean_test_score'])
    n_components = cv_results['param_reduce_dim__n_components']
    test_scores = cv_results['mean_test_score']
    componet_score_lst = list(zip(n_components, test_scores))
    # Eliminate (n_comp, test_score) pairs that do not fall into
    # range "best_score +/- std"
    candidates = [x for x in componet_score_lst
                  if test_score_upper >= x[1] >= test_score_lower]
    res = min(candidates, key=lambda x: x[0])
    # Find best_index_ given fewest PCA components and decent score
    best_index = [x for x, y in enumerate(componet_score_lst)
                  if y[0] == res[0] and y[1] == res[1]][0]
    return best_index


pipe = Pipeline([
        ('reduce_dim', PCA(random_state=42)),
        ('classify', LinearSVC(random_state=42)),
])

N_FEATURES_OPTIONS = [2, 4, 6, 8]
param_grid = [
    {
        'reduce_dim__n_components': N_FEATURES_OPTIONS
    }
]

grid = GridSearchCV(pipe, cv=3, n_jobs=1, param_grid=param_grid,
                    scoring='accuracy', refit=refit_callable)
digits = load_digits()
grid.fit(digits.data, digits.target)

n_components = grid.cv_results_['param_reduce_dim__n_components']
test_scores = grid.cv_results_['mean_test_score']

plt.figure()
COLORS = 'bgr'
plt.bar(n_components, test_scores, width=1.3, color='bgrc')

upper, lower = score_bounds(grid.cv_results_['mean_test_score'])
plt.axhline(upper, linestyle='--', color='.5', label='+/- 1 std')
plt.axhline(np.max(test_scores), linestyle='--', color='y',
            label='Best score')
plt.axhline(lower, linestyle='--', color='.5')

plt.title("Balance model complexity and cross-validated score")
plt.xlabel('Reduced number of features')
plt.ylabel('Digit classification accuracy')
plt.xticks(N_FEATURES_OPTIONS)
plt.ylim((0, 1.2))
plt.legend(loc='upper left')

best_index_ = grid.best_index_

print("The best_index_ is %d" % best_index_)
print("The n_components selected is %d" % N_FEATURES_OPTIONS[best_index_])
print("The corresponding accuracy score is %.2f" %
      grid.cv_results_['mean_test_score'][best_index_])
plt.show()

