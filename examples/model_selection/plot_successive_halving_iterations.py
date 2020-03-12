"""
Successive Halving Iterations
=============================

This example illustrates how a successive halving search (
:class:`~sklearn.model_selection.HalvingGridSearchCV` and
:class:`~sklearn.model_selection.HalvingRandomSearchCV`) selectively chooses
the best parameter combination out of multiple candidates.

At the first iteration, a small amount of resources is used. The resource here
is the number of samples that the estimators are trained on. All candidates are
evaluated.

At the second iteration, only the best half of the candidates is evaluated.
The number of allocated resources is doubled: candidates are evaluated on
twice as many samples.

This process is repeated until the last iteration, where only 2 candidates
are left. The best candidate is the candidate that has the best score at the
last iteration.
"""
import pandas as pd
from sklearn import datasets
import matplotlib.pyplot as plt
from scipy.stats import randint
import numpy as np

from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.ensemble import RandomForestClassifier


print(__doc__)

rng = np.random.RandomState(0)

X, y = datasets.make_classification(n_samples=700, random_state=rng)

clf = RandomForestClassifier(n_estimators=20, random_state=rng)

param_dist = {"max_depth": [3, None],
              "max_features": randint(1, 11),
              "min_samples_split": randint(2, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

rsh = HalvingRandomSearchCV(
    estimator=clf,
    param_distributions=param_dist,
    resource='n_samples',
    max_resources='auto',  # max_resources=n_samples
    n_candidates='auto',  # choose n_cdts so that last iter exhausts resources
    cv=5,
    ratio=2,
    random_state=rng)
rsh.fit(X, y)

results = pd.DataFrame(rsh.cv_results_)
results['params_str'] = results.params.apply(str)
mean_scores = results.pivot(index='iter', columns='params_str',
                            values='mean_test_score')
ax = mean_scores.plot(legend=False, alpha=.6)

r_i_list = results.groupby('iter')['resource_iter'].unique()
labels = ['{}\nn_samples={}\nn_candidates={}'
          .format(i, r_i_list[i][0], rsh.n_candidates_[i])
          for i in range(rsh.n_iterations_)]
ax.set_xticklabels(labels)
ax.set_title('Scores of candidates over iterations')
ax.set_ylabel('mean test score', fontsize=15)
ax.set_xlabel('iterations', fontsize=15)
plt.show()
