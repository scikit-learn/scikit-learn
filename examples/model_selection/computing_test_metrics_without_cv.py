"""
==================================================================
Estimating The Area Under the Curve (AUC) Without Cross-Validation
==================================================================

Example of computing the area under the curve (AUC) metric to evaluate
classifier output quality using bootstrap sampling.

AUC measures the area under a receiver operating characteristic (ROC) curve,
which measures the true positive rate vs false positive rate. Typically,
the larger the AUC the better. One typically would like to estimate the AUC
metric for any held-out test dataset. An option is to compute this using
cross-validation, but if one has limited data, then one can still compute
it using bootstrap [1]. For more information on the general bootstrapping
procedure [2].

This example shows how to compute an estimate of the AUC metric of held-out
test data using bootstrap sampling. We first compute an overly optimistic
AUC by computing it on the entire dataset. We then leverage a custom
cross-validation function to do bootstrapping sampling to estimate
the optimism in the AUC metric. Finally, we adjust for the optimism
of the original AUC value and get an estimate of the AUC metric.

Note: This procedure can be generalized to any metric, and not just
necessarily AUC.

.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :func:`sklearn.model_selection.cross_validate`,
             :ref:`sphx_glr_auto_examples_model_selection_plot_roc.py`,

.. topic:: References:

   [1] Gordon C. S. Smith, Shaun R. Seaman, Angela M. Wood, Patrick Royston, Ian R. White, Correcting for Optimistic Prediction in Small Data Sets, American Journal of Epidemiology, Volume 180, Issue 3, 1 August 2014, Pages 318–324, https://doi.org/10.1093/aje/kwu140
   [2] B. Efron "Bootstrap Methods: Another Look at the Jackknife," The Annals of Statistics, Ann. Statist. 7(1), 1-26, (January, 1979)
"""
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn import svm, datasets
import numpy as np
print(__doc__)


# %%
# Data IO and generation
# ----------------------
#
# Import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target
X, y = X[y != 2], y[y != 2]
n_samples, n_features = X.shape
n_classes = 2

# Add noisy features to make the problem harder
random_state = np.random.RandomState()
n_samples, n_features = X.shape
X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]

# %%
# Classification and ROC analysis

# Learn to predict each class against the other
classifier = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True,
                                 random_state=random_state))
y_score = classifier.fit(X, y).decision_function(X)

# Get the original estimate of the AUC
auc_est = roc_auc_score(y, y_score)
print(f'Original over-optimistic AUC is {auc_est}')

# %%
# Estimating the held-out test dataset AUC using bootstrap

# Define the boostrap 'cv' class and object.
# Here, the train indices are bootstrap samples of the dataset,
# while the test indices are the original dataset.


class BootstrapSplit():
    """A custom cross validator.

    This class is used to sample bootstraps of the
    dataset, which are just overlapping copies of
    X and y. All custom cross validators need
    to implement the ``get_n_splits`` and ``splits``
    function.
    """

    def __init__(self, n_splits, random_state):
        self.n_splits = n_splits
        self.random_state = random_state

    def get_n_splits(self, X, y, groups=None):
        return self.n_splits

    def split(self, X, y, group):
        n = X.shape[0]
        i = 1
        while i <= self.n_splits:
            orig_idx = np.arange(n)
            idx = self.random_state.choice(n, size=n, replace=True)

            # make sure the bootstrap samples contain at least
            # a sample from each class
            y_train_unique = np.unique(y[orig_idx])
            y_test_unique = np.unique(y[idx])
            if len(y_train_unique) == 1 or len(y_test_unique) == 1:
                continue
            yield idx, orig_idx
            i += 1


# The number of bootstrap samples, which will then sample
# with overlap the dataset. This will then create 100 bootstrap
# copies of the dataset.
n_splits = 100
cv = BootstrapSplit(n_splits=n_splits, random_state=random_state)

# %%
# Run classifier with bootstrap
scores = cross_validate(
    estimator=classifier,
    X=X,
    y=y,
    scoring='roc_auc',
    cv=cv,
    n_jobs=-1,
    return_train_score=True
)

# Now calculate the optimism as the average of the AUC over
# the bootstrap samples.
optimism = np.mean(scores['train_score'] - scores['test_score'])
print(f'Optimism is {optimism}')

# Now calculate the optimism adjusted AUC metric.
adjusted_auc = auc_est - optimism
print(f'The AUC computed by bootstrapping and adjusting for optimism: {adjusted_auc}')

# %%
# Now compare this with the answer using cross-validation.
# The number of splits used here splits the data into non-overlapping sets.
# Hence, these are usually some low-number like 5, or 10.
cv = StratifiedKFold(n_splits=5)
cv_scores = cross_validate(
    estimator=classifier,
    X=X,
    y=y,
    scoring='roc_auc',
    cv=cv,
    n_jobs=-1,
)
cv_auc = np.mean(scores['test_score'])
print(f'The AUC computed by cross-validation: {cv_auc}')
