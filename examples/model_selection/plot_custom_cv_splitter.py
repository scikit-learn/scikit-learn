"""
=================================================
Implementing a custom cross-validation splitter
=================================================

This example shows how to implement a custom cross-validation splitter by
subclassing the base cross-validator API used throughout
:mod:`sklearn.model_selection`. Writing your own splitter is useful when the
built-in splitters (:class:`~sklearn.model_selection.KFold`,
:class:`~sklearn.model_selection.StratifiedKFold`, etc.) do not match the
structure of your data, for example when you need custom grouping or
ordering logic for the folds.

A splitter only needs to implement two methods:

- ``split(X, y=None, groups=None)``, a generator yielding ``(train_idx,
  test_idx)`` pairs of indices, and
- ``get_n_splits(X=None, y=None, groups=None)``, returning the number of
  splitting iterations.

"""

import numpy as np

from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

# %%
# A minimal custom splitter
# --------------------------
#
# Below, ``CustomSplitter`` divides the data into a fixed number of
# contiguous folds, similar in spirit to :class:`~sklearn.model_selection.KFold`
# but written from scratch to illustrate the required interface.


class CustomSplitter:
    def __init__(self, n_folds=5):
        self.n_folds = n_folds

    def split(self, X, y=None, groups=None):
        idxs = np.arange(X.shape[0])
        splits = np.array_split(idxs, self.get_n_splits())
        for split_idx, split in enumerate(splits):
            train_idxs = np.concatenate(
                [split for idx, split in enumerate(splits) if idx != split_idx]
            )
            test_idxs = split
            yield train_idxs, test_idxs

    def get_n_splits(self, X=None, y=None, groups=None):
        return self.n_folds


# %%
# Using the custom splitter with ``cross_val_score``
# -----------------------------------------------------
#
# Because ``CustomSplitter`` follows the same ``split`` / ``get_n_splits``
# interface as the built-in splitters, it can be passed directly as the
# ``cv`` argument to any tool that accepts a cross-validation splitter, such
# as :func:`~sklearn.model_selection.cross_val_score`.

X, y = datasets.load_iris(return_X_y=True)

clf = LogisticRegression(random_state=42)
scores = cross_val_score(clf, X, y, cv=CustomSplitter(n_folds=5))

print(f"Scores per fold: {scores}")
print(f"Mean score: {scores.mean():.3f}")
