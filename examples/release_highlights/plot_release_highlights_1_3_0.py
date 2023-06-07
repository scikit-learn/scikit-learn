# flake8: noqa
"""
=======================================
Release Highlights for scikit-learn 1.3
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.3! Many bug fixes
and improvements were added, as well as some new key features. We detail
below a few of the major features of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <changes_1_3>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# HDBSCAN: hierarchical density-based clustering
# ----------------------------------------------
# Originally hosted in the scikit-learn-contrib repository, :class:`cluster.HDBSCAN`
# has been adpoted into scikit-learn. It's missing a few features from the original
# implementation which will be added in future releases.
# By performing :class:`cluster.DBSCAN` over varying epsilon values
# :class:`cluster.HDBSCAN` finds clusters of varying densities making it more robust to
# parameter selection than :class:`cluster.DBSCAN`. More details in the
# :ref:`User Guide <hdbscan>`.
import numpy as np
from sklearn.cluster import HDBSCAN
from sklearn.datasets import load_digits
from sklearn.metrics import v_measure_score

X, true_labels = load_digits(return_X_y=True)
print(f"number of digits: {len(np.unique(true_labels))}")

hdbscan = HDBSCAN(min_cluster_size=15).fit(X)
non_noisy_labels = hdbscan.labels_[hdbscan.labels_ != -1]
print(f"number of clusters: {len(np.unique(non_noisy_labels))}")

v_measure_score(true_labels[hdbscan.labels_ != 1], non_noisy_labels)

# %%
# TargetEncoder: a new category encoding strategy
# -----------------------------------------------
# Well suited for categorical features with high cardinality,
# :class:`preprocessing.TargetEncoder` encodes the categories based on a shrunk
# estimate of the average target values for observations belonging to that category.
# More details in the :ref:`User Guide <target_encoder>`.
import numpy as np
from sklearn.preprocessing import TargetEncoder
X = np.array([["cat"] * 30 + ["dog"] * 20 + ["snake"] * 38], dtype=object).T
y = [90.3] * 30 + [20.4] * 20 + [21.2] * 38

enc = TargetEncoder(random_state=0)
X_trans = enc.fit_transform(X, y)

enc.encodings_

# %%
# Missing values support in decision trees
# ----------------------------------------
# The classes :class:`tree.DecisionTreeClassifier` and
# :class:`tree.DecisionTreeRegressor` now support missing values. For each potential
# threshold on the non-missing data, the splitter will evaluate the split with all the
# missing values going to the left node or the right node.
# More details in the :ref:`User Guide <tree_missing_value_support>`.
import numpy as np
from sklearn.tree import DecisionTreeClassifier

X = np.array([0, 1, 6, np.nan]).reshape(-1, 1)
y = [0, 0, 1, 1]

tree = DecisionTreeClassifier(random_state=0).fit(X, y)
tree.predict(X)
