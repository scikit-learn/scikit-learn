"""
===============================================================================
Decision boundary of semi-supervised classifiers versus SVM on the Iris dataset
===============================================================================

This example compares decision boundaries learned by two semi-supervised
methods, namely class:`~semi_supervised.LabelSpreading` and
class:`~semi_supervised.SelfTrainingClassifier`, when varying the proportion of
labeled training data from small fractions up to the full dataset.

Both methods rely on RBF kernels: Label Spreading uses it by default, and
Self-training is paired here with class:`~svm.SVC` as base estimator (also
RBF-based by default) to allow a fair comparison. Self-training with 100%
labeled data is omitted since it is identical to training a fully supervised SVC
directly.
"""

# %%
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import load_iris
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.semi_supervised import LabelSpreading, SelfTrainingClassifier
from sklearn.svm import SVC

iris = load_iris()
X = iris.data[:, :2]
y = iris.target

rng = np.random.RandomState(42)
y_rand = rng.rand(y.shape[0])
y_10 = np.copy(y)
y_10[y_rand > 0.1] = -1  # set random samples to be unlabeled
y_30 = np.copy(y)
y_30[y_rand > 0.3] = -1

ls10 = (LabelSpreading().fit(X, y_10), y_10, "LabelSpreading with 10% labeled data")
ls30 = (LabelSpreading().fit(X, y_30), y_30, "LabelSpreading with 30% labeled data")
ls100 = (LabelSpreading().fit(X, y), y, "LabelSpreading with 100% labeled data")

base_classifier = SVC(gamma=0.5, probability=True, random_state=42)
st10 = (
    SelfTrainingClassifier(base_classifier).fit(X, y_10),
    y_10,
    "Self-training with 10% labeled data",
)
st30 = (
    SelfTrainingClassifier(base_classifier).fit(X, y_30),
    y_30,
    "Self-training with 30% labeled data",
)
rbf_svc = (base_classifier.fit(X, y), y, "SVC with rbf kernel (100% labeled data)")

tab10 = plt.get_cmap("tab10")
color_map = {cls: tab10(cls) for cls in np.unique(y)}
color_map[-1] = (1, 1, 1)
classifiers = (ls10, st10, ls30, st30, ls100, rbf_svc)

fig, axes = plt.subplots(nrows=3, ncols=2, sharex="col", sharey="row", figsize=(10, 12))
axes = axes.ravel()

for ax, (clf, y_train, title) in zip(axes, classifiers):
    DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        response_method="predict_proba",
        plot_method="contourf",
        ax=ax,
    )
    colors = [color_map[label] for label in y_train]
    ax.scatter(X[:, 0], X[:, 1], c=colors, edgecolor="black")
    ax.set_title(title)

fig.suptitle("Unlabeled points are colored white", y=1.01)
plt.tight_layout()
plt.show()

# %%
# We observe that the decision boundaries are already quite similar to those
# using the full labeled data available for training, even when using a very
# small subset of the labels.
