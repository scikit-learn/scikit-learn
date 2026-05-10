"""
=======================================================================
Plot the decision surface of decision trees trained on the iris dataset
=======================================================================

Plot the decision surface of a decision tree trained on pairs
of features of the iris dataset.

See :ref:`decision tree <tree>` for more information on the estimator.

For each pair of iris features, the decision tree learns decision
boundaries made of combinations of simple thresholding rules inferred from
the training samples.

We also show the tree structure of a model built on all of the features.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# First load the copy of the Iris dataset shipped with scikit-learn:
from sklearn.datasets import load_iris

iris = load_iris()


# %%
# Display the decision functions of trees trained on all pairs of features.
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.tree import DecisionTreeClassifier

for pairidx, pair in enumerate([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]):
    # We only take the two corresponding features
    X = iris.data[:, pair]
    y = iris.target

    # Train
    clf = DecisionTreeClassifier().fit(X, y)

    # Plot the decision boundary
    ax = plt.subplot(2, 3, pairidx + 1)
    plt.tight_layout(h_pad=0.5, w_pad=0.5, pad=2.5)
    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        response_method="predict",
        ax=ax,
        xlabel=iris.feature_names[pair[0]],
        ylabel=iris.feature_names[pair[1]],
        alpha=0.5,
    )

    # Plot the training points
    scatter = disp.ax_.scatter(
        X[:, 0],
        X[:, 1],
        c=y,
        cmap=ListedColormap(disp.multiclass_colors_),
        edgecolor="black",
        s=15,
    )

plt.suptitle("Decision surface of decision trees trained on pairs of features")
plt.figlegend(
    scatter.legend_elements()[0],
    iris.target_names,
    loc="lower center",
    ncols=len(iris.target_names),
)
_ = plt.axis("tight")

# %%
# Display the structure of a single decision tree trained on all the features
# together.
from sklearn.tree import plot_tree

plt.figure()
clf = DecisionTreeClassifier().fit(iris.data, iris.target)
plot_tree(clf, filled=True)
plt.title("Decision tree trained on all the iris features")
plt.show()
