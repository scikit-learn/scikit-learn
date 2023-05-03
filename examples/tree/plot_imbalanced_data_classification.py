"""
===================================================================
Decision Tree imbalanced data classification
===================================================================

In this example we demonstrate the affect of different split criteria
on decision tree classifier predictions on imbalanced data.

You can see that Gini index is biased towards the majority class while
Hellinger distance is sensitive to the ratio between the classes.

"""

# Import the necessary modules and libraries
import numpy as np
from sklearn import datasets
from sklearn.tree import DecisionTreeClassifier
import matplotlib.pyplot as plt

# Create imbalanced dataset
minority_class_ratio = 0.001
n_classes = 2
X, y = datasets.make_classification(
    n_samples=1000,
    n_features=2,
    n_informative=2,
    n_redundant=0,
    n_repeated=0,
    n_classes=n_classes,
    n_clusters_per_class=1,
    weights=[1 - minority_class_ratio],
    shuffle=False,
    random_state=0,
)

# Criteria to compare
criterions = ["gini", "entropy", "hellinger"]

# Plot parameters
fig = plt.figure(figsize=(10, 6))
max_subplots = len(criterions) * 2 - 1
plot_colors = ["darkgrey", "yellow"]
target_names = ["majority", "minority"]
markers = [".", "o"]

# Create mesh grid on feature space to draw classifier predictions
x_min, x_max = X[:, 0].min() - 0.2, X[:, 0].max() + 0.2
y_min, y_max = X[:, 1].min() - 2, X[:, 1].max() + 0.5
xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.1), np.arange(y_min, y_max, 0.01))

# Create plot per criterion
for sub_plot_idx, criterion in enumerate(criterions):
    # Add subplot to figure
    fig.add_subplot(1, max_subplots, sub_plot_idx * 2 + 1)

    # Train classifier
    clf = DecisionTreeClassifier(criterion=criterion)
    clf.fit(X, y)

    # Draw tree classifier prediction probability for minority class on feature space
    Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1].reshape(xx.shape)
    plt.contourf(xx, yy, Z, cmap=plt.cm.rainbow)
    plt.colorbar(label="minority class prediction probability")
    # Set the range for probability on the colorbar.
    plt.clim(0, 1)

    # Draw axis labels
    plt.xlabel("feature_1")
    plt.ylabel("feature_2")

    # Draw training points
    for i, color, marker in zip(range(n_classes), plot_colors, markers):
        idx = np.where(y == i)
        plt.scatter(
            X[:, 0][idx],
            X[:, 1][idx],
            c=color,
            label=target_names[i],
            cmap=plt.cm.RdYlBu,
            marker=marker,
            edgecolor="black",
        )

    # Draw subplot legend and title
    plt.legend(
        title="classes",
        handletextpad=0,
        loc="lower right",
        borderpad=0,
        scatterpoints=1,
        labelspacing=1,
    )
    plt.title(criterion)

    plt.plot()

plt.show()
