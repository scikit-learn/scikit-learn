"""
=================================================================================
Nearest Neighbors Classification with different values of n_neighbors and weights
=================================================================================

Sample usage of Nearest Neighbors classification.
It will plot the decision boundaries for each class.

This example show the effect of tweaking the parameters ``n_neighbors`` and ``weights``. 
The smaller ``n_neighbors`` is the noisier the boundaries are.
"""
print(__doc__)

from time import time
from sklearn import neighbors, datasets
from matplotlib.colors import ListedColormap
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# import some data to play with
iris = datasets.load_iris()

# we only take the first two features. We could avoid this ugly
# slicing by using a two-dim dataset
X = iris.data[:, :2]
y = iris.target

# Create the mesh
h = .02  # step size in the mesh
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
mesh = np.c_[xx.ravel(), yy.ravel()]

# Create color maps
cmap_light = ListedColormap(['orange', 'cyan', 'cornflowerblue'])
cmap_bold = ['darkorange', 'c', 'darkblue']

tested_parameters = {
    "n_neighbors": [1, 5, 15, 150],
    "weights": ["uniform", "distance"]
}

fig = plt.figure(figsize=(15, len(tested_parameters["n_neighbors"])*5))
fig.subplots_adjust(wspace=0.1, hspace=0.3)

for i, n_neighbors in enumerate(tested_parameters["n_neighbors"]):
    for j, weights in enumerate(tested_parameters["weights"]):
        starting_time = time()
        clf = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
        clf.fit(X, y)
        training_time = time() - starting_time

        # Plot the decision boundary. For that, we will assign a color to each
        # point in the mesh [x_min, x_max]x[y_min, y_max].
        Z = clf.predict(mesh)
        testing_time = time() - starting_time - training_time

        # Put the result into a color plot
        Z = Z.reshape(xx.shape)
        ax = plt.subplot(len(tested_parameters["n_neighbors"]), 2, 2*i+j+1)
        plt.contourf(xx, yy, Z, cmap=cmap_light)

        # Plot also the training points
        sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=iris.target_names[y],
                        palette=cmap_bold, alpha=1.0, edgecolor="black")
        plt.xlim(xx.min(), xx.max())
        plt.ylim(yy.min(), yy.max())
        plt.title(f"k={n_neighbors}, weights='{weights}'")
        plt.text(.03, .9,
                 f"training_time: {training_time*1000 : .3f}ms", transform=ax.transAxes)
        plt.text(.03, .8,
                 f"testing_time: {testing_time*1000 : .3f}ms", transform=ax.transAxes)
        plt.xlabel(iris.feature_names[0])
        plt.ylabel(iris.feature_names[1])

plt.suptitle(f"3-Class classification on {len(y)} samples", y=.92, fontsize=18)
plt.show()
