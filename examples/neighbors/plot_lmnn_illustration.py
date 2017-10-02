"""
==========================================
Large Margin Nearest Neighbor Illustration
==========================================

Illustration of Large Margin Nearest Neighbor.
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.neighbors import LargeMarginNearestNeighbor
from sklearn.neighbors.lmnn import _select_target_neighbors


print(__doc__)

n_neighbors = 1
random_state = 0

# Create a tiny data set of 9 samples from 3 classes
X, y = make_classification(n_samples=9, n_features=2, n_informative=2,
                           n_redundant=0, n_classes=3, n_clusters_per_class=1,
                           class_sep=1.0, random_state=random_state)

# Spread out the data so that a margin of 1 is visible
X *= 6

# Find the target neighbors
target_neighbors = _select_target_neighbors(X, y, n_neighbors)

# Plot the points in the original space
plt.figure()
ax = plt.gca()

# Draw the graph nodes
ax.scatter(X[:, 0], X[:, 1], s=300, c=y, alpha=0.4)
for i in range(X.shape[0]):
    ax.text(X[i, 0], X[i, 1], str(i), va='center', ha='center')

# Annotate the reference sample
ref_pos = X[3]
ref_text_pos = ref_pos - np.array([1, 1.5])
ax.annotate('reference', xy=ref_pos, xytext=ref_text_pos, color='k',
            style='italic')

# Annotate the target neighbor relationship
tn_pos = X[target_neighbors[3]][0]
tn_text_pos = tn_pos + np.array([-4.5, 1])
ax.annotate('target neighbor', xy=tn_pos, xytext=tn_text_pos, color='k',
            style='italic')

# Draw a circle with the radius touching the target neighbor
tn_edge = tn_pos - ref_pos
inner_radius = np.sqrt(np.dot(tn_edge, tn_edge))
tn_circle = plt.Circle(ref_pos, inner_radius, color='b', linestyle='dashed',
                       fill=False, alpha=0.6)
ax.add_artist(tn_circle)

# Draw an outer circle indicating the margin
margin = 1.
margin_circle = plt.Circle(ref_pos, inner_radius + margin, color='b',
                           linestyle='dashed', fill=False, alpha=0.6)
ax.add_artist(margin_circle)

# Annotate the margin
margin_color = 'orange'


def fill_between_circles(ax, center, radii, color=margin_color):
    n = 50
    theta = np.linspace(0, 2*np.pi, n, endpoint=True)
    xs = np.outer(radii, np.cos(theta)) + center[0]
    ys = np.outer(radii, np.sin(theta)) + center[1]

    # in order to have a closed area, the circles
    # should be traversed in opposite directions
    xs[1, :] = xs[1, ::-1]
    ys[1, :] = ys[1, ::-1]

    ax.fill(np.ravel(xs), np.ravel(ys), facecolor=color, alpha=0.1)


theta_ref = -3*np.pi / 4
vec_ref = np.array([np.cos(theta_ref), np.sin(theta_ref)])
p_inner = ref_pos + vec_ref * inner_radius
p_outer = p_inner + vec_ref * margin
margin_text_pos = np.array([-13, -1])
middle = (p_inner + p_outer) / 2
ax.annotate('margin', xy=middle, xytext=margin_text_pos, style='italic',
            arrowprops=dict(facecolor='black', arrowstyle='->',
                            connectionstyle="arc3,rad=-0.3"))
fill_between_circles(ax, X[3], [inner_radius, inner_radius + margin])

# Annotate the impostors (1, 4, 5, 7)
imp_centroid = (X[1] + X[4] + X[5] + X[7]) / 4
imp_arrow_dict = dict(facecolor='black', color='gray', arrowstyle='->',
                      connectionstyle="arc3,rad=0.3")
ax.annotate('', xy=X[1], xytext=imp_centroid, arrowprops=imp_arrow_dict)
ax.annotate('', xy=X[4], xytext=imp_centroid, arrowprops=imp_arrow_dict)
ax.annotate('', xy=X[5], xytext=imp_centroid, arrowprops=imp_arrow_dict)
ax.annotate('', xy=X[7], xytext=imp_centroid, arrowprops=imp_arrow_dict)
ax.text(imp_centroid[0] - 1, imp_centroid[1] + 1, 'impostors', color='k',
        style='italic')

# Make axes equal so that boundaries are displayed correctly as circles
plt.axis('equal')
ax.set_title("Input space")


# Learn an embedding with LargeMarginNearestNeighbor
lmnn = LargeMarginNearestNeighbor(n_neighbors=n_neighbors, max_iter=30,
                                  random_state=random_state)
lmnn = lmnn.fit(X, y)

# Plot the points after transformation with LargeMarginNearestNeighbor
plt.figure()
ax2 = plt.gca()

# Get the embedding and find the new nearest neighbors
X_embedded = lmnn.transform(X)

ax2.scatter(X_embedded[:, 0], X_embedded[:, 1], s=300, c=y, alpha=0.4)
for i in range(len(X)):
    ax2.text(X_embedded[i, 0], X_embedded[i, 1], str(i),
             va='center', ha='center')

# Draw a circle with the radius touching the target neighbor
tn_edge = X_embedded[3] - X_embedded[target_neighbors[3]][0]
inner_radius = np.sqrt(np.dot(tn_edge, tn_edge))
tn_circle = plt.Circle(X_embedded[3], inner_radius, color='b',
                       linestyle='dashed', fill=False, alpha=0.6)
ax2.add_artist(tn_circle)

# Draw an outer circle indicating the margin
margin_circle = plt.Circle(X_embedded[3], inner_radius + margin, color='b',
                           linestyle='dashed', fill=False, alpha=0.6)
ax2.add_artist(margin_circle)

# Fill the margin with color
fill_between_circles(ax2, X_embedded[3], [inner_radius, inner_radius + margin])

# Make axes equal so that boundaries are displayed correctly as circles
plt.axis('equal')
ax2.set_title("LMNN embedding")

plt.show()
