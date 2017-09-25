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
X, y = make_classification(9, n_features=2, n_informative=2, n_redundant=0,
                           n_classes=3, n_clusters_per_class=1,
                           class_sep=1.0, random_state=random_state)

X *= 6  # Spread out the data so that a margin of 1 is visible

# Put the result into a color plot
fig = plt.figure(figsize=(15, 6))

# Find the target neighbors
tn = _select_target_neighbors(X, y, n_neighbors)

# Plot the points in the original space
ax = fig.add_subplot(121)

# Draw the graph nodes
ax.scatter(X[:, 0], X[:, 1], s=300, c=y, alpha=0.4)
for point_id in range(len(X)):
    ax.text(X[point_id, 0], X[point_id, 1], str(point_id), va='center',
            ha='center')

# Annotate the reference sample
ref_pos = X[3]
ref_text_pos = ref_pos - np.array([1, 1.5])
ax.annotate('reference', xy=ref_pos, xytext=ref_text_pos, color='k',
            style='italic')

# Annotate the target neighbor relationship
tn_pos = X[tn[3]][0]
middle = (ref_pos + tn_pos) / 2
tn_text_pos = tn_pos + np.array([-4.5, 1])
ax.annotate('target neighbor', xy=tn_pos, xytext=tn_text_pos, color='k',
            style='italic')

# Draw a circle with the radius touching the target neighbor
tn_edge = tn_pos - ref_pos
r = np.sqrt(np.dot(tn_edge, tn_edge))
tn_circle = plt.Circle(ref_pos, r, color='b', linestyle='dashed',
                       fill=False, alpha=0.6)
ax.add_artist(tn_circle)

# Draw an outer circle indicating the margin
margin = 1.
margin_circle = plt.Circle(ref_pos, r + margin, color='b', linestyle='dashed',
                           fill=False, alpha=0.6)
ax.add_artist(margin_circle)

# Annotate the margin
p1 = ref_pos + np.array([0, -r])
p2 = p1 + np.array([0, -margin])
plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color='k', linestyle='-', linewidth=2)
ax.annotate('margin', xy=(p1 + p2) / 2, xytext=(-5, -5), style='italic',
            arrowprops=dict(facecolor='black', arrowstyle='->',
                            connectionstyle="arc3,rad=-0.3"))

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
ax.set_title("Original space")


# Learn an embedding with LargeMarginNearestNeighbor
lmnn = LargeMarginNearestNeighbor(n_neighbors=n_neighbors, max_iter=30,
                                  random_state=random_state)
lmnn = lmnn.fit(X, y)

# Plot the points after transformation with LargeMarginNearestNeighbor
ax2 = fig.add_subplot(122)

# Get the embedding and find the new nearest neighbors
X_embedded = lmnn.transform(X)

ax2.scatter(X_embedded[:, 0], X_embedded[:, 1], s=300, c=y, alpha=0.4)
for point_id in range(len(X)):
    ax2.text(X_embedded[point_id, 0], X_embedded[point_id, 1],
             str(point_id), va='center', ha='center')

# Draw a circle with the radius touching the target neighbor
tn_edge = X_embedded[3] - X_embedded[tn[3]][0]
r = np.sqrt(np.dot(tn_edge, tn_edge))
tn_circle = plt.Circle(X_embedded[3], r, color='b', linestyle='dashed',
                       fill=False, alpha=0.6)
ax2.add_artist(tn_circle)

# Draw an outer circle indicating the margin
margin_circle = plt.Circle(X_embedded[3], r + margin, color='b',
                           linestyle='dashed', fill=False, alpha=0.6)
ax2.add_artist(margin_circle)

# Make axes equal so that boundaries are displayed correctly as circles
plt.axis('equal')

ax2.set_title("LMNN embedding")

plt.suptitle('Large Margin Nearest Neighbor', fontweight='bold', fontsize=16)
plt.show()
