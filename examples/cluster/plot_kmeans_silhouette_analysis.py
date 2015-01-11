"""
===============================================================================
Silhouette analysis for sample data clustered using KMeans clustering algorithm
===============================================================================

Silhouette analysis can be used to study the separation distance between the
resulting clusters. The silhouette plot displays a measure of how close each
point in one cluster is to points in the neighboring clusters. This measure has
a range of [-1, 1]. Silhoette coefficients (as these values are referred to as)
near +1 indicate that the sample is far away from the neighboring clusters.
A value of 0 indicates that the sample is on or very close to the decision 
boundary between two neighboring clusters and negative values (upto -1) 
indicate that those samples might have been assigned to the wrong cluster.
"""

from __future__ import print_function

from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

print(__doc__)

# Generating the sample data from make_blobs
# This particular setting has one distict cluster and 3 clusters placed close
# together.
X, y = make_blobs(n_samples=500,
                  n_features=2,
                  centers=4,
                  cluster_std=1.0,
                  center_box=(-10.0, 10.0),
                  shuffle=True,
                  random_state=1)  # For reproducibility

range_n_clusters = [2, 4, 6]

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The n_clusters*10 are the additional samples to demarcate the space
    # between silhouette plots of individual clusters.
    ax1.set_ylim([0, len(X) + n_clusters * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    # This will hold the silhouette coefficient of all the clusters separated
    # by 0 samples.
    sorted_clustered_sample_silhouette_values = []

    for i in np.unique(cluster_labels):
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        # Add the ith_cluster_silhouette_values after sorting them
        ith_cluster_silhouette_values.sort()

        # The introduced 0 samples are to differentiate clearly between the
        # different clusters
        sorted_clustered_sample_silhouette_values += \
            ith_cluster_silhouette_values.tolist() + [0] * 10

    x_values = np.array(sorted_clustered_sample_silhouette_values)
    y_range = np.arange(len(X) + 10 * n_clusters)

    # Computing custom label coordinates for labeling the clusters
    # Plot the silhouette with the corresponding cluster color
    offset = 0
    for i in range(n_clusters):
        size_cluster_i = sum(cluster_labels == i)

        x = -0.05
        # Label them at the middle
        dy = size_cluster_i
        ax1.text(x, offset + 0.5 * dy, str(i))

        y_bottom = offset
        y_top = offset + dy

        color = cm.spectral(float(i) / n_clusters, 1)
        ax1.fill_betweenx(y_range, 0, x_values,
                          where=((y_range >= y_bottom) & (y_range < y_top)),
                          facecolor=color, edgecolor=color)
        # Compute the base offset for next plot
        offset += size_cluster_i + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")
    
    # The vertical line for average silhoutte score of all the values
    ax1.axvline(x = silhouette_avg, color = "red", linestyle = "--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    for k in range(len(X)):
        color = cm.spectral(float(cluster_labels[k]) / n_clusters, 1)
        ax2.scatter(X[k, 0], X[k, 1], marker='.', color=color)

    # Label the cluster centers with the cluster number for identification and
    # study of the corresponding silhouette plot.
    for c in clusterer.cluster_centers_:
        # Use the clusterer to know to which cluster number the current center
        # c belongs to
        i = clusterer.predict(c)[0]
        ax2.scatter(c[0], c[1], marker='o', c="white", alpha=1, s=200)
        ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50)

    ax2.set_title("The visualization of the clustered data.")
    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")
    plt.show()
