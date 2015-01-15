"""
===================================
Comparison of t-SNE methods and options
===================================

A comparison of t-SNE applied to the Olivetti faces data. The 'exact' and
'barnes-hut' methods are visualized, as well as init='pca' and init='random'.
In addition to varying method & init, the script also compares preprocessing
the input data by reducing the dimensionality via PCA and also varying the
perplexity.

For a discussion and comparison of these algorithms, see the
:ref:`manifold module page <manifold>`

"""
# Authors: Christopher Erick Moody <chrisemoody@gmail.com>
# License: BSD 3 clause

from sklearn.datasets import fetch_olivetti_faces
from sklearn.decomposition import RandomizedPCA
from sklearn import manifold
from matplotlib import offsetbox
import matplotlib.pyplot as plt
import numpy as np
from time import time


faces = fetch_olivetti_faces()

verbose = 5

#----------------------------------------------------------------------
# Scale and visualize the embedding vectors
def plot_embedding(dataset, X, title=None, ax=None):
    x_min, x_max = np.min(X, 0), np.max(X, 0)
    X = (X - x_min) / (x_max - x_min)
    y = dataset.target
    if ax is None:
        plt.figure(figsize=(18, 15))
        ax = plt.gca()
    for i in range(X.shape[0]):
        ax.text(X[i, 0], X[i, 1], str(dataset.target[i]),
                color=plt.cm.Set1(y[i] / 10.),
                fontdict={'weight': 'bold', 'size': 30})

    if hasattr(offsetbox, 'AnnotationBbox'):
        # only print thumbnails with matplotlib > 1.0
        shown_images = np.array([[1., 1.]])  # just something big
        for i in range(dataset.data.shape[0]):
            shown_images = np.r_[shown_images, [X[i]]]
            r, g, b, a = plt.cm.Paired(i * 1.0 / dataset.data.shape[0])
            gray = dataset.images[i]
            rgb_image = np.zeros((64, 64, 3))
            rgb_image[:, :, 0] = gray * r
            rgb_image[:, :, 1] = gray * g
            rgb_image[:, :, 2] = gray * b
            imagebox = offsetbox.AnnotationBbox(
                offsetbox.OffsetImage(rgb_image),
                X[i], frameon=False)
            ax.add_artist(imagebox)
    ax.set_xticks([]), ax.set_yticks([])
    if title is not None:
        ax.set_title(title, loc='left')


#----------------------------------------------------------------------
# Fit t-SNE on 50 PCA-reduced dimensions with random initialization.
# We will frequently want to reduce the dimensionality when we have many
# dimensions since the t-SNE complexity grows linearly with the dimensionality.
tsne = manifold.TSNE(n_components=2, init='random', random_state=3,
                     n_iter=1000, verbose=verbose, early_exaggeration=50.0,
                     learning_rate=100, perplexity=50)
rpca = RandomizedPCA(n_components=50)
time_start = time()
reduced = rpca.fit_transform(faces.data)
embedded = tsne.fit_transform(reduced)
time_end = time()
dt = time_end - time_start


title = ("Barnes-Hut t-SNE Visualization of Olivetti Faces\n" +
         "in %1.1f seconds on 50-dimensional PCA-reduced data ")
title = title % dt
plot_embedding(faces, embedded, title=title)


#-------------------------------------------------------------------------
# Fit t-SNE on all 4096 dimensions, but use PCA to initialize the embedding
# Initializing the embedding with PCA allows the final emebedding to preserve
# global structure leaving t-SNE to optimize the local structure
tsne = manifold.TSNE(n_components=2, init='pca', random_state=3,
                     n_iter=1000, verbose=verbose, early_exaggeration=50.0,
                     learning_rate=100, perplexity=5)
time_start = time()
embedded = tsne.fit_transform(faces.data)
time_end = time()
dt = time_end - time_start


title = ("Barnes-Hut t-SNE Visualization of Olivetti Faces\n" +
         "in %1.1f seconds & initialized embedding with PCA ")
title = title % dt
plot_embedding(faces, embedded, title=title)


#----------------------------------------------------------------------
# Fit t-SNE on 50 PCA-reduced dimensions with random initialization, but reduce
# the perplexity. Note that the increased perplexity roughly increases the
# number of neighbors, which effectively decreases the number of clusters.
# For the example here, theres about 5 images per person, which with a
# perplexity of 50 forces multiple face-clusters to come together
tsne = manifold.TSNE(n_components=2, init='random', random_state=3,
                     n_iter=1000, verbose=verbose, early_exaggeration=50.0,
                     learning_rate=100, perplexity=50)
rpca = RandomizedPCA(n_components=50)
time_start = time()
reduced = rpca.fit_transform(faces.data)
embedded = tsne.fit_transform(reduced)
time_end = time()
dt = time_end - time_start


title = ("Barnes-Hut t-SNE Visualization of Olivetti Faces\n" +
         "in %1.1f seconds on 50-dimensional PCA-reduced data\n" +
         "with perplexity increased from 5 to 50")
title = title % dt
plot_embedding(faces, embedded, title=title)


#----------------------------------------------------------------------
# Fit t-SNE by random embedding and the exact method. The exact method
# is similar to Barnes-Hut, but this visualization reinforces the idea
# that the two methods yield similar results
tsne = manifold.TSNE(n_components=2, init='random', random_state=3,
                     method='exact', n_iter=10000, verbose=verbose,
                     learning_rate=100, perplexity=5)
time_start = time()
embedded = tsne.fit_transform(faces.data)
time_end = time()
dt = time_end - time_start


title = ("Exact t-SNE Visualization of Olivetti Faces\n" +
         "in %1.1f seconds with random initialization")
title = title % dt
plot_embedding(faces, embedded, title=title)

plt.show()
