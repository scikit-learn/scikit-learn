# Authors: Christopher Erick Moody <chrisemoody@gmail.com>
# License: BSD 3 clause

from sklearn.datasets import fetch_olivetti_faces
from sklearn import manifold
from time import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn.decomposition import RandomizedPCA


faces = fetch_olivetti_faces()

verbose = 20

#----------------------------------------------------------------------
# Scale and visualize the embedding vectors
def plot_embedding(dataset, X, title=None, ax=None):
    x_min, x_max = np.min(X, 0), np.max(X, 0)
    X = (X - x_min) / (x_max - x_min)
    y = dataset.target
    if ax is None:
        ax = plt.subplot(111)
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

f, axs = plt.subplots(2, 2, figsize=(30, 30))

#----------------------------------------------------------------------
# Fit t-SNE on all 4096 dimensions with random initialization
tsne = manifold.TSNE(n_components=2, init='random', random_state=3,
                     n_iter=1000, verbose=verbose, early_exaggeration=50.0,
                     learning_rate=100, perplexity=50)
rpca = RandomizedPCA(n_components=50)
t0 = time()
reduced = rpca.fit_transform(faces.data)
embedded = tsne.fit_transform(reduced)
t1 = time()
dt = t1 - t0


t = ("Barnes-Hut t-SNE Visualization of Olivetti Faces digits\n" +
     "in %1.1f seconds & initialized embedding randomly")
t = t % dt
plot_embedding(faces, embedded, title=t, ax=axs[0, 0])


#-------------------------------------------------------------------------
# Fit t-SNE on all 4096 dimensions, but use PCA to initialize the embedding
tsne = manifold.TSNE(n_components=2, init='pca', random_state=3,
                     n_iter=1000, verbose=verbose, early_exaggeration=50.0,
                     learning_rate=100, perplexity=50)
t0 = time()
embedded = tsne.fit_transform(faces.data)
t1 = time()
dt = t1 - t0


t = ("Barnes-Hut t-SNE Visualization of Olivetti Faces digits\n" +
     "in %1.1f seconds & initialized embedding with PCA ")
t = t % dt
plot_embedding(faces, embedded, title=t, ax=axs[0, 1])


#----------------------------------------------------------------------
# Fit t-SNE on PCA-reduced dataset from 4096 dimensions to 50
tsne = manifold.TSNE(n_components=2, init='random', random_state=3,
                     n_iter=1000, verbose=verbose, early_exaggeration=50.0,
                     learning_rate=100, perplexity=50)
rpca = RandomizedPCA(n_components=50)
t0 = time()
reduced = rpca.fit_transform(faces.data)
embedded = tsne.fit_transform(reduced)
t1 = time()
dt = t1 - t0


t = ("Barnes-Hut t-SNE Visualization of Olivetti Faces digits\n" +
     "in %1.1f seconds & randomized PCA preprocessing")
t = t % dt
plot_embedding(faces, embedded, title=t, ax=axs[1, 0])


#----------------------------------------------------------------------
# Fit t-SNE using PCA and the exact method
tsne = manifold.TSNE(n_components=2, init='random', random_state=3,
                     method='exact', n_iter=10000, verbose=verbose,
                     learning_rate=100, perplexity=50)
t0 = time()
embedded = tsne.fit_transform(faces.data)
t1 = time()
dt = t1 - t0


t = ("Barnes-Hut t-SNE Visualization of Olivetti Faces digits\n" +
     "in %1.1f seconds with random initialization and exact method")
t = t % dt
plot_embedding(faces, embedded, title=t, ax=axs[1, 1])

plt.show()
