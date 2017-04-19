"""
=============================================================================
 t-SNE: The effect of various perplexity values on the shape
=============================================================================

An illustration of t-SNE on the two concentric circles and the S-curve
datasets for different perplexity values.

We observe the tendency of higher variations for small perplexity values and
clearer shapes as the value increases.

The size, the distance and the shape of clusters may vary upon initialization,
perplexity values and does not always convey a meaning.

As shown below, t-SNE for higher perplexities finds meaningful topology of
two concentric circles, however the size and the distance of the circles varies
slightly from the original. In contrary to two circles dataset, the shapes
visually diverge from S-curve topology on S-curve dateset even for larger
perplexity values.

For further details, "How to Use t-SNE Effectively"
http://distill.pub/2016/misread-tsne/ provides a good discussion of the
effects of various parameters, as well as interactive plots to explore
those effects.
"""

# License: BSD 3 clause

print(__doc__)

import matplotlib.pyplot as plt

from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D
from sklearn import manifold, datasets
from time import time

# Next line to silence pyflakes. This import is needed.
Axes3D

n_samples = 500
n_components = 2
subplot_index = 251
perplexities = [5, 50, 100, 150]

X, y = datasets.make_circles(n_samples=n_samples, factor=.5, noise=.05)

red = y == 0
green = y == 1

fig = plt.figure(figsize=(15, 8))
ax = fig.add_subplot(subplot_index)
ax.scatter(X[red, 0], X[red, 1], c="r")
ax.scatter(X[green, 0], X[green, 1], c="g")
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

for i, perplexity in enumerate(perplexities):
    subplot_index += 1
    ax = fig.add_subplot(subplot_index)

    t0 = time()
    tsne = manifold.TSNE(n_components=n_components, init='random',
                         random_state=0, perplexity=perplexity)
    Y = tsne.fit_transform(X)
    t1 = time()

    plt.title("%s %d (%.2g sec)" % ("Perplexity ", perplexity, t1 - t0))
    ax.scatter(Y[red, 0], Y[red, 1], c="r")
    ax.scatter(Y[green, 0], Y[green, 1], c="g")
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')

X, color = datasets.samples_generator.make_s_curve(n_samples, random_state=0)

subplot_index += 1
ax = fig.add_subplot(subplot_index, projection='3d')
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=plt.cm.Spectral)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.view_init(4, -72)

for i, perplexity in enumerate(perplexities):
    subplot_index += 1
    if subplot_index < 260:
        ax = fig.add_subplot(subplot_index)
    else:
        ax = fig.add_subplot(2, 5, 10)

    t0 = time()
    tsne = manifold.TSNE(n_components=n_components, init='random',
                         random_state=0, perplexity=perplexity)
    Y = tsne.fit_transform(X)
    t1 = time()

    plt.title("%s %d (%.2g sec)" % ("Perplexity ", perplexity, t1 - t0))
    plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=plt.cm.Spectral)
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')

plt.show()
